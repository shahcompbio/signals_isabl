#' Read and Process Haplotype Files
#'
#' This function reads haplotype data from a list of file paths, applies optional filters, and formats the data for further analysis.
#'
#' @param paths A character vector of file paths to haplotype files. Each file is expected to contain columns such as `cell_id`, `chromosome`, `start`, `end`, `allele_id`, `hap_label`, and `readcount`.
#' @param cols_to_keep A character vector of column names to keep from the input data. Default columns include `cell_id`, `chr`, `start`, `end`, `allele_id`, `hap_label`, and `readcount`.
#' @param cell_ids A character vector of `cell_id`s to filter the data by. If `NULL`, no filtering by cell is applied. Default is `NULL`.
#' @param bin_ids A character vector of bin IDs to filter the data by. Each bin ID should follow the format `chr_start_end`. If `NULL`, no filtering by bins is applied. Default is `NULL`.
#' @param format_haplotypes Logical, if `TRUE` (default), the function converts the haplotype data into a wide format suitable for downstream analysis.
#'
#' @details 
#' If `bin_ids` is provided, the function extracts the bin size from the first `bin_id`, which is expected to follow the format `chr_start_end`. It uses this bin size to adjust the start and end positions of the data before filtering by bin IDs. 
#' If `format_haplotypes` is `TRUE`, the data is reshaped into wide format with separate columns for different alleles' read counts and a total count column is added.
#'
#' @return A `data.table` containing the processed haplotype data. The columns include `cell_id`, `chr`, `start`, `end`, `hap_label`, and read counts for each allele, among others.
#'
read_haplotypes_dlp <- function(paths,
                                cols_to_keep = c("cell_id", "chr", "start", "end", "allele_id", "hap_label", "readcount"),
                                cell_ids = NULL,
                                bin_ids = NULL,
                                format_haplotypes = TRUE){
  options("scipen" = 20)
  if (!is.null(bin_ids)){
    #get bin size
    hmmcopybinsize <- as.numeric(strsplit(bin_ids[1], "_")[[1]][3]) - as.numeric(strsplit(bin_ids[1], "_")[[1]][2]) + 1
  }
  
  hapslist <- list()
  for (x in paths){
    #read in haplotypes file
    haps_ <- data.table::fread(x, colClasses = list(character = c("chromosome", "cell_id"), 
                                                    integer = c("start", "end", "hap_label", "readcount", "allele_id")))
    if (!is.null(cell_ids)){
      #filter out cells not in HMMcopy data
      haps_ <- haps_[cell_id %in% cell_ids]
    }
    if (!is.null(bin_ids)){
      #filter out bins not in HMMcopy data
      haps_ <- haps_ %>% 
        .[, start := floor(start / hmmcopybinsize) * hmmcopybinsize + 1] %>%
        .[, end := start + hmmcopybinsize - 1] %>%
        .[, hbinid := paste(chromosome, start, end, sep = "_")] %>%
        .[hbinid %in% bin_ids] %>% 
        .[, hbinid := NULL]
    }

    if (dim(haps_)[1] == 0){
      message(paste0("No data in", x,  " post filtering, moving to next sample"))
      next
    }
    
    if (format_haplotypes) {
      #convert to wide format required by signals
      haps_ <- haps_ %>%
        .[, allele_id := paste0("allele", allele_id)] %>%
        data.table::dcast(., ... ~ allele_id, value.var = "readcount", fill = 0L, fun.aggregate = sum) %>%
        .[, lapply(.SD, sum), by = .(cell_id, chromosome, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
        .[, totalcounts := allele1 + allele0]
    }
    
    hapslist[[x]] <- haps_
    rm(haps_)
  }
  
  #merge all data files
  hapsdata <- data.table::rbindlist(hapslist, use.names = T)
  
  #change column name for consistency
  hapsdata <- dplyr::rename(hapsdata, chr = chromosome)
  
  rm(hapslist)
  
  return(hapsdata)
}

#' Read and Process Copy Number Data
#'
#' This function reads copy number data from a list of file paths, applies optional filters based on cell metrics, and returns a filtered dataset.
#'
#' @param cnpaths A character vector of file paths to copy number data files. Each file should contain columns such as `cell_id`, `chr`, `start`, `end`, `map`, `copy`, and `state`.
#' @param metricspaths A character vector of file paths to metrics data files. Metrics data is used for filtering cells based on quality and other criteria.
#' @param sample_ids A character vector of sample IDs to associate with the metrics files. If `NULL`, no sample IDs are assigned. Default is `NULL`.
#' @param filtercells Logical, if `TRUE` (default), cells are filtered based on quality and contamination criteria from the metrics data.
#' @param cols_to_keep A character vector of column names to retain from the copy number data. Default columns include `cell_id`, `chr`, `start`, `end`, `map`, `copy`, and `state`.
#' @param mappability A numeric threshold for filtering regions based on mappability. Only regions with mappability scores greater than this value are retained. If `NULL`, no mappability filter is applied. Default is `NULL`.
#' @param s_phase_filter Logical, if `TRUE`, cells in S-phase are filtered out. Default is `FALSE`.
#' @param quality_filter A numeric threshold for filtering cells based on quality. Only cells with quality scores higher than this value are retained. Default is `0.75`.
#' @param filter_reads A numeric threshold for filtering cells based on the total number of mapped reads. Default is `0`.
#'
#' @details 
#' The function reads metrics data and applies several filtering criteria, including quality, contamination, and mappability. Cells that pass the filters are used to subset the copy number data. The final output includes both the filtered copy number data and the associated cell metrics.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cn}{A `data.frame` containing the filtered copy number data.}
#'   \item{metrics}{A `data.frame` containing the filtered cell metrics.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- read_copynumber_dlp(cnpaths = c("cn_file1.csv", "cn_file2.csv"),
#'                               metricspaths = c("metrics_file1.csv", "metrics_file2.csv"),
#'                               sample_ids = c("sample1", "sample2"),
#'                               filtercells = TRUE)
#' }
#'
read_copynumber_dlp <- function(cnpaths,
                                metricspaths,
                                sample_ids = NULL,
                                filtercells = TRUE,
                                cols_to_keep = c("cell_id", "chr", "start", "end", "map", "copy", "state"),
                                mappability = 0.99,
                                s_phase_filter = FALSE,
                                quality_filter = 0.75,
                                filter_reads = 0){
  
  metricslist <- list()
  i <- 1
  for (x in metricspaths){
    if (!is.null(sample_ids)) {
      metricslist[[x]] <- data.table::fread(x)[, sample_id := sample_ids[i]]
      i <- i + 1
    } else {
      metricslist[[x]] <- data.table::fread(x)
      metricslist[[x]]$file <- x
    }
  }
  
  #merge metrics files
  metricsdata <- data.table::rbindlist(metricslist, use.names = T, fill = T)
  
  if (filtercells == TRUE){
    
    message(paste0("Total number of cells pre-filtering: ", dim(metricsdata)[1]))

    #find cells that pass filtering criteria
    cells_to_keep <- metricsdata %>%
      .[quality > quality_filter] %>%
      .[is_contaminated == FALSE] %>%
      .[total_mapped_reads > filter_reads] %>%
      .[is_control == FALSE]
    
    if (s_phase_filter == TRUE){
      cells_to_keep <- cells_to_keep[is_s_phase == FALSE]
    }
    
    metricsdata <- metricsdata[cell_id %in% cells_to_keep$cell_id]
    message(paste0("Total number of cells post-filtering: ", dim(metricsdata)[1]))
  }
  
  # read in copy number data, only keeping filtered cells
  cnlist <- list()
  for (x in cnpaths){
    cn_ <- data.table::fread(x,
                             colClasses = list(character = c("chr", "cell_id"),
                                               integer = c("start", "end", "reads", "state"),
                                               numeric = c("map", "gc", "copy")))
    cn_ <- cn_[cell_id %in% metricsdata$cell_id]
    cnlist[[x]] <- cn_
  }
  
  #bind copy number data
  cndata <- data.table::rbindlist(cnlist, use.names = T, fill = T)
  cndata <- cndata[, ..cols_to_keep]
  
  if (!is.null(mappability)) {
    cndata <- cndata[map > mappability]
  }
  
  rm(cnlist)
  
  return(list(cn = cndata %>% as.data.frame(), metrics = metricsdata %>% as.data.frame()))
}

#' Read and Process Copy Number Quality Control (QC) Metrics
#'
#' This function reads cell-level QC metrics from a list of file paths, applies optional filters based on quality and other criteria, and returns the filtered metrics data.
#'
#' @param metricspaths A character vector of file paths to QC metrics files. Each file should contain columns such as `cell_id`, `quality`, `is_contaminated`, `is_s_phase`, and `is_control`.
#' @param sample_ids A character vector of sample IDs to associate with the metrics files. If `NULL`, no sample IDs are assigned. Default is `NULL`.
#' @param filtercells Logical, if `TRUE` (default), cells are filtered based on quality, contamination, and S-phase status from the metrics data.
#' @param mappability A numeric threshold for filtering regions based on mappability. Not used in this function but kept for consistency. Default is `0.99`.
#' @param s_phase_filter Logical, if `TRUE`, cells in S-phase are filtered out. Default is `FALSE`.
#' @param quality_filter A numeric threshold for filtering cells based on quality. Only cells with quality scores higher than this value are retained. Default is `0.75`.
#' @param filter_reads A numeric threshold for filtering cells based on the total number of mapped reads. Not currently used in this function. Default is `0`.
#'
#' @details 
#' The function reads metrics data and applies several filtering criteria, including quality, contamination, and S-phase status. If `filtercells` is set to `TRUE`, the function filters out cells that do not meet these criteria.
#'
#' @return A `data.frame` containing the filtered metrics data.
#'
read_copynumber_dlpqc <- function(metricspaths,
                                sample_ids = NULL,
                                filtercells = TRUE,
                                mappability = 0.99,
                                s_phase_filter = FALSE,
                                quality_filter = 0.75,
                                filter_reads = 0.1e6){
  
  metricslist <- list()
  i <- 1
  for (x in metricspaths){
    if (!is.null(sample_ids)) {
      #read in metrics data
      metricslist[[x]] <- data.table::fread(x)[, sample_id := sample_ids[i]]
      i <- i + 1
    } else {
      metricslist[[x]] <- data.table::fread(x)
      metricslist[[x]]$file <- x
    }
  }
  
  metricsdata <- data.table::rbindlist(metricslist, use.names = T, fill = T)
  
  if (filtercells == TRUE){
    
    message(paste0("Total number of cells pre-filtering: ", dim(metricsdata)[1]))
    
    #find cells that pass filtering criteria
    cells_to_keep <- metricsdata %>%
      .[quality > quality_filter] %>%
      .[is_contaminated == FALSE] %>%
      .[is_s_phase == FALSE] %>%
      .[total_mapped_reads > filter_reads] %>%
      .[is_control == FALSE]
    
    if (s_phase_filter){
      cells_to_keep <- cells_to_keep[is_s_phase == FALSE]
    }
    
    metricsdata <- metricsdata[cell_id %in% cells_to_keep$cell_id]
    message(paste0("Total number of cells post-filtering: ", dim(metricsdata)[1]))
  }

  
  return(metrics = metricsdata %>% as.data.frame())
}

#' Filter Out Normal (Diploid) Cells
#'
#' This function filters out normal or diploid cells from copy number data, QC metrics, and haplotype data based on a specified threshold for non-diploid fractions and cell ploidy.
#'
#' @param cn A `data.table` containing copy number data. The table should have columns such as `cell_id` and `state`.
#' @param qc A `data.frame` or `data.table` containing QC metrics for cells. It should have at least a `cell_id` column.
#' @param haps A `data.table` containing haplotype data. It should also have a `cell_id` column for filtering purposes.
#' @param diploidcutoff A numeric value representing the cutoff for filtering cells based on the fraction of non-diploid states. Cells with a fraction of non-diploid states greater than this value are retained. Default is `0.05`.
#'
#' @details
#' The function computes the fraction of each cell that is diploid (state == 2) and filters cells based on the specified `diploidcutoff`. It calculates the fraction of the cell's genome in the most common state (ploidy), removes cells that are likely diploid or miscalled, and retains non-diploid cells for further analysis.
#'
#' If `diploidcutoff` is set to `0.0`, the function removes miscalled diploid cells where the most common state is greater than 2. Otherwise, it finds tumor cells by filtering out cells with a non-diploid fraction greater than the cutoff and a ploidy fraction below an outlier cutoff.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{cn}{A filtered `data.table` containing the copy number data for non-diploid cells.}
#'   \item{qc}{A `data.frame` containing QC metrics for non-diploid cells, with additional columns for ploidy cutoff and tumor cell identification.}
#'   \item{haps}{A filtered `data.table` containing haplotype data for non-diploid cells.}
#' }
#'
filter_normal_cells <- function(cn,
                               qc,
                               haps,
                               diploidcutoff = 0.05){
  
  message(paste0("Number of cells (pre filtering for normal cells): ", length(unique(cn$cell_id))))
  
  #find the fraction of each cell that is diploid and fraction that == cell ploidy
  cn_ploidy <- cn[, list(state_mean = mean(state), 
                         state_var = var(state), 
                         frac_nondiploid = sum(state != 2) / .N, 
                         pct_cell_ploidy = sum(state == signals:::Mode(state)) / .N, 
                         mode_state = signals:::Mode(state)), by = "cell_id"] %>% 
    .[order(pct_cell_ploidy, decreasing = TRUE)]
  
  qc <- qc %>% 
    left_join(cn_ploidy %>% dplyr::select(cell_id, frac_nondiploid, pct_cell_ploidy, state_mean), by = "cell_id")
  
  #find cells that are probably diploid or that are diploid + misscalled ploidy 
  cn_ploidy <- cn_ploidy[pct_cell_ploidy < (1 - diploidcutoff)]
  cutoff <- mean(cn_ploidy$pct_cell_ploidy) + 3 * sd(cn_ploidy$pct_cell_ploidy) #compute outlier cutoff
  message(paste0("Fraction ploidy cutoff: ", round(cutoff, 3)))
  
  if (diploidcutoff == 0.0){
    message("Removing misscalled diploid cells")
    cells_to_remove <- cn_ploidy[pct_cell_ploidy > 0.95 & mode_state > 2]
    non_diploidcells <- cn_ploidy[!(cell_id %in% cells_to_remove$cell_id)]
  } else{
    message("Finding tumour cells")
    non_diploidcells <- cn_ploidy[frac_nondiploid > diploidcutoff & pct_cell_ploidy < cutoff]
  }
  
  #filter out diploid cells
  cn <- cn[cell_id %in% non_diploidcells$cell_id]
  qc <- qc %>% 
    dplyr::mutate(pct_cell_ploidy_cutoff = cutoff, is_tumour_cell = cell_id %in% non_diploidcells$cell_id)
  haps <- haps[cell_id %in% unique(cn$cell_id)]
  
  message(paste0("Total number of cells after filtering diploid or misscalled ploidy: ", length(unique(cn$cell_id))))
  
  return(list(cn = cn, qc = qc, haps = haps))
}

#' Identify Tumor Cells Based on Copy Number Data
#'
#' This function identifies tumor cells by clustering copy number data and calculating the fraction of non-diploid states for each cell. Cells with a non-diploid fraction above a specified cutoff are classified as tumor cells.
#'
#' @param cn A `data.table` containing copy number data. The table should include columns such as `cell_id` and `state` to assess diploidy.
#' @param diploidcutoff A numeric value representing the cutoff for classifying a cell as a tumor cell based on its fraction of non-diploid states. Cells with a fraction of non-diploid states greater than this cutoff are classified as tumor cells. Default is `0.05`.
#'
#' @details
#' The function performs clustering on the copy number data using UMAP and DBSCAN, then computes the consensus copy number states for each cluster. For each cell, it calculates the fraction of the genome that is not in the diploid state (state != 2). Cells with a non-diploid fraction greater than the `diploidcutoff` are marked as tumor cells.
#'
#' @return A `data.table` with the following columns:
#' \describe{
#'   \item{cell_id}{The identifier for each cell.}
#'   \item{frac_nondiploid}{The fraction of the genome that is not in the diploid state for each cell.}
#'   \item{tumor_clone_id_clustering}{The clone ID derived from UMAP clustering.}
#'   \item{is_tumor_cell}{Logical indicating whether the cell is classified as a tumor cell (`TRUE`) or not (`FALSE`).}
#' }
#'
find_tumor_cells <- function(cn, diploidcutoff = 0.05){
  cl <- umap_clustering(cn, field = "copy", minPts = 5, umapmetric = "correlation")
  cn_clones <- consensuscopynumber(cn, cl$clustering)
  tumor_cell <- cn_clones %>% as.data.table(.) %>% 
    .[, list(frac_nondiploid = sum(state != 2) / .N), by = "cell_id"] %>% 
    .[, is_tumor_cell := ifelse(frac_nondiploid < diploidcutoff, FALSE, TRUE)] %>% 
    rename(clone_id = cell_id) %>% 
    full_join(cl$clustering) %>% 
    select(cell_id, frac_nondiploid, clone_id, is_tumor_cell) %>% 
    rename(tumor_clone_id_clustering = clone_id)
  return(tumor_cell)
}

#' Call Haplotype-Specific Copy Number (HSCN)
#'
#' This function calls haplotype-specific copy number (HSCN) states from copy number data, quality control (QC) metrics, and haplotype data. It optionally filters out normal (diploid) cells before inferring HSCN states.
#'
#' @param cn A `data.table` containing copy number data. The table should have columns such as `cell_id`, `state`, `chr`, `start`, `end`.
#' @param qc A `data.table` containing QC metrics for the cells, including `cell_id`.
#' @param haps A `data.table` containing haplotype data, which also includes `cell_id`.
#' @param diploidcutoff A numeric value representing the cutoff for classifying a cell as diploid based on its fraction of non-diploid states. Default is `0.05`.
#' @param mincells The minimum number of cells required to infer HSCN states. Default is `8`.
#' @param filternormalcells Logical, if `TRUE`, the function filters out normal (diploid) cells based on the `diploidcutoff`. Default is `FALSE`.
#' @param maskedbins A character vector specifying bins to mask (exclude) from the analysis. Default is `NULL`.
#' @param global_phasing_for_balanced Logical, if `TRUE`, global phasing (all cells) is applied for balanced chromosomes. Default is `FALSE`.
#' @param chrs_for_global_phasing A character vector of chromosomes for which global phasing should be applied. Default is `NULL`.
#' @param frachaps A numeric value specifying the fraction of haplotypes to consider. Default is `0.8`.
#' @param ncores The number of cores to use for parallel processing. Default is `1`.
#' @param female Logical, indicating whether the sample is female, which affects chromosome X analyses. Default is `TRUE`.
#'
#' @details
#' The function first checks if normal (diploid) cells should be filtered using the `find_tumor_cells` function. 
#' If normal cells are filtered, the copy number, QC, and haplotype data are restricted to tumor cells.
#' The function then calls haplotype-specific copy number states using the `callHaplotypeSpecificCN` function. 
#' If fewer than 10 cells remain after filtering, the function returns an error. If valid data exist, the HSCN results are returned, including updated QC information.
#'
#' @return A list containing HSCN data, including the inferred copy number states, and the updated QC metrics (`qc_per_cell`).
#'
callhscn <- function(cn,
                     qc,
                     haps,
                     diploidcutoff = 0.05,
                     mincells = 8,
                     filternormalcells = FALSE,
                     maskedbins = NULL,
                     global_phasing_for_balanced = FALSE,
                     chrs_for_global_phasing = NULL,
                     frachaps = 0.8,
                     selftransitionprob = 0.95,
                     chr_cell_list = NULL,
                     ncores = 1,
                     female = TRUE){

  cn <- as.data.table(cn)
  qc <- as.data.table(qc)
  haps <- as.data.table(haps)
  
  if (filternormalcells == FALSE){
    message(paste0("Number of cells: ", length(unique(cn$cell_id))))
  } else if (filternormalcells == TRUE){
    message(paste0("Number of cells (pre filtering for normal cells): ", length(unique(cn$cell_id))))
    x <- find_tumor_cells(cn, diploidcutoff)
    cells_to_keep <- x[is_tumor_cell == TRUE]$cell_id
    cn <- cn[cell_id %in% cells_to_keep]
    qc <- qc[cell_id %in% cells_to_keep]
    haps <- haps[cell_id %in% cells_to_keep]
    message(paste0("Total number of cells after filtering diploid or misscalled ploidy: ", length(unique(cn$cell_id))))
  }
  
  if (length(unique(cn$cell_id)) < 10){
    stop("Number of cells is < 10, unable to run signals")
    hscn <- list(data = NULL, message = "Number of cells is less than 10, unable to run signals")
    return(hscn)
  }

  if (length(unique(cn$cell_id)) == 0){
    stop("No cells returning NULL")
    return(NULL)
  }

  message("Infer HSCN")
  hscn <- callHaplotypeSpecificCN(cn,
                                  haps,
                                  maskedbins = maskedbins,
                                  chrs_for_global_phasing = chrs_for_global_phasing,
                                  global_phasing_for_balanced = global_phasing_for_balanced,
                                  selftransitionprob = selftransitionprob,
                                  likelihood = "auto",
                                  mincells = mincells,
                                  cluster_per_chr = TRUE, 
                                  chr_cell_list = chr_cell_list,
                                  ncores = ncores,
                                  female = female)
  
  #add library and sample IDs to qc
  qc$sample_id <- unlist(signals:::get_library_labels(qc$cell_id, idx=1))
  qc$library_id <- unlist(signals:::get_library_labels(qc$cell_id, idx=2))

  #add qc to signals object
  hscn$qc_per_cell <- dplyr::left_join(hscn$qc_per_cell, qc, by = "cell_id")
  
  return(hscn)
}

callascn <- function(hscn,
                     qc,
                     haps,
                     diploidcutoff = 0.05,
                     ncores = 1){
  
  message("Infer ASCN")
  ascn <- callAlleleSpecificCNfromHSCN(hscn, ncores = ncores)
  ascn <- filtercn(ascn)
  ascn$qc_per_cell <- dplyr::left_join(ascn$qc_per_cell, hscn$qc_per_cell)
  
  return(ascn)
}

extra_qc_annotations <- function(res){
  fracLOH <- res$data %>% 
    as.data.table() %>% 
    .[, list(fracLOH = sum(LOH == "LOH") / .N), by = c("chr", "cell_id")] %>% 
    .[chr != "Y"] %>% 
    pivot_wider(names_from = "chr", values_from = "fracLOH", names_prefix = "fracLOH_")
  
  cn_ploidy <- as.data.table(res$data)[, list(state_mean = mean(state),
                         state_var = var(state), 
                         frac_nondiploid = sum(state != 2) / .N, 
                         pct_cell_ploidy = sum(state == signals:::Mode(state)) / .N), by = "cell_id"] %>% 
    .[order(pct_cell_ploidy, decreasing = TRUE)]
  
  res$qc_per_cell <- left_join(res$qc_per_cell, fracLOH, by = "cell_id") %>% 
    left_join(cn_ploidy, by = "cell_id")
  return(res)
}


#' Run the SIGNALS Pipeline for HSCN and QC Analysis
#'
#' This function executes the SIGNALS pipeline for haplotype-specific copy number (HSCN) analysis and quality control (QC) on copy number, QC, and haplotype data. It generates visualizations such as heatmaps and QC plots.
#'
#' @param args A list of arguments containing paths to the required data and additional parameters:
#' - `hmmcopyreads`: Path(s) to the copy number data files.
#' - `hmmcopyqc`: Path(s) to the QC metrics files.
#' - `allelecounts`: Path to the allele counts (haplotype) data.
#' - `cell_list`: Optional, path to a file containing a list of cell IDs to retain.
#' - `maskbins`: Optional, path to a file with masked bins (excluded from the analysis).
#' - `chrs_for_global_phasing`: Optional, a comma-separated string of chromosomes for global phasing.
#' - `sex`: String indicating the sex of the sample ("MALE" or "FEMALE").
#' - `mincells`: Minimum number of cells required for HSCN analysis.
#' - `diploidcutoff`: Cutoff for defining diploid cells based on non-diploid fractions.
#' - `ncores`: Number of cores to use for parallel processing.
#' - `filternormalcells`: Logical, if `TRUE`, filters normal (diploid) cells before HSCN analysis.
#' - `qcplot`: File path to save the generated QC plots.
#' - `csvfile`: File path to save the resulting HSCN data as a CSV file.
#' - `qccsvfile`: File path to save the QC data as a CSV file.
#' - `Rdatafile`: File path to save the full results as an RDS file.
#' - `heatmap`: File path to save the HSCN heatmap.
#' - `heatmapraw`: File path to save the raw copy and BAF heatmaps.
#' - `maxcellsplotting`: Maximum number of cells to plot in heatmaps.
#' - `mappability`: Optional numeric threshold for filtering bins by mappability (only bins with map > this value are kept). Default `NULL` (no filter).
#'
#' @details
#' The function begins by reading in copy number data, QC metrics, and haplotype data, optionally filtering the cell list and applying masking of bins. 
#' It then calls haplotype-specific copy number (HSCN) using the `callhscn` function, applies QC filters, and generates various QC plots.
#' The results, including inferred HSCN states and quality metrics, are saved in CSV and RDS formats. 
#' Heatmaps showing cell clustering and HSCN states are generated and saved.
#'
#' @return
#' The function does not return any values but saves the results and plots to the specified output files.
#'
run_signals <- function(args){
  
  cndata <- read_copynumber_dlp(cnpaths = args$hmmcopyreads,
                                metricspaths = args$hmmcopyqc,
                                mappability = args$mappability)
  cell_ids <- unique(cndata$metrics$cell_id)
  message(paste0("Number of cells: ", length(cell_ids)))
  message(paste0("Number of bins: ", dim(distinct(cndata$cn, chr, start))[1]))
  
  if (!is.null(args$cell_list)){
    message(paste0("Filtering for cells in ", args$cell_list))
    cell_list <- data.table::fread(args$cell_list, header = FALSE)
    cndata$cn <- filter(cndata$cn, cell_id %in% cell_list$V1)
    cndata$metrics <- filter(cndata$metrics, cell_id %in% cell_list$V1)
    message(paste0("Number of cells: ", length(cell_ids)))
  }
  
  bin_ids <- as.data.table(cndata$cn)[, c("chr", "start", "end")] %>%
    unique(., by = c("chr", "start", "end")) %>%
    .[, binid := paste(chr, start, end, sep = "_")] %>%
    .$binid
  
  message("Read in haplotypes data")
  haplotypes <- read_haplotypes_dlp(args$allelecounts, 
                                    cell_ids = cell_ids, 
                                    bin_ids = bin_ids, 
                                    format_haplotypes = TRUE)
  
  message("Read in maskbins")
  if (!is.null(args$maskbins)){
    maskedbins <- fread(args$maskbins)
  } else{
    maskedbins <- NULL
  }
  
  message("Set up parameters for global phasing")
  if (!is.null(args$chrs_for_global_phasing)){
    chrs_for_global_phasing <- strsplit(args$chrs_for_global_phasing, ",")[[1]]
    global_phasing_for_balanced <- TRUE
  } else{
    chrs_for_global_phasing <- NULL
    global_phasing_for_balanced <- FALSE
  }
  
  if (args$sex == "FEMALE"){
    female <- TRUE
  } else{
    female <- FALSE
  }
  
  message("Process chr_cell_list if provided")
  chr_cell_list_processed <- NULL
  if (!is.null(args$chr_cell_list)){
    message(paste0("Reading chr_cell_list from: ", args$chr_cell_list))
    chr_cell_df <- data.table::fread(args$chr_cell_list)
    
    # Check if required columns exist
    if (!all(c("chr", "cell_id") %in% colnames(chr_cell_df))){
      stop("chr_cell_list file must contain columns 'chr' and 'cell_id'")
    }
    
    # Get unique chromosomes from CN data
    cn_chromosomes <- unique(cndata$cn$chr)
    chr_cell_chromosomes <- unique(chr_cell_df$chr)
    
    # Check if chromosomes in chr_cell_list match those in CN data
    missing_in_cn <- setdiff(chr_cell_chromosomes, cn_chromosomes)
    if (length(missing_in_cn) > 0){
      stop(paste0("Chromosomes in chr_cell_list not found in CN data: ", paste(missing_in_cn, collapse = ", ")))
    }
    
    # Convert to named list format required by callHaplotypeSpecificCN
    chr_cell_list_processed <- split(chr_cell_df$cell_id, chr_cell_df$chr)
    # Convert to character names as expected by signals
    names(chr_cell_list_processed) <- as.character(names(chr_cell_list_processed))
    
    message(paste0("Processed chr_cell_list for ", length(chr_cell_list_processed), " chromosomes"))
    for (chr in names(chr_cell_list_processed)){
      message(paste0("  Chr ", chr, ": ", length(chr_cell_list_processed[[chr]]), " cells"))
    }
  }
  
  message("Call haplotype specific copy number")
  res <- callhscn(cndata$cn,
                  cndata$metrics,
                  haplotypes,
                  frachaps = 0.8,
                  maskedbins = maskedbins,
                  mincells = args$mincells,
                  diploidcutoff = args$diploidcutoff,
                  ncores = args$ncores,
                  filternormalcells = args$filternormalcells,
                  chrs_for_global_phasing = chrs_for_global_phasing,
                  global_phasing_for_balanced = global_phasing_for_balanced,
                  selftransitionprob = args$selftransitionprob,
                  chr_cell_list = chr_cell_list_processed,
                  female = female)
  res <- extra_qc_annotations(res)
  
  print(res)
  
  # message("Call allele specific copy number")
  # res2 <- callascn(res, ncores = args$ncores)
  # res2$likelihood <- res$likelihood
  # print(res2)
  
  message("Make QC plots")
  g1 <- plotBAFperstate(res$data %>% filter(totalcounts > 9))+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  g2 <- plot_variance_state(res$data %>% filter(totalcounts > 9))+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  g3 <- res$qc_per_cell %>% 
    mutate(Ploidy = paste0(ploidy)) %>% 
    ggplot(aes(x = average_distance, fill = Ploidy)) + 
    geom_histogram(alpha = 0.7, bins = 50) +
    xlab("Average distance raw BAF\nto expected BAF per cell") +
    ylab("Counts") +
    theme(legend.position = "none") +
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  g4 <- res$qc_per_cell %>% 
    mutate(Ploidy = paste0(ploidy)) %>% 
    ggplot(aes(x = totalhapcounts / 1e6, fill = Ploidy)) + 
    geom_histogram(alpha = 0.7, bins = 50) +
    xlab("Total number of hap counts per cell (Millions)") +
    ylab("Counts") +
    theme(legend.position = "none") +
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  g5 <- res$qc_per_cell %>% 
    mutate(Ploidy = paste0(ploidy)) %>% 
    ggplot(aes(x = pct_cell_ploidy, fill = Ploidy)) + 
    geom_histogram( alpha = 0.7, bins = 50) +
    xlab("Fraction of genome = ploidy per cell") +
    ylab("Counts") +
    theme(legend.position = c(0.1, 0.9)) +
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  g6 <- plot_clusters_used_for_phasing(res) +
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    )
  
  gall <- cowplot::plot_grid(cowplot::plot_grid(g1, g2, align = "h", axis = "tb"), cowplot::plot_grid(g3, g4, g5, ncol = 3), g6, ncol = 1, rel_heights = c(0.25, 0.25, 0.5))
  
  cowplot::save_plot(args$qcplot,
                     gall,
                     base_width = 15, base_height = 20)
  
  #filter out cells where BAF and inferred state indicate possible errors
  cell_ids <- unique(res$data$cell_id)
  message(paste0("Number of cells (pre final QC filtering): ", length(cell_ids)))
  res <- filtercn(res)
  cell_ids <- unique(res$data$cell_id)
  message(paste0("Number of cells (after final QC filtering): ", length(cell_ids)))
  
  message("Cluster data")
  cl <- umap_clustering(res$data, 
                        field = "copy", 
                        hscn = FALSE,
                        umapmetric = "correlation",
                        minPts = round(0.03 * length(cell_ids)))
  
  res$qc_per_cell <- dplyr::left_join(res$qc_per_cell, cl$clustering %>% dplyr::select(cell_id, clone_id), by = "cell_id")
  
  message("Write csv file")
  fwrite(x = res$data, file = args$csvfile)
  fwrite(x = res$qc_per_cell, file = args$qccsvfile)
  
  message("Write Rdata file")
  saveRDS(file = args$Rdatafile, object = list(hscn = res, cl = cl))
  
  message("Make heatmaps")
  if (dim(cl$clustering)[1] > args$maxcellsplotting){
    message(paste0(dim(cl$clustering)[1], " total cells, downsampling to ", args$maxcellsplotting, " for plotting"))
    cells <- dplyr::sample_n(cl$clustering, args$maxcellsplotting) %>% pull(cell_id)
  } else{
    cells <- cl$clustering %>% pull(cell_id)
  }
  
  if (length(unique(cl$clustering$clone_id)) > 20){
    show_clone_label <- FALSE
  } else{
    show_clone_label <- TRUE
  }
  
  
  h1 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                    tree = ape::keep.tip(cl$tree, cells),
                    reorderclusters = T,
                    show_clone_label = show_clone_label,
                    clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                    plottree = FALSE,
                    plotcol = "state",
                    str_to_remove = "SPECTRUM-OV-")
  h2 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                    tree = ape::keep.tip(cl$tree, cells),
                    reorderclusters = T,
                    clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                    plottree = FALSE,
                    plotcol = "state_phase", 
                    show_clone_label = F, 
                    show_library_label = F,
                    str_to_remove = "SPECTRUM-OV-")
  
  png(args$heatmap, width = 25, height = 12, units = "in", res = 300)
  print(ComplexHeatmap::draw(h1 + h2,
                              ht_gap = unit(0.6, "cm"),
                              column_title = paste0("Total number of cells: ", length(unique(res$data$cell_id)), "\nNumber of cells plotted: ", length(cells)),
                              column_title_gp = grid::gpar(fontsize = 20),
                              heatmap_legend_side = "bottom",
                              annotation_legend_side = "bottom",
                              show_heatmap_legend = TRUE))
  dev.off()
  
  h3 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                    tree = ape::keep.tip(cl$tree, cells),
                    reorderclusters = T,
                    show_clone_label = show_clone_label,
                    clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                    plottree = FALSE,
                    plotcol = "copy",
                    str_to_remove = "SPECTRUM-OV-")
  
  h4 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                    tree = ape::keep.tip(cl$tree, cells),
                    reorderclusters = T,
                    clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                    plottree = FALSE,
                    plotcol = "BAF",
                    show_clone_label = F,
                    show_library_label = F,
                    str_to_remove = "SPECTRUM-OV-")
  
  png(args$heatmapraw, width = 25, height = 12, units = "in", res = 300)
  print(ComplexHeatmap::draw(h3 + h4,
                              ht_gap = unit(0.6, "cm"),
                              column_title = paste0("Total number of cells: ", length(unique(res$data$cell_id)), "\nNumber of cells plotted: ", length(cells)),
                              column_title_gp = grid::gpar(fontsize = 20),
                              heatmap_legend_side = "bottom",
                              annotation_legend_side = "bottom",
                              show_heatmap_legend = TRUE))
  dev.off()
  
  message("Finished")
}


main <- function(){
  
  message(paste0("signals version: ", packageVersion("signals")))
  
  parser <- ArgumentParser()
  
  parser$add_argument("--hmmcopyreads", default="character", nargs = "+",
                      help = "hmmcopy reads files")
  parser$add_argument("--filternormalcells", action="store_true", default=FALSE,
                      help="Optionally filter out normal cells")
  parser$add_argument("--hmmcopyqc", default=NULL, type="character", nargs = "+",
                      help="hmmcopy QC files")
  parser$add_argument("--allelecounts", default=NULL, type="character", nargs = "+",
                      help="count haps files")
  parser$add_argument("--ncores", default=NULL, type="integer",
                      help="Number of cores in signals inference")
  parser$add_argument("--csvfile", default=NULL, type="character",
                      help="output csvfile")
  parser$add_argument("--qccsvfile", default=NULL, type="character",
                      help="output csvfile")
  parser$add_argument("--Rdatafile", default=NULL, type="character",
                      help="output Rdata file")
  parser$add_argument("--heatmap", default=NULL, type="character",
                      help="heatmap plot")
  parser$add_argument("--heatmapraw", default=NULL, type="character",
                      help="heatmap plot of raw data")
  parser$add_argument("--qcplot", default=NULL, type="character",
                      help="QC plot")
  parser$add_argument("--diploidcutoff", default=0.05, type="double",
                      help="Cutoff for filtering diploid cells, this is the fraction of the genome that is non diploid")
  parser$add_argument("--maxcellsplotting", default=2500, type="integer",
                      help="Max number of cells to plot in the heatmap")
  parser$add_argument("--mincells", default=8, type="integer",
                      help="Min number of cells for per chromosome clustering")
  parser$add_argument("--selftransitionprob", default=0.95, type="double",
                      help="Self transition probability for HMM")
  parser$add_argument("--maskbins", default=NULL, type="character",
                      help="HMMcopy style table with bins to be masked")
  parser$add_argument("--chr_cell_list", default=NULL, type="character",
                      help="List cells to use for phasing for each chromosome, needs to be the path to a csv with two columns, chr and cell_id")
  parser$add_argument("--chrs_for_global_phasing", default=NULL, type="character",
                      help="Chromosomes to phase using global phasing for diploid regions")
  parser$add_argument("--sex", default="FEMALE", type="character",
                      help="Patient sex, either FEMALE or MALE, this controls what happens with chromosome X")
  parser$add_argument("--cell_list", default=NULL, type="character",
                      help="List of cells in a text file, one cell_id per row, no header. Only use this subset of cells.")
  parser$add_argument("--mappability", default=NULL, type="double",
                      help="Mappability threshold; only bins with map > this value are kept. Default NULL (no filter).")
  args <- parser$parse_args()
  
  print(args)
  #saveRDS(args, file = "args.rdata")
  
  run_signals(args)
}

library(argparse)
library(dplyr)
library(tidyr)
library(data.table)
library(signals)
#devtools::load_all("/home/william1/bin/R/signals")
library(ggplot2)
library(cowplot)

main()
print(sessionInfo())
