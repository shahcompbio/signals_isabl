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
find_tumor_cells <- function(cn, diploidcutoff = 0.05, seed = NULL){
  cl <- umap_clustering(cn, field = "copy", minPts = 5, umapmetric = "correlation",
                        seed = seed)
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

#' Phase-free imbalance per cell per unit
#'
#' Fits a symmetric two-component mixture to each cell's block counts: reads are
#' Binomial(n, q) or Binomial(n, 1-q) with equal prior, so the unknown block
#' phase is marginalised out. `q` is the allelic ratio and `lr` the evidence
#' against balance (q = 0.5).
#'
#' Use `lr`, not `q`, to rank cells: a thin cell can reach q = 0.03 on noise.
#' `lr` accounts for depth.
#'
#' Unlike the depth-stratified Tarone score this is ABSOLUTE, not relative to
#' other cells - which matters on samples where imbalance is common, because
#' there a relative score compares each imbalanced cell against other imbalanced
#' cells and reports that nothing stands out.
phase_free_imbalance <- function(blocks, mincounts = 2, minblocks = 30){
  mix <- function(a, n, q) sum(log(0.5 * dbinom(a, n, q) + 0.5 * dbinom(a, n, 1 - q)))
  blocks[totalcounts >= mincounts, {
    o <- optimize(function(q) -mix(allele0, totalcounts, q), c(0.02, 0.5))
    .(q = o$minimum, lr = -o$objective - mix(allele0, totalcounts, 0.5), blocks = .N)
  }, by = .(cell_id, unit)][blocks >= minblocks]
}

#' Seed the block phasing for the joint optimisation
#'
#' Per unit, take the block-level majority of the cell with the most imbalance
#' evidence; blocks it does not cover fall back to the all-cell majority. The
#' seed only has to give SOME cell a non-zero weight so the vote can bootstrap -
#' seeded from mirrored phasing every cell looks balanced after the HMM
#' discretises, every weight is zero, and the update has nothing to work with.
seed_phasing <- function(blocks, imb, key = c("chr","start","end","hap_label")){
  glob <- blocks[, .(s0 = sum(allele0), s1 = sum(allele1)), by = key]
  glob[, phi := as.integer(s0 >= s1)]
  glob <- glob[, c(key, "phi"), with = FALSE]
  best <- imb[order(-lr), .SD[1], by = unit]
  seeded <- rbindlist(lapply(seq_len(nrow(best)), function(i){
    x <- blocks[unit == best$unit[i] & cell_id == best$cell_id[i] & totalcounts > 0]
    x[, phi := as.integer(allele0 > allele1)]
    x[, c(key, "phi"), with = FALSE] }))
  rbind(seeded, glob[!seeded, on = key])
}

#' One joint-optimisation update of the block phases
#'
#' With phi_h = 1 meaning "raw allele0 is the A side", the log-likelihood ratio
#' for phi_h is
#'
#'   logLR_h = sum_c (allele0_ch - allele1_ch) * w_cb
#'
#' a vote over EVERY cell with no threshold and no cell selection. A cell that
#' is balanced in bin b contributes exactly zero; a cell holding the opposite
#' parent has w < 0, so its counts flip sign and it reinforces the vote instead
#' of cancelling - which is what the previous pooled-majority approach could not
#' do without explicitly orienting a hand-picked subset.
#'
#' The weight takes its SIGN from the per-bin HMM state and its MAGNITUDE from
#' the cell's phase-free ratio. The HMM's hard Viterbi call over-states
#' imbalance when overdispersion is high - on the noisier of two test samples it
#' called 326 of 369 cells imbalanced - and those spurious weights swamp the few
#' real ones. The phase-free ratio stays near 0.5 for a genuinely balanced cell
#' whatever the HMM said.
#'
#' @return `data.table` of the block key plus `phi` and `llr`; `abs(llr)` is the
#'   evidence for that block's phase, in log units.
joint_phase_update <- function(blocks, states, imb, loherror = 0.035,
                               key = c("chr","start","end","hap_label")){
  cl <- function(x) pmin(pmax(x, loherror), 1 - loherror)
  st <- as.data.table(states)[!is.na(A) & !is.na(B) & (A + B) > 0,
                              .(cell_id, chr, start, sgn = sign(A - B))]
  im <- imb[, .(cell_id, unit, mag = abs(log(cl(q) / (1 - cl(q)))))]
  b <- blocks[st, on = .(cell_id, chr, start), nomatch = 0]
  b <- b[im, on = .(cell_id, unit), nomatch = 0]
  b[, wt := sgn * mag]
  b[, .(llr = sum((allele0 - allele1) * wt)), by = key][
    , phi := as.integer(llr > 0)][]
}

#' Re-key a phasing list onto the detection units
#'
#' The first pass may phase by chromosome while detection runs per arm. In that
#' case every arm inherits its chromosome's cell choice, which is exactly the
#' situation being tested: "given the cells signals actually used here, is there
#' a better set for this arm?"
#'
#' @param phasing `hscn$phasing`, keyed by chromosome or by unit.
#' @param units The unit names detection will use.
expand_phasing_to_units <- function(phasing, units){
  stats::setNames(lapply(units, function(u){
    if (!is.null(phasing[[u]])) phasing[[u]] else phasing[[sub("[pq]$", "", u)]]
  }), units)
}

#' Is the installed signals new enough to honour phasebyarm?
#'
#' `phasebyarm` is a script-level flag but arm-level phasing is a package
#' change, and the script is bind-mounted while the package is baked into the
#' image - so the flag can take effect against an installed package that still
#' ignores it. In that combination `get_cells_per_chr_local()` keys the cell
#' list by chromosome while `phase_haplotypes_bychr()` filters on chromosome
#' arm, nothing matches, and the run dies much later with an opaque
#' "'no' is of type NULL" from `fifelse()` after the block join returns 0 rows.
#'
#' Probed behaviourally rather than by version, so a dev install is judged on
#' what it does. `ncells_for_clustering = 1` takes the branch that gives every
#' cell its own clone, so this needs neither uwot nor real data.
#'
#' `get_cells_per_chr_local` is internal to signals, so it has to be fetched
#' from the namespace - calling it unqualified fails with "could not find
#' function" and would make this report every build as unsupported.
#'
#' @param fn The function to probe. Defaults to signals' internal one; injected
#'   in tests so both a fixed and an old implementation can be exercised.
phasebyarm_supported <- function(fn = NULL){
  fail <- function(reason){ out <- FALSE; attr(out, "reason") <- reason; out }
  if (is.null(fn)){
    fn <- tryCatch(utils::getFromNamespace("get_cells_per_chr_local", "signals"),
                   error = function(e) NULL)
    if (is.null(fn))
      return(fail("get_cells_per_chr_local not found in the signals namespace"))
  }
  ascn <- data.table::data.table(
    cell_id = rep(c("a", "b"), each = 2), chr = "1",
    start = c(1, 150000001, 1, 150000001),          # chr1 p and q (hg19 cen ~121.5Mb)
    end = c(5e5, 150500000, 5e5, 150500000),
    state = 2L, LOH = "NO", state_BAF = 0.5, balance = 1L)
  haps <- data.table::copy(ascn)[, `:=`(hap_label = seq_len(.N), allele1 = 7L, allele0 = 2L,
                                        state = NULL, LOH = NULL, state_BAF = NULL,
                                        balance = NULL)]
  cl <- try(suppressMessages(
    fn(ascn, haps, ncells_for_clustering = 1, phasebyarm = TRUE)), silent = TRUE)
  if (inherits(cl, "try-error"))
    return(fail(paste0("probe call errored: ", conditionMessage(attr(cl, "condition")))))
  if (length(cl) == 0) return(fail("probe returned an empty cell list"))
  if (!all(grepl("[pq]$", names(cl))))
    return(fail(paste0("cell list is keyed by [", paste(names(cl), collapse = ", "),
                       "] rather than by chromosome arm")))
  TRUE
}

#' Label rows with the phasing unit
#'
#' The unit must match what `phase_haplotypes_bychr()` keys on: the chromosome,
#' or the chromosome arm when `phasebyarm = TRUE`. Everything downstream - the
#' imbalance score, the flagging, the re-phasing - works on `unit`, so switching
#' granularity is one argument rather than a parallel code path.
#'
#' @param d A data.frame with `chr` and `start`.
#' @param by `"chr"` (default) or `"arm"`.
#' @return `d` as a data.table with a `unit` column.
add_unit <- function(d, by = c("chr", "arm")){
  by <- match.arg(by)
  d <- as.data.table(d)
  if (by == "arm") d[, unit := paste0(chr, coord_to_arm(chr, start))] else d[, unit := chr]
  d[]
}

#' Per-cell, per-unit allelic imbalance from UNPHASED haplotype counts
#'
#' Tarone's Z tests whether a cell's per-block allele counts are overdispersed
#' relative to Binomial(n, 0.5). An allelically imbalanced cell is overdispersed
#' at the block level regardless of how the blocks are phased, so this detects
#' imbalance without any phasing, clustering or HMM. Only the ranking within a
#' chromosome is meaningful - Z is not calibrated as a p-value, because
#' genome-wide overdispersion inflates the baseline for every cell.
#'
#' @param haplotypes Wide-format haplotypes: `cell_id`, `chr`, `hap_label`,
#'   `allele0`, `allele1`, `totalcounts`.
#' Z grows with sequencing depth, so raw Z is NOT comparable between cells: a
#' shallow cell and a deep cell with identical allelic ratios get very different
#' Z. Each cell is therefore scored against cells of similar depth on the same
#' chromosome (robust z-score within a depth stratum), and it is `score`, not
#' `Z`, that should be compared or thresholded.
#'
#' @param mincounts Minimum counts for a block to contribute. Blocks with a
#'   single read carry no dispersion information. Default 2.
#' @param nstrata Number of depth strata per unit. Default 10.
#' @param by Phasing unit: `"chr"` (default) or `"arm"`. Scoring per chromosome
#'   dilutes a cell that is imbalanced over only one arm, which is how a
#'   chromosome like 4 - p-arm LOH in a cluster of cells, q-arm imbalance in one
#'   other cell - passes the chromosome-level check while the q arm is phased
#'   off noise.
#'
#' @return A `data.table` with `cell_id`, `unit`, `chr`, `Z`, `counts`,
#'   `nblocks`, `stratum` and `score`.
cell_imbalance_Z <- function(haplotypes, mincounts = 2, nstrata = 10, by = c("chr", "arm")){
  h <- add_unit(haplotypes, by)
  z <- h[totalcounts >= mincounts,
    {
      p <- sum(allele0) / sum(totalcounts)
      S <- sum((allele0 - totalcounts * p)^2 / (p * (1 - p)))
      .(Z = (S - sum(totalcounts)) / sqrt(2 * sum(totalcounts * (totalcounts - 1))),
        counts = sum(totalcounts), nblocks = .N, chr = chr[1])
    }, by = .(cell_id, unit)][is.finite(Z)]
  z[, stratum := cut(counts,
                     breaks = unique(quantile(counts, seq(0, 1, length.out = nstrata + 1))),
                     include.lowest = TRUE, labels = FALSE), by = unit]
  # Floor the per-stratum spread at a fraction of the unit-wide spread. A
  # stratum whose cells happen to have near-identical Z otherwise collapses the
  # MAD to ~0 and the ratio explodes, letting a cell with no dispersion at all
  # clear the threshold.
  z[, mad_unit := mad(Z), by = unit]
  z[, score := (Z - median(Z)) / pmax(mad(Z), 0.25 * mad_unit, 1e-6), by = .(unit, stratum)]
  z[, mad_unit := NULL]
  z[]
}

#' Pick the cells to phase one unit with
#'
#' Only cells whose depth-stratified imbalance score clears `min_score`, capped
#' at `kmax`. There is deliberately no padding to a minimum count: pooling
#' near-balanced cells into the vote actively destroys the phasing, so a unit
#' with no qualifying cell returns none and is left alone rather than re-phased
#' off noise.
select_phasing_cells <- function(Zcell, myunit, min_score = 10, kmax = 10){
  z <- Zcell[unit == myunit][order(-score)]
  if (nrow(z) == 0) return(character(0))
  head(z[score >= min_score]$cell_id, kmax)
}

#' Flag phasing units where the clustering missed the imbalanced cells
#'
#' `callHaplotypeSpecificCN` phases each chromosome from an hdbscan cluster, and
#' `minPts` is floored at `mincells`, so an event carried by only a handful of
#' cells can never be isolated into its own cluster. When that happens the
#' chromosome gets phased off read noise in near-balanced cells and every cell's
#' A/B assignment flips from bin to bin. The signature is that the cells signals
#' chose rank far below the most imbalanced cells available.
#'
#' @param hscn Result of `callHaplotypeSpecificCN` (needs `$phasing`).
#' @param Zcell Output of `cell_imbalance_Z`.
#' @param min_score Only flag chromosomes that have a convincingly imbalanced
#'   cell to begin with. Where no cell stands out from others of similar depth,
#'   no phasing is recoverable from the data and re-phasing would just fit noise.
#'   Default 10.
#' @param min_gap Flag when the best available score exceeds the best chosen
#'   cell's score by this margin. Default 10.
#'
#' `names(hscn$phasing)` are already the units signals phased on - chromosomes,
#' or arms when `phasebyarm = TRUE` - so this works on either granularity with
#' no branching, provided `Zcell` was built with the matching `by`.
#'
#' @return A `data.table` with one row per unit and a `flagged` column.
flag_failed_phasing <- function(hscn, Zcell, min_score = 10, min_gap = 10){
  rbindlist(lapply(names(hscn$phasing), function(cc){
    z <- Zcell[unit == cc][order(-score)]
    used <- unique(hscn$phasing[[cc]])
    zu <- z[cell_id %in% used]
    if (nrow(z) == 0) return(NULL)
    data.table(unit = cc,
               ncells_used = length(used),
               score_best = z$score[1],
               score_used_max = if (nrow(zu)) max(zu$score) else NA_real_,
               best_used_rank = if (nrow(zu)) min(match(zu$cell_id, z$cell_id)) else NA_integer_,
               n_candidates = z[score >= min_score, .N],
               ncells = nrow(z))
  }))[, gap := score_best - pmax(score_used_max, 0)][
      , flagged := is.finite(gap) & score_best >= min_score & gap >= min_gap][]
}

#' Put candidate phasing cells into a common orientation
#'
#' Which parent a cell over-represents is arbitrary and differs between cells:
#' two cells can both be strongly imbalanced over the same region while
#' favouring opposite haplotypes. Pooling their raw counts then makes them
#' cancel rather than reinforce, and the pooled majority becomes a patchwork of
#' whichever cell happens to dominate locally.
#'
#' So before pooling, orient every candidate against the highest-scoring one by
#' pairwise agreement on shared blocks, and flip those that disagree. Cells that
#' genuinely carry the opposite parent then contribute constructively instead of
#' eroding the signal.
#'
#' @param hc Haplotype rows for one unit.
#' @param cells Candidate cell ids, highest score first.
#' @param key Block key columns.
#' @return Named integer vector of +1 / -1, one per cell.
orient_cells <- function(hc, cells, key){
  sgn <- stats::setNames(rep(1L, length(cells)), cells)
  if (length(cells) < 2) return(sgn)
  maj <- function(cc){
    x <- hc[cell_id == cc & totalcounts >= 2, c(key, "allele0", "allele1"), with = FALSE]
    x[, m := as.integer(allele0 > allele1)]
    x[, c(key, "m"), with = FALSE]
  }
  ref <- maj(cells[1])
  for (cc in cells[-1]){
    j <- merge(ref, maj(cc), by = key, suffixes = c(".r", ".o"))
    if (nrow(j) < 10) next          # too little overlap to judge; leave as is
    if (mean(j$m.r == j$m.o) < 0.5) sgn[cc] <- -1L
  }
  sgn
}

#' Rebuild haplotype block phasing for the flagged units
#'
#' Unflagged units keep signals' own phase labels verbatim, so they
#' reproduce the original run exactly. Flagged units are re-phased from the
#' most imbalanced cells, with a fallback to all-cell pooling for blocks those
#' cells happen not to cover - `phase_haplotypes_bychr` omits uncovered blocks
#' entirely, and `format_haplotypes` then drops them from *every* cell, so
#' without the fallback a small phasing set silently discards a large share of
#' the chromosome's counts.
#'
#' @param by Phasing unit, matching the `by` used for `Zcell`.
#' @return List with `phased_haplotypes` (pass to `callHaplotypeSpecificCN`),
#'   `cells` (named list of the cells used per re-phased unit) and `rephased`.
rebuild_phasing <- function(hscn, haplotypes, Zcell, flagged,
                            min_score = 10, kmax = 10, by = c("chr", "arm")){
  by <- match.arg(by)
  key <- c("chr", "start", "end", "hap_label")
  haps <- add_unit(haplotypes, by)
  newrows <- list(); used <- list()
  for (cc in flagged){
    cells <- select_phasing_cells(Zcell, cc, min_score, kmax)
    if (length(cells) == 0){
      message(paste0("  ", cc, ": no cell clears the imbalance threshold, ",
                     "keeping the original phasing"))
      next
    }
    hc <- haps[unit == cc]
    sgn <- orient_cells(hc, cells, key)
    sub <- hc[cell_id %in% cells][, s := sgn[cell_id]][
      , .(a0 = sum(data.table::fifelse(s > 0L, allele0, allele1)),
          a1 = sum(data.table::fifelse(s > 0L, allele1, allele0))), by = key]
    glob <- hc[, .(a0 = sum(allele0), a1 = sum(allele1)), by = key]
    if (any(sgn < 0L))
      message(paste0("  ", cc, ": flipped ", sum(sgn < 0L), " of ", length(cells),
                     " phasing cells into a common orientation"))
    ph <- rbind(sub, glob[!sub, on = key])
    newrows[[cc]] <- ph[, .(chr, start, end, hap_label,
                            phase = fifelse(a0 < a1, "allele0", "allele1"))]
    used[[cc]] <- cells
    message(paste0("  ", cc, ": re-phased from ", length(cells), " cells, ",
                   nrow(ph), " blocks (", nrow(sub), " covered by those cells)"))
  }
  # Units we did not actually re-phase keep signals' own labels verbatim, so
  # they reproduce the original run exactly. Exclude by UNIT, not by chr: with
  # arms, `rephased` holds "4q" while haplotype_phasing only has chr, so
  # excluding on chr would keep the 4q rows as well and duplicate them.
  rephased <- names(newrows)
  base <- add_unit(hscn$haplotype_phasing, by)[!unit %in% rephased,
                                               c(key, "phase"), with = FALSE]
  list(phased_haplotypes = as.data.frame(rbind(base, rbindlist(newrows))),
       cells = used, rephased = rephased)
}

#' Minimum-breakpoint orientation for one cell on one chromosome
#'
#' Which allele is called A and which is B is an arbitrary label on two otherwise
#' unidentifiable sets of haplotype blocks. Swapping the label within a set of
#' bins leaves the emission likelihood *exactly* unchanged, because signals'
#' emission probabilities are symmetric: p(B) + p(state - B) = 1 for every state
#' and every B, so dbetabinom(b, n, p) == dbetabinom(n - b, n, 1 - p). The label
#' is therefore free to choose, and the only thing that distinguishes one choice
#' from another is how many transitions it implies.
#'
#' Two rules, applied in order:
#'  1. Within a run of constant total copy number the orientation is held fixed
#'     (`pmax`/`pmin`). A phase switch inside a segment - `2|1` to `1|2` at
#'     constant total 3 - is not identifiable from a single cell: it is the same
#'     model with the labels of the second part swapped. So it is never emitted.
#'  2. Across segment boundaries the relative orientation is chosen to minimise
#'     the number of *allele-specific* breakpoints, counting which alleles change
#'     rather than by how much. `2|1` -> `0|1` changes only A and costs 1;
#'     `2|1` -> `1|0` changes both and costs 2, so the former is preferred.
#'     Exact, by dynamic programming over two orientations per segment.
#'
#' @param state Integer vector of total copy number, in genomic order.
#' @param A,B Integer vectors of allele-specific copy number, in genomic order.
#' @return List with relabelled `A` and `B`.
relabel_min_breakpoints <- function(state, A, B, free = TRUE){
  n <- length(state)
  free <- rep_len(free, n)
  seg <- data.table::rleid(state)
  ns <- max(seg)
  # Decide free/pinned per COPY-NUMBER SEGMENT, by majority of its bins, and
  # mirror whole segments. Mirroring per bin against a scattered free set
  # forces some bins to A >= B while their pinned neighbours keep B > A, so
  # every free/pinned boundary becomes an A|B flip: on chr1q that produced 93
  # flips in a cell, all 93 of them sitting exactly on such a boundary.
  segfree <- vapply(seq_len(ns), function(j) mean(free[seg == j]) > 0.5, logical(1))
  mirrorable <- segfree[seg]
  Am <- ifelse(mirrorable, pmax(A, B), A)
  Bm <- ifelse(mirrorable, pmin(A, B), B)
  if (ns > 1 && any(segfree)){
    first_idx <- match(seq_len(ns), seg)
    last_idx  <- n - match(seq_len(ns), rev(seg)) + 1L
    av <- function(i, s) if (s == 1L) Am[i] else Bm[i]
    bv <- function(i, s) if (s == 1L) Bm[i] else Am[i]
    allowed <- function(j) if (segfree[j]) 1:2 else 1L
    dp <- matrix(Inf, ns, 2); bk <- matrix(1L, ns, 2)
    dp[1, allowed(1)] <- 0
    for (j in 2:ns){
      for (s in allowed(j)){
        opts <- vapply(allowed(j - 1L), function(sp)
          dp[j - 1L, sp] + (av(last_idx[j-1L], sp) != av(first_idx[j], s)) +
                            (bv(last_idx[j-1L], sp) != bv(first_idx[j], s)),
          numeric(1))
        dp[j, s] <- min(opts); bk[j, s] <- allowed(j - 1L)[which.min(opts)]
      }
    }
    orient <- integer(ns); orient[ns] <- which.min(dp[ns, ])
    for (j in ns:2) orient[j - 1L] <- bk[j, orient[j]]
    flip <- (orient == 2L)[seg]
    out <- list(A = ifelse(flip, Bm, Am), B = ifelse(flip, Am, Bm))
    # A global flip has identical cost and no meaning, so fix it by convention -
    # but only when nothing is pinned, since a pinned segment defines the frame.
    if (all(segfree)){
      big <- which.max(tabulate(seg, ns))
      if (sum(out$B[seg == big]) > sum(out$A[seg == big]))
        out <- list(A = out$B, B = out$A)
    }
    return(out)
  }
  list(A = Am, B = Bm)
}

#' Enforce the relabelling invariant across an hscn object
#'
#' `state_AS_phased` should differ from `state_AS` only where cross-cell evidence
#' actually pins the phase. Where it does not, the difference between the two
#' columns is accumulated noise that reads as real allelic switching. This walks
#' every (cell, chromosome) whose phase has no *external* support - no other cell
#' on that chromosome carries detectable allelic imbalance - and replaces the
#' labelling with the minimum-breakpoint one.
#'
#' A cell that is itself the only support for a chromosome counts as unsupported:
#' its own consistency there is circular, since the phasing was derived from it.
#'
#' Allele counts are swapped alongside the calls so that the returned table stays
#' internally consistent (BAF agrees with the called states). Note this makes the
#' per-cell BAF-vs-expected comparison uninformative for relabelled rows - but
#' that comparison is a function of the labelling anyway, so it was never
#' evidence about the copy number call. Use the labelling-invariant imbalance
#' score from [cell_imbalance_Z()] for that.
#'
#' @param hscn Result of `callHaplotypeSpecificCN`.
#' @param Zcell Output of [cell_imbalance_Z()].
#' @param min_score Imbalance score above which a cell counts as supporting the
#'   phase of its chromosome. Default 10.
#' @return The `hscn` object with `data` relabelled and a logical
#'   `phase_relabelled` column marking the affected rows.
relabel_unsupported_phase <- function(hscn, Zcell = NULL, min_score = 10,
                                     by = c("chr", "arm"), free_bins = NULL){
  by <- match.arg(by)
  d <- data.table::as.data.table(hscn$data)
  stopifnot(all(c("cell_id","chr","start","state","A","B") %in% names(d)))
  d <- d[order(cell_id, chr, start)]

  # `free_bins` (chr, start) lets the caller say directly which bins have an
  # undetermined phase. The joint optimisation supplies it from per-block |llr|,
  # which is the actual evidence; the older path derives it from a per-cell
  # score threshold, which is a proxy.
  if (!is.null(free_bins)){
    d <- add_unit(d, by)
    fb <- unique(as.data.table(free_bins)[, .(chr, start)])[, free := TRUE]
    d[, free := FALSE][fb, on = .(chr, start), free := TRUE]
    return(relabel_core(d, hscn, by))
  }

  # Support is evaluated per UNIT and passed to the DP as a per-bin flag, so a
  # chromosome with one supported arm still gets its unsupported arm relabelled.
  # Being conservative here (whole chromosome free or not) meant a single
  # supported arm blocked relabelling of the other one entirely.
  sup <- unique(Zcell[score >= min_score, .(cell_id, unit)])
  # deliberately NOT named `nsup`: data.table resolves a bare symbol in `i`
  # against the target's own columns first, so a join table sharing a column
  # name with `d` can silently be read as that column instead
  nsup_tbl <- sup[, .(nsup = .N), by = unit]
  d <- add_unit(d, by)
  d[, nsup := 0L][nsup_tbl, on = "unit", nsup := i.nsup]
  d[, is_sup := FALSE][sup, on = .(cell_id, unit), is_sup := TRUE]
  # A unit absent from Zcell has no haplotype data at all, which is not the same
  # as having data that nothing supports - there is nothing to phase and nothing
  # to relabel, so pin it. Acrocentric p-arms are the case in point: 14p has
  # bins in the CN table and zero haplotype counts, and relabelling them just
  # introduced a label break at the 14p/14q junction in every cell.
  units_with_data <- unique(Zcell$unit)
  d[, has_data := unit %in% units_with_data]
  # free when the unit has data and no OTHER cell supports its phase
  d[, free := has_data & (nsup - as.integer(is_sup)) == 0L]
  return(relabel_core(d, hscn, by))
}

#' Apply the minimum-breakpoint relabelling given a per-bin `free` flag
relabel_core <- function(d, hscn, by){
  runs <- function(x) sum(x != data.table::shift(x, fill = ""), na.rm = TRUE)
  before <- d[free == TRUE, .(n = runs(state_AS_phased)), by = .(cell_id, chr)][, sum(n)]

  d[, B0 := B]
  d[!is.na(A) & !is.na(B) & !is.na(state),
    c("A", "B") := relabel_min_breakpoints(state, A, B, free), by = .(cell_id, chr)]

  # A run of bins with no haplotype counts carries no evidence, and its A/B
  # values were filled from neighbours by `fillmissing` anyway. Left alone they
  # keep whatever orientation the HMM emitted, which need not match the segment
  # they adjoin - 13p is 36 count-free bins whose labels disagreed with 13q
  # across the centromere, costing one segment in 72 cells.
  #
  # This MUST run after the DP: done before, the DP then re-orients the free
  # segments and strands the count-free bins it had just been matched to.
  # Invents nothing - there is no data in these bins either way.
  if ("totalcounts" %in% names(d)){
    d[, informative := !is.na(totalcounts) & totalcounts > 0]
    d[, c("A", "B") := {
      if (any(informative) && any(!informative)){
        idx <- which(informative)
        nearest <- idx[max.col(-abs(outer(seq_len(.N), idx, "-")), ties.method = "first")]
        # Align the ORIENTATION of a count-free bin to its nearest informative
        # neighbour by swapping its own A and B - do not copy the neighbour's
        # values. Copying them took A/B from a bin with a different total copy
        # number, which broke A + B == state (2484 bins on the 725-cell sample,
        # e.g. state 2 receiving A=4,B=1 from a state-5 neighbour). A swap
        # preserves the bin's own total by construction.
        maj_nb <- sign(A[nearest] - B[nearest])
        maj_me <- sign(A - B)
        needswap <- !informative & maj_nb != 0 & maj_me != 0 & maj_nb != maj_me
        needswap[is.na(needswap)] <- FALSE
        .(ifelse(needswap, B, A), ifelse(needswap, A, B))
      } else .(A, B)
    }, by = .(cell_id, chr)]
    d[, informative := NULL]
  }

  d[, swapped := free & !is.na(B) & !is.na(B0) & B != B0]

  if (all(c("alleleA", "alleleB") %in% names(d))){
    d[swapped == TRUE, c("alleleA", "alleleB") := .(alleleB, alleleA)]
    # keep NA as NA rather than turning an empty bin into NaN
    if ("totalcounts" %in% names(d))
      d[, BAF := data.table::fifelse(is.na(totalcounts) | totalcounts == 0,
                                     NA_real_, alleleB / totalcounts)]
  }

  d <- add_states(d)
  after <- d[free == TRUE, .(n = runs(state_AS_phased)), by = .(cell_id, chr)][, sum(n)]
  message(paste0("Relabelled ", d[free == TRUE, uniqueN(paste(cell_id, chr))],
                 " (cell, chromosome) pairs with no external phase support; ",
                 "A/B runs ", before, " -> ", after,
                 " (", d[, sum(swapped, na.rm = TRUE)], " bins relabelled)"))
  d[, phase_relabelled := free]
  drop <- intersect(c("unit","nsup","is_sup","has_data","free","B0","swapped"), names(d))
  if (length(drop)) d[, (drop) := NULL]
  hscn$data <- as.data.frame(d)
  hscn
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
                     phased_haplotypes = NULL,
                     ncores = 1,
                     female = TRUE,
                     seed = NULL,
                     viterbiver = "cpp",
                     phasebyarm = FALSE,
                     allele_aware = FALSE){

  cn <- as.data.table(cn)
  qc <- as.data.table(qc)
  haps <- as.data.table(haps)
  
  if (filternormalcells == FALSE){
    message(paste0("Number of cells: ", length(unique(cn$cell_id))))
  } else if (filternormalcells == TRUE){
    message(paste0("Number of cells (pre filtering for normal cells): ", length(unique(cn$cell_id))))
    x <- find_tumor_cells(cn, diploidcutoff, seed = seed)
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
                                  phased_haplotypes = phased_haplotypes,
                                  ncores = ncores,
                                  female = female,
                                  seed = seed,
                                  viterbiver = viterbiver,
                                  phasebyarm = phasebyarm,
                                  allele_aware = allele_aware)
  
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
#' - `quality_filter`: Optional numeric threshold for filtering bins by quality (only bins with quality > this value are kept). Default 0.75.
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
                                mappability = args$mappability,
                                quality_filter = args$quality_filter)
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
  
  if (isTRUE(args$phasebyarm)){
    supported <- phasebyarm_supported()
    if (!isTRUE(supported)){
      stop("--phasebyarm requires a signals build in which get_cells_per_chr_local() ",
           "keys the cell list by chromosome arm (>= 0.17.0), otherwise the cell list ",
           "and the phasing filter disagree and the phasing table comes back empty.\n",
           "  installed signals: ", as.character(utils::packageVersion("signals")), "\n",
           "  probe result:      ", attr(supported, "reason"), "\n",
           "Rebuild the image, or drop --phasebyarm.", call. = FALSE)
    }
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
                  female = female,
                  seed = args$seed,
                  viterbiver = args$viterbiver,
                  phasebyarm = args$phasebyarm,
                  allele_aware = args$allele_aware)

  # Phasing unit, needed by both phasing paths and by the relabelling that
  # follows them - defined here rather than inside a branch.
  unit_by <- args$rephase_unit

  # Joint cross-cell optimisation of the block phases (roadmap C2). Replaces
  # the detect-flag-rebuild path entirely: no score threshold, no cell
  # selection, no orientation forcing. Every cell votes, weighted by its own
  # local imbalance, so cells carrying opposite parents reinforce rather than
  # cancel, and |llr| per block is the confidence.
  if (isTRUE(args$joint_phasing) && is.null(chr_cell_list_processed)){
    key <- c("chr","start","end","hap_label")
    haps_used <- add_unit(as.data.table(haplotypes)[cell_id %in% unique(res$data$cell_id)],
                          args$rephase_unit)
    message("Fitting phase-free imbalance per cell per unit...")
    imb <- phase_free_imbalance(haps_used)
    message(paste0("  ", nrow(imb), " (cell, unit) pairs; ",
                   imb[lr > 100, .N], " with strong imbalance evidence"))
    phi <- seed_phasing(haps_used, imb, key)
    prev <- NULL
    for (it in seq_len(args$joint_iters)){
      message(paste0("Joint phasing iteration ", it, " of ", args$joint_iters, "..."))
      res_it <- callhscn(cndata$cn, cndata$metrics, haplotypes,
                         frachaps = 0.8, maskedbins = maskedbins,
                         mincells = args$mincells, diploidcutoff = args$diploidcutoff,
                         ncores = args$ncores, filternormalcells = args$filternormalcells,
                         chrs_for_global_phasing = chrs_for_global_phasing,
                         global_phasing_for_balanced = global_phasing_for_balanced,
                         selftransitionprob = args$selftransitionprob,
                         phased_haplotypes = as.data.frame(phi[, c(key, "phi"), with = FALSE][
                           , .(chr, start, end, hap_label,
                               phase = fifelse(phi == 1L, "allele1", "allele0"))]),
                         female = female, seed = args$seed,
                         viterbiver = args$viterbiver, phasebyarm = args$phasebyarm,
                         allele_aware = args$allele_aware)
      upd <- joint_phase_update(haps_used, res_it$data, imb, loherror = res_it$loherror, key)
      if (!is.null(prev)){
        chg <- merge(prev[, c(key, "phi"), with = FALSE],
                     upd[, c(key, "phi"), with = FALSE], by = key)
        message(paste0("  blocks flipped since last iteration: ",
                       chg[phi.x != phi.y, .N], " of ", nrow(chg)))
        if (chg[phi.x != phi.y, .N] == 0) { res <- res_it; prev <- upd; break }
      }
      prev <- upd
      phi <- rbind(upd[, c(key, "phi"), with = FALSE],
                   seed_phasing(haps_used, imb, key)[!upd, on = key][, c(key,"phi"), with=FALSE])
      res <- res_it
    }
    # callHaplotypeSpecificCN returns phasing = NULL whenever phased_haplotypes
    # is supplied, and plot_clusters_used_for_phasing() iterates names() of it -
    # an empty list reaches cowplot::plot_grid(plotlist = list()) and kills the
    # run at the QC-plot stage, AFTER the inference and BEFORE the csv is
    # written. C2 has no "cells used for phasing" (every cell votes), so report
    # the cells carrying the most imbalance evidence in each unit.
    keep <- unique(res$data$cell_id)
    res$phasing <- stats::setNames(
      lapply(sort(unique(imb$unit)), function(u)
        head(imb[unit == u & cell_id %in% keep][order(-lr)]$cell_id, 5)),
      sort(unique(imb$unit)))
    res$phasing <- res$phasing[lengths(res$phasing) > 0]

    phase_llr <- prev[, .(chr, start, end, hap_label, llr)]
    res$phase_confidence <- phase_llr[, .(chr, start, end, hap_label, llr,
                                          abs_llr = abs(llr))]
    message("Per-unit phase confidence (median |logLR| per block):")
    pc <- add_unit(prev, args$rephase_unit)[, .(blocks = .N,
            median_abs_llr = round(median(abs(llr)), 2),
            frac_determined = round(mean(abs(llr) >= 5), 3)), by = unit]
    print(pc[order(frac_determined)])
    if (!is.null(args$rephase_report)) fwrite(pc, args$rephase_report)
  } else if (isTRUE(args$rephase_by_imbalance) && is.null(chr_cell_list_processed)){
    message("Checking per-chromosome phasing against unphased allelic imbalance")
    haps_used <- as.data.table(haplotypes)[cell_id %in% unique(res$data$cell_id)]
    # Detection granularity is independent of how the first pass clustered.
    # (unit_by is set above, before either branch.)
    # Arm-level CLUSTERING is worse - an arm has far fewer bins for UMAP, so the
    # first-pass cell choice degrades badly - but arm-level DETECTION is what
    # catches a chromosome like 4, where one arm is phased correctly and the
    # other is not. So default to a chromosome-level first pass with per-arm
    # detection, and let the two be set separately.
    Zcell <- cell_imbalance_Z(haps_used, by = unit_by)
    res_units <- res
    res_units$phasing <- expand_phasing_to_units(res$phasing, sort(unique(Zcell$unit)))
    diag <- flag_failed_phasing(hscn = res_units, Zcell = Zcell,
                                min_score = args$rephase_min_score,
                                min_gap = args$rephase_min_gap)
    print(diag[order(-gap)])
    if (!is.null(args$rephase_report)){
      fwrite(diag, args$rephase_report)
    }

    flagged <- diag[flagged == TRUE]$unit
    # guard against the column name and the call site drifting apart: reading a
    # column that does not exist yields NULL, which silently looks like
    # "nothing to re-phase" and skips the entire second pass
    stopifnot(length(flagged) == diag[flagged == TRUE, .N])
    rp <- if (length(flagged) > 0){
      message(paste0("Re-phasing ", length(flagged), " ",
                     ifelse(unit_by == "arm", "arm(s): ", "chromosome(s): "),
                     paste(flagged, collapse = ", ")))
      rebuild_phasing(res, haps_used, Zcell, flagged,
                      min_score = args$rephase_min_score,
                      kmax = args$rephase_maxcells,
                      by = unit_by)
    } else NULL

    if (is.null(rp) || length(rp$rephased) == 0){
      message("No chromosome needs re-phasing, keeping the original phasing")
    } else {
      res_rp <- callhscn(cndata$cn,
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
                         phased_haplotypes = rp$phased_haplotypes,
                         female = female,
                         seed = args$seed,
                         viterbiver = args$viterbiver,
                         phasebyarm = args$phasebyarm,
                         allele_aware = args$allele_aware)
      # callHaplotypeSpecificCN returns phasing = NULL when phased_haplotypes is
      # supplied; restore it so plot_clusters_used_for_phasing() still works and
      # reports the cells actually used.
      res_rp$phasing <- res$phasing
      for (cc in names(rp$cells)) res_rp$phasing[[cc]] <- rp$cells[[cc]]
      res_rp$rephasing_diagnostics <- diag
      res <- res_rp
    }

  }

  # Where no cross-cell evidence pins the phase, the A/B label is arbitrary:
  # replace it with the labelling that implies the fewest allele-specific
  # breakpoints. Likelihood-neutral by construction.
  #
  # This sits OUTSIDE the two phasing branches: it previously lived inside the
  # detect-flag-rebuild branch, so --joint_phasing skipped it entirely and the
  # unphaseable regions kept their raw flipping.
  if (isTRUE(args$relabel_unsupported_phase)){
    if (exists("phase_llr", inherits = FALSE)){
      # joint path: gate on the actual per-block evidence rather than a proxy
      # Judge a bin by the FRACTION of its blocks that are individually
      # determined, not by the summed |logLR|: a sum grows with block count, so
      # a bin of 30 weak blocks clears any fixed total while none of its blocks
      # is actually informative.
      fb <- phase_llr[, .(frac = mean(abs(llr) >= args$phase_llr_min)),
                      by = .(chr, start)][frac < args$phase_min_frac, .(chr, start)]
      nb <- phase_llr[, uniqueN(paste(chr, start))]
      message(paste0("Relabelling bins whose phase is undetermined (<",
                     round(100 * args$phase_min_frac), "% of blocks reach |logLR| ",
                     args$phase_llr_min, "): ", nrow(fb), " of ", nb, " bins"))
      res <- relabel_unsupported_phase(res, by = unit_by, free_bins = fb)
    } else if (exists("Zcell", inherits = FALSE)){
      message("Applying the relabelling invariant to unsupported chromosomes")
      res <- relabel_unsupported_phase(res, Zcell, min_score = args$rephase_min_score,
                                       by = unit_by)
    }
  }

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
  
  g6 <- tryCatch(plot_clusters_used_for_phasing(res) +
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      legend.background = element_rect(fill = "white"), # get rid of legend bg
      legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
    ),
    error = function(e){
      warning("plot_clusters_used_for_phasing() failed: ", conditionMessage(e),
              " - continuing without that panel", call. = FALSE)
      ggplot() + theme_void() })

  # Size the device from the number of phasing panels instead of a fixed 20in.
  # Phasing by arm roughly doubles the panel count (23 chromosomes -> up to 46
  # units), and a fixed height simply halves the space each row gets: at 6
  # columns the phasing grid goes from 4 rows of 2.5in to 7 rows of 1.4in, which
  # is what squashes the CN/BAF panels.
  phase_rows <- attr(g6, "phasing_rows")
  if (is.null(phase_rows)) phase_rows <- 4      # blank fallback panel
  h_top <- 5                                    # BAF- and variance-per-state
  h_mid <- 5                                    # the three histograms
  h_phase <- max(2.5 * phase_rows, 6)
  gall <- cowplot::plot_grid(
    cowplot::plot_grid(g1, g2, align = "h", axis = "tb"),
    cowplot::plot_grid(g3, g4, g5, ncol = 3),
    g6, ncol = 1, rel_heights = c(h_top, h_mid, h_phase))

  cowplot::save_plot(args$qcplot,
                     gall,
                     base_width = 15, base_height = h_top + h_mid + h_phase)
  
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
                        minPts = round(0.03 * length(cell_ids)),
                        seed = args$seed)
  
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
  parser$add_argument("--rephase_unit", default="chr", type="character",
                      help="Granularity for detecting and correcting bad phasing: 'chr' (default) or 'arm'. Independent of --phasebyarm, which changes how the FIRST pass clusters. Arm-level clustering is worse (fewer bins per unit for UMAP), so the useful combination is a chromosome-level first pass with --rephase_unit arm, which catches cases like chromosome 4 where one arm phases correctly and the other does not.")
  parser$add_argument("--phasebyarm", action="store_true", default=FALSE,
                      help="Phase each chromosome arm separately rather than the whole chromosome. Use where imbalance is confined to one arm: on the test sample chromosome 4 has p-arm LOH in a cluster of cells and q-arm imbalance in a different single cell, so whole-chromosome phasing gets the q arm wrong. Also switches the imbalance detection to per-arm scoring. Requires signals >= 0.17.0.")
  parser$add_argument("--allele_aware", action="store_true", default=FALSE,
                      help=paste("Score HMM transitions by how many alleles change",
                                 "(1[dA!=0] + 1[dB!=0]) instead of a flat matrix over B.",
                                 "Identical wherever total CN is constant; at a copy number",
                                 "change a one-allele move is preferred, so 2|1 -> 0|1 beats",
                                 "2|1 -> 1|0, and 1|1 -> 1|2 is no longer penalised relative",
                                 "to 1|1 -> 2|1."))
  parser$add_argument("--viterbiver", default="cpp", type="character",
                      help="Viterbi implementation: 'cpp' (default) or 'R'. Append _legacy ('cpp_legacy', 'R_legacy') to reproduce the incorrect backtrace used up to signals 0.16.0, which did not return the maximum a posteriori path and emitted spurious single-bin state changes. Requires signals >= 0.17.0.")
  parser$add_argument("--seed", default=NULL, type="integer",
                      help="Random seed. Phasing picks which cells phase each chromosome via UMAP + hdbscan and random subsampling, so without a seed two runs on identical input can phase a chromosome differently and return different haplotype-specific states. Set this to make a run reproducible. Requires signals >= 0.16.0.")
  parser$add_argument("--joint_phasing", action="store_true", default=FALSE,
                      help="Jointly optimise the block phases with the per-cell states (roadmap C2). Replaces --rephase_by_imbalance: every cell votes weighted by its own local imbalance, with no score threshold and no cell selection, so cells carrying opposite parental haplotypes reinforce rather than cancel. Emits per-block |logLR| as phase confidence.")
  parser$add_argument("--phase_llr_min", default=5, type="double",
                      help="Per-BLOCK |logLR| at which a block counts as determined; see --phase_min_frac for how bins are judged. Only used with --joint_phasing. Default 5.")
  parser$add_argument("--phase_min_frac", default=0.5, type="double",
                      help="A bin counts as having a determined phase when at least this fraction of its blocks reach |logLR| >= --phase_llr_min. Below it the A/B label is arbitrary and the minimum-breakpoint relabelling applies. Only used with --joint_phasing. Default 0.5.")
  parser$add_argument("--joint_iters", default=4, type="integer",
                      help="Maximum joint-phasing iterations; stops early once no block flips. Each iteration is a full HSCN inference. Default 4 (converged in 3 on both test samples).")
  parser$add_argument("--rephase_by_imbalance", action="store_true", default=FALSE,
                      help="After the first pass, detect chromosomes where the per-chromosome clustering selected cells with little allelic imbalance, and re-phase those chromosomes from the most imbalanced cells instead. Runs the HSCN inference a second time.")
  parser$add_argument("--rephase_min_score", default=10, type="double",
                      help="Minimum depth-stratified imbalance score for a cell to be used for phasing, and for a chromosome to be re-phased at all. Below this the cell is not distinguishable from others of similar sequencing depth, so re-phasing would fit noise. Default 10.")
  parser$add_argument("--rephase_min_gap", default=10, type="double",
                      help="Re-phase when the best available imbalance score exceeds the best score among the cells signals chose by this margin. Default 10.")
  parser$add_argument("--rephase_maxcells", default=10, type="integer",
                      help="Cap on cells used to re-phase a chromosome. Adding near-balanced cells to the vote degrades the phasing, so keep this small. Default 10.")
  parser$add_argument("--relabel_unsupported_phase", action="store_true", default=FALSE,
                      help="For (cell, chromosome) pairs where no other cell carries detectable allelic imbalance, replace the arbitrary A/B labelling with the one implying the fewest allele-specific breakpoints. Likelihood-neutral. Requires --rephase_by_imbalance (it reuses the imbalance scores). Adds a phase_relabelled column to the output.")
  parser$add_argument("--rephase_report", default=NULL, type="character",
                      help="Optional csv path for the per-chromosome phasing diagnostics table.")
  parser$add_argument("--mappability", default=NULL, type="double",
                      help="Mappability threshold; only bins with map > this value are kept. Default NULL (no filter).")
  parser$add_argument("--quality_filter", default=0.75, type="double",
                      help="Quality filter threshold; only bins with quality > this value are kept. Default 0.75.")
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
