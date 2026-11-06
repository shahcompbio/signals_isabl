#!/usr/bin/env Rscript
# Tests for the phasing helpers in run_signals.R. These live outside the package
# so they are not covered by tests/testthat. Run with:
#
#   Rscript test-run-signals-functions.R
#
# Needs data.table, dplyr and the signals package (for coord_to_arm), so in
# practice run it inside the signals image.

suppressMessages({library(data.table); library(dplyr); library(signals)})

want <- c("add_unit", "cell_imbalance_Z", "select_phasing_cells",
          "flag_failed_phasing", "rebuild_phasing", "relabel_min_breakpoints",
          "relabel_unsupported_phase", "expand_phasing_to_units", "orient_cells",
          "phase_free_imbalance", "seed_phasing", "joint_phase_update", "relabel_core")
# find run_signals.R beside this file, so the tests work from any directory
here <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script <- if (length(here)) file.path(dirname(here[1]), "run_signals.R") else "run_signals.R"
for (e in parse(script))
  if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
      as.character(e[[2]])[1] %in% want) eval(e)

fails <- 0L
ok <- function(cond, msg) {
  cat(if (isTRUE(cond)) "  PASS  " else "  FAIL  ", msg, "\n")
  if (!isTRUE(cond)) fails <<- fails + 1L
}

# --- synthetic chr4: p-arm imbalance in a cluster, q-arm imbalance in one cell -
set.seed(1)
ncells <- 40L
cells <- sprintf("s-l-R%02d-C%02d", seq_len(ncells), seq_len(ncells))
starts <- c(seq(1, by = 5e5, length.out = 40),           # p arm  (< 49.66Mb)
            seq(60000001, by = 5e5, length.out = 40))    # q arm
nblk <- 6L
haps <- CJ(cell_id = cells, start = starts, blk = seq_len(nblk))
haps[, `:=`(chr = "4", end = start + 5e5 - 1,
            hap_label = as.integer(factor(paste(start, blk))))]
parm <- haps$start < 49660117
loh_cluster <- cells[1:8]                       # imbalanced on p
q_cell      <- cells[20]                        # imbalanced on q only
# the true phase is a property of the BLOCK, shared by every cell
psi_tbl <- data.table(hap_label = sort(unique(haps$hap_label)))
psi_tbl[, psi := sample(0:1, .N, replace = TRUE)]
psi <- psi_tbl$psi[match(haps$hap_label, psi_tbl$hap_label)]
frac <- fifelse(haps$cell_id %in% loh_cluster & parm, 0.97,
        fifelse(haps$cell_id == q_cell & !parm, 0.97, 0.5))
p0 <- ifelse(psi == 1, frac, 1 - frac)
haps[, totalcounts := rpois(.N, 12)]
haps[, allele0 := rbinom(.N, totalcounts, p0)][, allele1 := totalcounts - allele0]
haps <- haps[totalcounts > 0][, blk := NULL]

# --- add_unit ----------------------------------------------------------------
ok(identical(sort(unique(add_unit(haps, "chr")$unit)), "4"),
   "add_unit(by='chr') gives one unit per chromosome")
ok(identical(sort(unique(add_unit(haps, "arm")$unit)), c("4p", "4q")),
   "add_unit(by='arm') splits at the centromere")

# --- scoring -----------------------------------------------------------------
Zc <- suppressMessages(cell_imbalance_Z(haps, by = "chr"))
Za <- suppressMessages(cell_imbalance_Z(haps, by = "arm"))
ok(all(Zc$unit == Zc$chr), "by='chr' leaves unit == chr")
ok(setequal(unique(Za$unit), c("4p", "4q")), "by='arm' scores each arm")

top_p <- Za[unit == "4p"][order(-score)][score >= 10]$cell_id
top_q <- Za[unit == "4q"][order(-score)][score >= 10]$cell_id
ok(length(top_p) > 0 && all(top_p %in% loh_cluster),
   "4p candidates are exactly the p-arm cluster")
ok(identical(top_q, q_cell), "4q candidate is the single q-arm cell")
ok(!q_cell %in% top_p, "the q-arm cell is not selected for 4p")

# --- selection ---------------------------------------------------------------
ok(identical(select_phasing_cells(Za, "4q"), q_cell),
   "select_phasing_cells picks the q-arm cell for 4q")
ok(length(select_phasing_cells(Za, "4q", min_score = 1e6)) == 0,
   "an unreachable threshold selects nobody rather than falling back")

# --- flagging ----------------------------------------------------------------
# mimic the failure: every arm inherits the chromosome's cell choice
hscn <- list(phasing = list("4p" = loh_cluster, "4q" = loh_cluster),
             haplotype_phasing = unique(haps[, .(chr, start, end, hap_label)])[
               , phase := "allele1"])
d <- flag_failed_phasing(hscn, Za)
ok(isTRUE(d[unit == "4q"]$flagged), "4q is flagged when phased from the p-arm cluster")
ok(isFALSE(d[unit == "4p"]$flagged), "4p is not flagged - its cells are already the best")

# --- rebuilding --------------------------------------------------------------
rp <- suppressMessages(rebuild_phasing(hscn, haps, Za, d[flagged == TRUE]$unit, by = "arm"))
ph <- as.data.table(rp$phased_haplotypes)
ok(identical(rp$rephased, "4q"), "only 4q is rebuilt")
ok(identical(rp$cells[["4q"]], q_cell), "4q is rebuilt from the q-arm cell")
ok(nrow(ph) == nrow(unique(ph, by = c("chr", "start", "end", "hap_label"))),
   "no duplicate block keys after rebuilding by arm")
ok(nrow(ph) == uniqueN(haps$hap_label), "every block still has exactly one phase")
ok(!anyNA(ph) && all(ph$phase %in% c("allele0", "allele1")), "phases are well formed")

# the rebuilt 4q phase should recover the true block orientation
qblocks <- unique(haps[start >= 49660117]$hap_label)
chk <- merge(ph[hap_label %in% qblocks, .(hap_label, phase)], psi_tbl, by = "hap_label")
acc <- max(mean(chk$phase == "allele0" & chk$psi == 1 | chk$phase == "allele1" & chk$psi == 0),
           mean(chk$phase == "allele0" & chk$psi == 0 | chk$phase == "allele1" & chk$psi == 1))
ok(acc > 0.9, sprintf("rebuilt 4q phase matches the true orientation (%.0f%%)", 100 * acc))

# --- joint cross-cell optimisation (C2) --------------------------------------
# The case the pooled-majority approach cannot handle: two cells imbalanced over
# the same region with OPPOSITE parents. Pooled raw they cancel; oriented by
# hand they can be forced together; here they must simply come out opposite.
local({
  set.seed(3)
  cl <- sprintf("s-l-R%02d-C%02d", 1:30, 1:30)
  st <- seq(1, by = 5e5, length.out = 40)
  bb <- data.table::CJ(cell_id = cl, start = st, k = 1:8)
  bb[, `:=`(chr = "4", end = start + 5e5 - 1,
            hap_label = as.integer(factor(paste(start, k))))]
  tru <- data.table(hap_label = sort(unique(bb$hap_label)))
  tru[, psi := sample(0:1, .N, replace = TRUE)]
  bb[, ps := tru$psi[match(hap_label, tru$hap_label)]]
  bb[, frac := fifelse(cell_id == cl[1], 0.8, fifelse(cell_id == cl[2], 0.2, 0.5))]
  bb[, totalcounts := rpois(.N, 14)]
  bb[, allele0 := rbinom(.N, totalcounts, fifelse(ps == 1, frac, 1 - frac))]
  bb[, allele1 := totalcounts - allele0]
  bb <- add_unit(bb[totalcounts > 0], "chr")
  kk <- c("chr", "start", "end", "hap_label")

  imb <- phase_free_imbalance(bb)
  ok(nrow(imb) == 30, "phase_free_imbalance: one row per cell")
  ok(all(imb[cell_id %in% cl[1:2]]$lr > imb[!cell_id %in% cl[1:2], max(lr)]),
     "the imbalanced cells rank top by EVIDENCE, not by point estimate")
  ok(max(abs(imb[cell_id %in% cl[1:2]]$q - 0.2)) < 0.06,
     "q recovers the ratio for both cells regardless of which parent")

  phi0 <- seed_phasing(bb, imb, kk)
  ok(nrow(phi0) == data.table::uniqueN(bb$hap_label), "the seed covers every block")

  mk <- function(phi){
    x <- bb[phi, on = kk, nomatch = 0]
    x[, `:=`(aA = fifelse(phi == 1L, allele0, allele1),
             aB = fifelse(phi == 1L, allele1, allele0))]
    y <- x[, .(A = sum(aA), B = sum(aB)), by = .(cell_id, chr, start)]
    y[, `:=`(A = round(2 * A / (A + B)))][, B := 2L - A][] }
  upd <- joint_phase_update(bb, mk(phi0), imb, 0.035, kk)
  ok(nrow(upd) == data.table::uniqueN(bb$hap_label), "the update scores every block")
  j <- merge(upd[, .(hap_label, phi)], tru, by = "hap_label")
  acc <- max(mean(j$phi == j$psi), mean(j$phi != j$psi))
  ok(acc > 0.9, sprintf("one update recovers the true block phase (%.0f%%)", 100 * acc))

  r <- bb[upd[, c(kk, "phi"), with = FALSE], on = kk, nomatch = 0]
  r[, `:=`(A = fifelse(phi == 1L, allele0, allele1),
           B = fifelse(phi == 1L, allele1, allele0))]
  r <- r[, .(baf = sum(B) / sum(A + B)), by = cell_id]
  s1 <- r[cell_id == cl[1]]$baf; s2 <- r[cell_id == cl[2]]$baf
  ok((s1 - 0.5) * (s2 - 0.5) < 0,
     sprintf("cells with opposite parents land on OPPOSITE sides (%.2f vs %.2f)", s1, s2))
  ok(abs(r[!cell_id %in% cl[1:2], mean(baf)] - 0.5) < 0.05,
     "balanced cells stay at BAF ~ 0.5 and contribute nothing")
  ok(all(c("llr", "phi") %in% names(upd)), "abs(llr) is emitted per block as confidence")
})

# --- pooling cells that carry opposite parents -------------------------------
# Two cells can both be strongly imbalanced over the same region while favouring
# opposite haplotypes. Pooled raw, they cancel and the majority becomes a
# patchwork; oriented first, they reinforce.
key <- c("chr", "start", "end", "hap_label")
qb <- unique(haps[start >= 49660117]$hap_label)
opp <- cells[21]                                  # second q-arm cell, opposite parent
hq <- copy(haps)
flip <- hq$cell_id == opp & hq$hap_label %in% qb
hq[flip, c("allele0", "allele1") := .(allele1, allele0)]
psi_q <- psi_tbl[hap_label %in% qb]
pflip <- psi_q$psi[match(hq$hap_label, psi_q$hap_label)]
sel <- hq$cell_id == opp & hq$hap_label %in% qb
hq[sel, totalcounts := rpois(.N, 12)]
hq[sel, allele0 := rbinom(.N, totalcounts, ifelse(pflip[sel] == 1, 0.03, 0.97))]
hq[sel, allele1 := totalcounts - allele0]
hq <- hq[totalcounts > 0]

two <- c(q_cell, opp)
sgn <- orient_cells(hq[hap_label %in% qb], two, key)
ok(length(sgn) == 2 && sgn[[1]] == 1L, "the reference cell keeps orientation +1")
ok(sgn[[2]] == -1L, "a cell carrying the opposite parent is flipped")

pool <- function(signs){
  x <- hq[hap_label %in% qb & cell_id %in% names(signs)][, s := signs[cell_id]]
  x[, .(a0 = sum(fifelse(s > 0L, allele0, allele1)),
        a1 = sum(fifelse(s > 0L, allele1, allele0))), by = key][
     , .(hap_label, phase = fifelse(a0 < a1, "allele0", "allele1"))]
}
acc <- function(ph){
  j <- merge(ph, psi_q, by = "hap_label")
  max(mean(j$phase == "allele0" & j$psi == 1 | j$phase == "allele1" & j$psi == 0),
      mean(j$phase == "allele0" & j$psi == 0 | j$phase == "allele1" & j$psi == 1))
}
a_raw <- acc(pool(setNames(rep(1L, 2), two)))
a_ori <- acc(pool(sgn))
ok(a_ori > a_raw + 0.15,
   sprintf("orienting recovers the true phase far better (%.0f%% vs %.0f%% unoriented)",
           100 * a_ori, 100 * a_raw))
ok(length(orient_cells(hq, q_cell, key)) == 1, "a single candidate needs no orientation")

# --- wiring between flagging and rebuilding ----------------------------------
# These exist because the two were tested in isolation while the line joining
# them read the wrong column name: flag_failed_phasing returns `unit`, the call
# site asked for `$chr`, that yielded NULL, and the whole second pass was
# silently skipped for an entire run.
ok("unit" %in% names(d), "flag_failed_phasing returns a `unit` column")
ok(!"chr" %in% names(d), "and not a `chr` column the caller might reach for")
flagged_units <- d[flagged == TRUE]$unit
ok(length(flagged_units) == d[flagged == TRUE, .N],
   "extracting flagged units yields one name per flagged row (NULL would be 0)")
rp2 <- suppressMessages(rebuild_phasing(hscn, haps, Za, flagged_units, by = "arm"))
ok(identical(rp2$rephased, "4q"), "the extracted names drive rebuild_phasing end to end")

# --- re-keying a chromosome-keyed phasing list onto arms ----------------------
by_chr <- list("4" = loh_cluster)
ex <- expand_phasing_to_units(by_chr, c("4p", "4q"))
ok(setequal(names(ex), c("4p", "4q")), "expand_phasing_to_units keys by unit")
ok(identical(ex[["4p"]], loh_cluster) && identical(ex[["4q"]], loh_cluster),
   "each arm inherits its chromosome's cells")
ok(identical(expand_phasing_to_units(list("4p" = "a", "4q" = "b"), c("4p","4q"))[["4q"]], "b"),
   "an already unit-keyed list is passed through unchanged")
# and detection on the re-keyed list still finds 4q
d2 <- flag_failed_phasing(list(phasing = ex), Za)
ok(isTRUE(d2[unit == "4q"]$flagged) && isFALSE(d2[unit == "4p"]$flagged),
   "chromosome-level first pass + arm-level detection flags 4q only")

# --- minimum-breakpoint relabelling ------------------------------------------
r <- relabel_min_breakpoints(state = c(3,3,1,1), A = c(2,1,1,0), B = c(1,2,0,1))
ok(identical(paste0(r$A, "|", r$B), c("2|1","2|1","0|1","0|1")),
   "relabel prefers 2|1 -> 0|1 (one allele changes) over 2|1 -> 1|0")
r2 <- relabel_min_breakpoints(state = rep(3, 4), A = c(2,1,2,1), B = c(1,2,1,2))
ok(length(unique(paste0(r2$A, "|", r2$B))) == 1,
   "relabel removes phase flips within a constant-CN segment")

# pinned segments: a supported arm must keep its labels while an unsupported one
# is still cleaned up. Being conservative here (all-or-nothing per chromosome)
# meant one supported arm blocked relabelling of the other entirely.
st <- c(2,2,3,3,3,3)                      # segment 1 = pinned, segment 2 = free
Ai <- c(2,0,1,2,1,2); Bi <- c(0,2,2,1,2,1)
rp3 <- relabel_min_breakpoints(st, Ai, Bi, free = c(FALSE,FALSE,TRUE,TRUE,TRUE,TRUE))
ok(identical(rp3$A[1:2], Ai[1:2]) && identical(rp3$B[1:2], Bi[1:2]),
   "pinned bins keep their labels exactly")
ok(length(unique(paste0(rp3$A[3:6], "|", rp3$B[3:6]))) == 1,
   "the free segment is still made internally consistent")
rp4 <- relabel_min_breakpoints(st, Ai, Bi, free = FALSE)
ok(identical(rp4$A, Ai) && identical(rp4$B, Bi),
   "with nothing free, relabelling is a no-op")
rp5 <- relabel_min_breakpoints(c(3,3,1,1), c(2,1,1,0), c(1,2,0,1), free = TRUE)
ok(identical(paste0(rp5$A, "|", rp5$B), c("2|1","2|1","0|1","0|1")),
   "free = TRUE still reproduces the unpinned result")

# --- gating relabelling on evidence instead of a score proxy ------------------
# The joint optimisation supplies free bins from per-block |logLR|. This block
# previously lived inside the detect-flag-rebuild branch, so --joint_phasing
# skipped relabelling entirely and unphaseable regions kept their raw flipping.
local({
  dd <- data.table(cell_id = rep(c("a","b"), each = 6), chr = "4",
                   start = rep(seq(1, by = 5e5, length.out = 6), 2))
  dd[, `:=`(end = start + 5e5 - 1, state = 3L,
            A = c(2L,1L,2L,1L,2L,2L, 2L,2L,1L,2L,1L,2L),
            B = c(1L,2L,1L,2L,1L,1L, 1L,1L,2L,1L,2L,1L),
            totalcounts = 30L, alleleA = 20L, alleleB = 10L)]
  dd[, state_AS_phased := paste0(A, "|", B)]
  nruns <- function(x) sum(x != data.table::shift(x, fill = x[1]))

  all_free <- unique(dd[, .(chr, start)])
  o1 <- as.data.table(suppressMessages(relabel_unsupported_phase(
    list(data = as.data.frame(dd)), by = "chr", free_bins = all_free)$data))
  ok(o1[cell_id == "a", nruns(state_AS_phased)] == 0 &&
     o1[cell_id == "b", nruns(state_AS_phased)] == 0,
     "free_bins covering everything removes every arbitrary flip")

  o2 <- as.data.table(suppressMessages(relabel_unsupported_phase(
    list(data = as.data.frame(dd)), by = "chr",
    free_bins = data.table(chr = character(), start = numeric()))$data))
  ok(identical(o2$state_AS_phased, dd$state_AS_phased),
     "empty free_bins is an exact no-op")

  part <- unique(dd[start < 1.2e6, .(chr, start)])
  o3 <- as.data.table(suppressMessages(relabel_unsupported_phase(
    list(data = as.data.frame(dd)), by = "chr", free_bins = part)$data))
  ok(all(o3[start >= 1.2e6]$state_AS_phased == dd[start >= 1.2e6]$state_AS_phased),
     "bins with a determined phase keep their labels exactly")
  ok(sum(o3$phase_relabelled) == nrow(part) * 2,
     "only the undetermined bins are relabelled")
})

# --- a scattered free mask must not create flips ------------------------------
# Deciding free/pinned per BIN and then mirroring per bin forces free bins to
# A >= B while their pinned neighbours keep B > A, so every boundary becomes an
# A|B flip. With the per-bin gate the free set is scattered, and on chr1q this
# turned four cleanly-phased cells into 88-94 segments each - every flip sitting
# exactly on a free/pinned boundary.
local({
  n <- 60
  x <- data.table(cell_id = "a", chr = "1", start = seq(1, by = 5e5, length.out = n))
  x[, `:=`(end = start + 5e5 - 1, state = 3L, A = 2L, B = 1L)]   # one clean 2|1 segment
  set.seed(9)
  scattered <- runif(n) < 0.4                      # free set with many runs
  nruns <- function(v) sum(v != c(v[1], v[-length(v)]))
  r <- relabel_min_breakpoints(x$state, x$A, x$B, free = scattered)
  ok(nruns(paste0(r$A, "|", r$B)) == 0,
     "a scattered free mask over one clean segment produces no flips")
  ok(nruns(scattered) > 10, "the mask really is scattered (sanity check)")

  # and a segment that is mostly free is still made consistent
  y <- copy(x)[, `:=`(A = rep(c(2L, 1L), n / 2), B = rep(c(1L, 2L), n / 2))]
  r2 <- relabel_min_breakpoints(y$state, y$A, y$B, free = rep(TRUE, n))
  ok(nruns(paste0(r2$A, "|", r2$B)) == 0,
     "an all-free flipping segment is still collapsed to one orientation")

  # a mostly-pinned segment is left alone entirely
  r3 <- relabel_min_breakpoints(y$state, y$A, y$B, free = c(rep(TRUE, 5), rep(FALSE, n - 5)))
  ok(identical(r3$A, y$A) && identical(r3$B, y$B),
     "a segment that is mostly pinned is left untouched, not part-mirrored")
})

# --- units with no haplotype data at all --------------------------------------
# Acrocentric p-arms (13p, 14p, 15p, 22p) have CN bins and zero haplotype
# counts. "no cell supports this unit" and "this unit has no data" are different
# things: relabelling the second introduced a label break at the p/q junction in
# every cell.
cnbins <- rbind(
  data.table(cell_id = rep(cells[1:4], each = 4), chr = "14",
             start = rep(c(1, 500001, 20000001, 20500001), 4)),
  NULL)
cnbins[, `:=`(end = start + 5e5 - 1, state = 3L,
              A = rep(c(2L,1L,2L,1L), 4), B = rep(c(1L,2L,1L,2L), 4),
              totalcounts = rep(c(0L,0L,40L,40L), 4),
              alleleA = 0L, alleleB = 0L)]
cnbins[, state_AS_phased := paste0(A, "|", B)]
# haplotype data exists only on 14q (start >= 20Mb here); 14p has none
Zsyn <- data.table(cell_id = cells[1:4], unit = "14q", chr = "14", score = 0)
hs <- list(data = as.data.frame(cnbins))
out <- suppressMessages(relabel_unsupported_phase(hs, Zsyn, min_score = 10, by = "arm"))
od <- as.data.table(out$data)
ok(all(od[start < 20000000]$phase_relabelled == FALSE),
   "bins in a unit with no haplotype data are left alone")
ok(all(od[start >= 20000000]$phase_relabelled == TRUE),
   "bins in an unsupported unit that does have data are still relabelled")

# A pinned run of count-free bins must not keep an internal flip that the HMM
# happened to emit there (13p: 36 bins 0|1 then 3 bins 1|0). It should inherit
# from the nearest informative bin - no evidence either way, so no break.
nb <- data.table(cell_id = "a", chr = "14",
                 start = seq(1, by = 5e5, length.out = 6))
nb[, `:=`(end = start + 5e5 - 1, state = 3L,
          A = c(1L, 1L, 2L, 2L, 2L, 2L),      # internal flip in the count-free run
          B = c(2L, 2L, 1L, 1L, 1L, 1L),
          totalcounts = c(0L, 0L, 0L, 40L, 40L, 40L),
          alleleA = 0L, alleleB = 0L)]
nb[, state_AS_phased := paste0(A, "|", B)]
Znb <- data.table(cell_id = "a", unit = "14q", chr = "14", score = 0)
o2 <- as.data.table(suppressMessages(
  relabel_unsupported_phase(list(data = as.data.frame(nb)), Znb, by = "arm")$data))[order(start)]
ok(length(unique(o2$state_AS_phased)) == 1,
   "a count-free run inherits one orientation instead of keeping an internal flip")
ok(all(o2[totalcounts > 0]$state_AS_phased == nb[totalcounts > 0]$state_AS_phased),
   "and the informative bins are untouched")


# --- relabel_core must not break A + B == state -------------------------------
# The count-free bin alignment used to copy A/B from the nearest informative
# bin, which imported a different total copy number.
cat("\n-- relabel_core preserves total copy number --\n")
{
  d <- data.table::data.table(
    cell_id = "c1", chr = "1", start = seq(1, by = 5e5, length.out = 6),
    state   = c(5L, 5L, 2L, 2L, 3L, 3L),
    A       = c(4L, 4L, 1L, 1L, 2L, 2L),
    B       = c(1L, 1L, 1L, 1L, 1L, 1L),
    totalcounts = c(30L, 25L, 0L, 0L, 20L, 18L),
    alleleA = c(20L, 17L, 0L, 0L, 13L, 12L),
    alleleB = c(10L,  8L, 0L, 0L,  7L,  6L),
    free = TRUE)
  d[, end := start + 5e5 - 1]
  d[, BAF := ifelse(totalcounts > 0, alleleB / totalcounts, NA_real_)]
  d[, state_AS_phased := paste0(A, "|", B)]
  d[, state_min := pmin(A, B)]
  out <- relabel_core(data.table::copy(d), hscn = NULL, by = "chr")
  ok(all(out$A + out$B == out$state),
        "A + B == state in every bin after relabel_core")
  ok(all(out$state == d$state), "state itself is unchanged")
  ok(all(out$A >= 0) && all(out$B >= 0), "no negative alleles")
  # the count-free bins (3,4) keep their own total of 2, not the neighbour's 5
  cf <- out$totalcounts == 0
  ok(all(out$A[cf] + out$B[cf] == 2L),
     "count-free bins keep their own total copy number, not the neighbour's")
}

cat("\nfailures:", fails, "\n")
quit(status = if (fails > 0L) 1L else 0L)
