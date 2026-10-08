# Phasing improvements in `run_signals.R`

Work done on poorly phased chromosomes — the case where only one or two cells in
a sample are non-diploid over a region, and those cells come out with noisy,
alternating `1|2`/`2|1` calls.

Validated on two DLP+ samples:

* **S1** `SHAH_H002216_T01_01_DLP01` — 725 cells, overdispersion ~0.028
* **S2** `SHAH_H003831_T01_01_DLP01` — 369 cells, overdispersion ~0.077

Everything here is **off by default**; the script's behaviour without new flags
is unchanged except for the two fixes noted under *Fixes*.

---

## Recommended invocation

```
Rscript run_signals.R \
  --seed 42 \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --rephase_report confidence.csv \
  ...
```

This is the configuration that scored best on both samples. Note what is **not**
in it: `--relabel_unsupported_phase`. See *Do not use*.

---

## How phasing quality was measured

Segment counts and mean-BAF-versus-expected are **not** used, because both are
functions of the A/B labelling and so reward aggressive relabelling regardless of
whether a call is right. On S1 the two metrics point in opposite directions: the
most correct configuration has the *most* segments.

The metric used instead is a **phase-free pairwise test**. For two cells that are
both imbalanced over a region, it asks whether their allele counts are better
explained by the two cells favouring the same parental haplotype or opposite
ones, marginalising the unknown block phase:

```
logLR = sum_h log[ 0.5 B(a_ih; n_ih, q_i) B(a_jh; n_jh, q_j)
                 + 0.5 B(a_ih; n_ih, 1-q_i) B(a_jh; n_jh, 1-q_j) ]
      - sum_h log[ 0.5 B(a_ih; n_ih, q_i) B(a_jh; n_jh, 1-q_j)
                 + 0.5 B(a_ih; n_ih, 1-q_i) B(a_jh; n_jh, q_j) ]
```

where `q` is each cell's phase-free allelic ratio. It uses no phasing
information, so its answer cannot move when the phasing changes. `logLR` scales
with block count and depth, so it is calibrated by a permutation that breaks the
block pairing while keeping both marginals; the reported statistic is `z`, and
`|z| >= 5` counts as determined.

On S1 this gives 1,751 well-determined cell pairs across 20 chromosome arms; on
S2, 1,631 across 13.

---

## Changes

### `--seed`
Phasing has three stochastic steps: the subsampling in `min_cells` that sets the
cluster size, the UMAP embedding that picks which cells phase each chromosome,
and the subsampling in the beta-binomial fit. Left unseeded, two runs on
identical input can select different phasing cells and return different
haplotype-specific states. On S1, chromosome 6 ranked 226th by imbalance in one
run and 2nd in the next on the same data.

Threaded through the script's own `find_tumor_cells` as well as into signals.

### `--rephase_unit {chr,arm}`
The unit over which phasing is assessed and re-derived. `arm` matters where
imbalance is confined to one arm: on S1, chromosome 4 has 17 cells with p-arm
LOH that phase the p arm correctly but are balanced on q, so q-arm blocks were
phased off noise.

Two analysis traps worth knowing, both of which produced wrong conclusions here:

* **Never summarise per chromosome when phasing per arm.** A genuine phase
  switch at the centromere averages to BAF ~0.5 across the whole chromosome. On
  S2 chromosome 19 this looked like failed phasing; per arm, both cells were
  cleanly separated and correct.
* **Check arm boundaries against the cytoband map.** A hand-picked centromere
  coordinate (19.5Mb for chr13, whose centromere is ~18Mb) misread a p/q
  junction as an internal flip.

### `--joint_phasing`, `--joint_iters`, `--phase_llr_min`, `--phase_min_frac`
The main addition. Rather than picking a subset of cells to phase each unit,
**every** cell votes on every haplotype block, weighted by its own imbalance:

```
LLR_h = sum_c (a0_ch - a1_ch) * sgn_c * |logit(q_c)|
```

`q_c` is cell `c`'s phase-free allelic ratio over the unit, so a balanced cell
has `sgn = 0` and contributes nothing, and a cell in LOH contributes heavily.
Iterated with the inference to a fixed point (`--joint_iters`, typically
converging in 2–4).

This is what fixed the original complaint: a region where only one or two cells
are imbalanced is now phased from those cells' counts rather than from a cluster
chosen by a feature that cannot see them.

Correctness, against the phase-free verdict:

| configuration | S1 | S2 |
|---|---|---|
| before | — | 1601/1631 (98%) |
| arm re-phasing + relabelling | 94% | 89% |
| joint phasing + relabelling | 99.5% | — |
| **joint phasing, no relabelling** | **1748/1748 (100%)** | **1630/1631 (100%)** |

`--rephase_report` writes a per-unit confidence table (`blocks`,
`median_abs_llr`, `frac_determined`), which is how to tell a well-phased unit
from one that merely looks phased.

### `--allele_aware`
Passes through to signals' allele-aware HMM transition cost. The HMM state is
the B allele copy number with `A = state - B`, so a flat matrix over B makes "B
unchanged" the cheap self-transition and gives "A unchanged" no standing:
`1|1 -> 2|1` was cheap while `1|1 -> 1|2` was penalised, purely because of how
the state space is parameterised. Scoring by how many alleles change is
symmetric in A and B.

Left correctness unchanged on both samples and reduced allele-specific segments
by 1.3–1.5%. Requires a signals build with the `allele_aware` argument.

### `--phasebyarm`
Passes through to signals' own `phasebyarm`. Note this was unusable in signals
before 0.18.0 — the per-unit cell list was keyed by chromosome while the
haplotype filter used chromosome arm, so `"4p" == "4"` never matched and zero
blocks were phased. The script guards against running it on a build that
silently ignores it.

---

## Fixes

Two bugs fixed in code this script owns. Both are on by default because the old
behaviour was simply wrong.

* **`relabel_core` broke total copy number.** The alignment of bins with no
  haplotype counts copied `A`/`B` from the nearest informative bin, importing
  that bin's total: a `state 2` bin could receive `A=4, B=1` from a state-5
  neighbour, affecting 2,484 bins on S1. The intent was only to align
  *orientation*, which is a swap, not a copy. Now swaps the bin's own `A`/`B`,
  preserving `A + B == state` by construction.
* **The QC plot could kill a completed run.** `plot_clusters_used_for_phasing`
  errored when `$phasing` was absent (which happens whenever
  `phased_haplotypes` is supplied), after the inference had finished but
  *before* the csv was written. Now wrapped in `tryCatch` and `$phasing` is
  repopulated.

---

## Do not use

### `--relabel_unsupported_phase` (and `--rephase_min_score`, `--rephase_min_gap`, `--rephase_maxcells`)

Where nothing pins the phase, the A/B label is arbitrary, so this picks the
labelling with the fewest allele-specific breakpoints — exactly minimising
`1[dA != 0] + 1[dB != 0]` by dynamic programming over the two orientations per
copy number segment, which is likelihood-neutral.

**The criterion is sound; the policy is not.** Mirroring undetermined bins to
`A >= B` makes two cells that gained *opposite* parents look identical. Measured
on S1, on opposite-parent pairs specifically:

| | opposite-parent pairs correct |
|---|---|
| with relabelling | 90% |
| without | **100%** |

All 8 disagreements in the relabelled run were opposite-parent pairs reported as
same-parent — on 1p, 2q, 5q, 9q, 18p, at `|z|` up to 22.6 over 1,345 blocks. The
otherwise-identical run without relabelling gets all 8 right. S2 agrees: 60% with,
100% without.

The flags are retained so the result is reproducible, not because they are
recommended. The principled replacement is to put the same preference inside the
likelihood where evidence can outweigh it — which is what `--allele_aware` does.

### `--rephase_by_imbalance` and the detect-flag-rebuild path
Superseded by `--joint_phasing`, which needs no score threshold, no cell
selection and no orientation forcing. Retained only for reproducing earlier runs.

---

## Known limitation: regions with too little evidence

Some regions cannot be phased by any method here, and it is worth being able to
tell. On S1, chromosome 8 looks wrong on the heatmap in every configuration
tried. The cause is not a phasing failure:

* only **2 of 726** cells are imbalanced on chr8
* both carry real signal (one fits `q = 0.050`, genuine LOH, `lr = 79.9`)
* but counts from voting cells on 8p total **41,470** against **3,608,059** on
  13q — roughly 1%
* worse, the cell with the strong imbalance contributes only 266 counts over 76
  blocks, a median of 3 per block

That is about 1 nat per block, roughly 2.7:1 odds — far below determination. No
transition cost or voting scheme manufactures information. The right handling is
to read `--rephase_report` and treat such a unit as undetermined rather than
trusting its labels.
