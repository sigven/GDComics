# HRD (scarHRD-style) scores for TCGA-BRCA, OV, PRAD, PAAD

Exploratory run of `pcgr/hrd.py` (PCGR `dev`, v2.3.3) on GDC ASCAT3 allele-specific
segments (GRCh38), annotated with BRCA1/2 copy-number state and somatic LoF mutations.

Run from `output/release46_20260810/hrd_tcga/` (scripts read/write the working directory):

1. `Rscript ../../../code/hrd_tcga/export.R`  - segments_{project}.tsv (not kept, ~regenerable) + sample_meta.tsv
2. `python3 ../../../code/hrd_tcga/run_hrd.py` - hrd_scores_raw_<project>.tsv (~6 min for all four cohorts, 8 procs; per-project results are cached)
3. `python3 ../../../code/hrd_tcga/merge.py`   - tcga_hrd_scores.tsv (+ summary to stdout)
4. `Rscript ../../../code/hrd_tcga/plot.R`     - hrd_tcga_cohorts.png

`run_hrd.py` imports `pcgr.hrd` from a hard-coded PCGR checkout path and uses the 20261001
bundle's `chromsize.grch38.tsv`; adjust both at the top of the script if they move.
BRCA1/2 status: CN state of the overlapping segment (HOMDEL / HEMDEL / cnLOH / intact) and
somatic LoF mutations only - no germline variants or promoter methylation.
