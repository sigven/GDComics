## Select biomarker-rich TCGA sample sets (one per tumor type) for PCGR testing
## Output sample lists can be used directly with scripts/prepare_tcga_pcgr_batch.py
source('code/sample_selection.R')

gdc_release <- "release46_20260810"

sample_sets <- select_tcga_samples(
  tumors = c("STAD", "LUAD", "PRAD", "BRCA", "COAD",
             "PAAD", "SKCM", "BLCA", "UCEC", "HNSC"),
  biomarkers_tsv = file.path(
    here::here(), "data-raw", "sample_lists", "tcga_biomarkers_curated.tsv"),
  output_dir = file.path(here::here(), "output"),
  gdc_release = gdc_release,
  sample_list_dir = file.path(
    here::here(), "data-raw", "sample_lists", gdc_release),
  n_samples = 10,
  n_controls = 1,
  min_purity = 0.4,
  target_purity = 0.7,
  tmb_hypermutated = 20,
  max_hypermutated = 2
)
