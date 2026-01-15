## GDC omics data retrieval scripts
source('code/utils.R')
source('code/gdc.R')
source('code/clinical.R')
source('code/mutation.R')
source('code/msi.R')
source('code/cna.R')
source('code/rnaseq.R')

gdc_projects <- paste0(
  "TCGA-", c(
  "ACC","BLCA","BRCA",
  "CESC","CHOL","COAD",
  "DLBC","ESCA","GBM",
  "HNSC","KICH","KIRC",
  "KIRP","LAML","LGG",
  "LIHC","LUAD","LUSC",
  "MESO","OV","PAAD",
  "PCPG","PRAD","READ",
  "SARC","SKCM","STAD",
  "TGCT","THCA","THYM",
  "UCEC","UCS","UVM"
))

msi_report_template_quarto <-
  file.path(here::here(),"code","msi_classifier.qmd")

gdc_release <- "release45_20251204"
output_dir <- file.path(
  here::here(), "output"
)

if(!dir.exists(
  file.path(output_dir, gdc_release))){
  dir.create(
    file.path(output_dir, gdc_release))
}

data_raw_dir <- file.path(
  here::here(), "data-raw"
)
gdc_biospecimen_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "biospecimen"
  )

gencode_xref <- get_gencode_xref(
  data_raw_dir = data_raw_dir
)


#####--- Clinical -----#####
tcga_clinical_info <- gdc_tcga_clinical(
  gdc_release = gdc_release,
  overwrite = F,
  tumor_samples_only = T,
  output_dir = output_dir)


####--- Mutation -----####
tcga_calls <- gdc_tcga_snv(
  gdc_projects = gdc_projects,
  output_dir = output_dir,
  gdc_release = gdc_release,
  data_raw_dir = data_raw_dir,
  gencode_xref = gencode_xref,
  overwrite = F
)

####--- Gene mutation rate ----####
tcga_gene_mutation_rate <- calculate_tcga_gene_mutation_rate2(
  gdc_projects = gdc_projects,
  gdc_release = gdc_release,
  output_dir = output_dir
)

#####--- VCF output -----#####
write_tcga_vcf2(
  tcga_calls = tcga_calls,
  output_dir = output_dir,
  gdc_release = gdc_release)


####--- TMB -----####
tmb_stats <- calculate_sample_tmb2(
  tcga_calls = tcga_calls,
  t_depth_min = 30,
  t_vaf_min = 0.05,
  gdc_release = gdc_release,
  overwrite = TRUE,
  output_dir = output_dir)

####--- MSI status ----####
msi_data <- gdc_tcga_msi(
  output_dir = output_dir,
  gdc_release = gdc_release,
  data_raw_dir = data_raw_dir,
  overwrite = FALSE)

####--- CNA -----#####
for(project in gdc_projects){
  tcga_cna <- gdc_tcga_cna(
    gdc_project = project,
    output_dir = output_dir,
    gdc_release = gdc_release,
    data_raw_dir = data_raw_dir,
    gencode_xref = gencode_xref,
    overwrite = TRUE
  )
}

####--- RNASeq -----#####
for(project in gdc_projects){
  calls <- gdc_tcga_rnaseq(
    gdc_project = project,
    output_dir = output_dir,
    gdc_release = gdc_release,
    data_raw_dir = data_raw_dir,
    gencode_xref = gencode_xref,
    overwrite = F
  )
}
