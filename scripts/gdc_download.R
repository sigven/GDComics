source(file.path(here::here(),'code','gdc.R'))

data_raw_dir <-
  file.path(
    here::here(),
    "data-raw"
  )

gdc_biospecimen_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "biospecimen"
  )

gdc_maf_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "maf"
  )

gdc_cna_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "cna"
  )

gdc_rna_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "rnaseq"
  )

if(!dir.exists(gdc_biospecimen_cache_path)){
  dir.create(
    gdc_biospecimen_cache_path,
    recursive = T
  )
}


if(!dir.exists(gdc_maf_cache_path)){
  dir.create(
    gdc_maf_cache_path,
    recursive = T
  )
}

if(!dir.exists(gdc_rna_cache_path)){
  dir.create(
    gdc_rna_cache_path,
    recursive = T
  )
}

if(!dir.exists(gdc_cna_cache_path)){
  dir.create(
    gdc_cna_cache_path,
    recursive = T
  )
}
#####--- RNA -----#####
# download_rna_tsv(
#   gdc_rna_cache_path =
#     file.path(
#       data_raw_dir,
#       "GDCdata",
#       "rnaseq"
#     ),
#   tcga_projects = tcga_projects
# )

#####--- CNA -----#####
#
download_cna_tsv(
   gdc_cna_cache_path =
     gdc_cna_cache_path,
   tcga_projects = tcga_projects
 )

#####--- MAFs -----#####
# download_ssm_maf(
#   gdc_maf_cache_path =
#     gdc_maf_cache_path,
#   tcga_projects = tcga_projects
# )

#####--- Biospecimen XML -----#####
# download_tcga_biospecimen_xml(
#   gdc_biospecimen_cache_path =
#     gdc_biospecimen_cache_path
# )
