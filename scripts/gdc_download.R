source(file.path(
  here::here(),
  'code','gdc.R'))

download_RNA <- FALSE
download_CNA <- FALSE
download_MAF <- FALSE
download_biospecimen <- FALSE
download_clinical_supp <- TRUE

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

gdc_supplement_cache_path <-
  file.path(
    data_raw_dir,
    "GDCdata",
    "clinical_supplement"
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

if(!dir.exists(gdc_supplement_cache_path)){
  dir.create(
    gdc_supplement_cache_path,
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
if(download_RNA){
  download_rna_tsv(
    gdc_rna_cache_path =
      file.path(
        data_raw_dir,
        "GDCdata",
        "rnaseq"
      ),
    gdc_projects = gdc_projects
  )
}

#####--- CNA -----######
if(download_CNA){
  download_cna_tsv(
     gdc_cna_cache_path =
       gdc_cna_cache_path,
     gdc_projects = gdc_projects
   )
}

#####--- MAFs -----#####
if(download_MAF){
  download_ssm_maf(
    gdc_maf_cache_path =
      gdc_maf_cache_path,
    gdc_projects = gdc_projects
  )
}

#####--- Biospecimen XML -----#####
if(download_biospecimen){
  download_biospecimen_xml(
    gdc_biospecimen_cache_path =
      gdc_biospecimen_cache_path
  )
}

#####--- Clinical Supplement XML -----#####
if(download_clinical_supp){
  download_clinical_supplement(
    gdc_supplement_cache_path =
      gdc_supplement_cache_path
  )
}
