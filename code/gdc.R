



# tcga_msi_raw <- GenomicDataCommons::files() |>
#   GenomicDataCommons::filter(
#      ~ data_type == "Aligned Reads"
#   ) |>
#   GenomicDataCommons::select(
#     c("file_id",
#       "submitter_id",
#       "data_type",
#       "data_category",
#       "data_format",
#       "cases.case_id",
#       "cases.submitter_id",
#       "cases.samples.submitter_id",
#       "cases.samples.sample_id",
#       "cases.samples.sample_type",
#       "analysis.workflow_type",
#       "analysis.workflow_version",
#       "analysis.analysis_type",
#       "analysis.created_datetime",
#       "msi_score",
#       "msi_status")) |>
#   GenomicDataCommons::response_all()
#
# msi_df <- tcga_msi_raw$results |>
#   tibble::as_tibble() |>
#   dplyr::mutate(
#     case_df = purrr::map(cases, get_first_df),
#     sample_df = purrr::map(case_df, get_first_sample)) |>
#   dplyr::transmute(
#     file_id,
#     file_submitter_id = submitter_id,
#     msi_score,
#     msi_status,
#     analysis_workflow_type = analysis$workflow_type,
#     analysis_type = analysis$analysis_type,
#     analysis_created_datetime = analysis$created_datetime,
#     analysis_workflow_version = analysis$workflow_version,
#     case_id      = purrr::map_chr(
#       case_df, ~ .x$case_id %||% NA_character_),
#     bcr_patient_barcode = purrr::map_chr(
#       case_df, ~ .x$submitter_id %||% NA_character_),
#     sample_id      = purrr::map_chr(
#       sample_df, ~ get_scalar(.x, "sample_id")),
#     tumor_sample_barcode = purrr::map_chr(
#       sample_df, ~ get_scalar(.x, "submitter_id")),
#     sample_type    = purrr::map_chr(
#       sample_df, ~ get_scalar(.x, "sample_type"))
#   ) |>
#   as.data.frame() |>
#   dplyr::filter(
#     !is.na(msi_score) &
#       stringr::str_detect(bcr_patient_barcode,"TCGA-"))
# get_first_df <- function(x) {
#   if (is.data.frame(x)) return(x)
#   if (is.list(x) && length(x) > 0 && is.data.frame(x[[1]])) return(x[[1]])
#   NULL
# }
# get_first_sample <- function(case_df) {
#   if (!is.data.frame(case_df)) return(NULL)
#   if (!"samples" %in% names(case_df)) return(NULL)
#
#   s <- case_df$samples[[1]]
#   if (is.data.frame(s)) return(s)
#   if (is.list(s) && length(s) > 0 && is.data.frame(s[[1]])) return(s[[1]])
#   NULL
# }
# get_scalar <- function(x, field) {
#   if (is.null(x)) return(NA_character_)
#   v <- x[[field]]
#   if (length(v) == 0) return(NA_character_)
#   as.character(v[[1]])
# }


#' Function that downloads GDC Gene Level Copy Number Estimates
#' using ASCAT3 workflow (TCGA projects only)
#'
#' @param gdc_cna_cache_path character, path to GDC MAF download cache
#' @param tcga logical, whether to download TCGA MAF files only
#' @param tcga_projects character vector, TCGA project IDs
#'
#'
download_cna_tsv <- function(
    tcga = TRUE,
    gdc_cna_cache_path = NULL,
    tcga_projects = NULL){

  assertthat::assert_that(
    !is.null(gdc_cna_cache_path)
  )
  assertthat::assert_that(
    !is.null(tcga_projects)
  )
  assertthat::assert_that(
    dir.exists(gdc_cna_cache_path))

  cna_GDCfilesResults <-
    GenomicDataCommons::files() |>

    ## filter on Masked Somatic Mutation / MAF files
    GenomicDataCommons::filter(
      ~ data_format == "TSV" &
        data_type == "Gene Level Copy Number" &
        analysis.workflow_type ==
        'ASCAT3') |>
    GenomicDataCommons::select(
      c("file_id",
        "data_type",
        "data_format",
        "data_category",
        "submitter_id",
        "analysis.workflow_type",
        "cases.project.project_id")) |>
    GenomicDataCommons::results_all()

  ## make data frame from complex nested list
  gdc_cna_df <-
    as.data.frame(
      dplyr::bind_rows(
        cna_GDCfilesResults$cases, .id = "id")) |>
    dplyr::mutate(project_id = project$project_id) |>
    dplyr::select(-c("project")) |>
    ## only keep TCGA projects
    dplyr::filter(
      project_id %in% tcga_projects
    ) |>
    dplyr::inner_join(
      data.frame(
        workflow_type = cna_GDCfilesResults$analysis$workflow_type,
        data_format = cna_GDCfilesResults$data_format,
        id = cna_GDCfilesResults$id,
        file_id = cna_GDCfilesResults$file_id,
        submitter_id = cna_GDCfilesResults$submitter_id,
        data_type = cna_GDCfilesResults$data_type
      ),
      by = "id"
    ) |>
    #TCGA-5P-A9KF-01A
    dplyr::mutate(
      tumor_sample_barcode = stringr::str_match(
        submitter_id,
        "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z0-9]{1}")[,1]
      )

  ## Download MAF files for each TCGA project
  for(project in tcga_projects){

    ## Create project directory in GDC cache path
    if(!(dir.exists(
      file.path(
        gdc_cna_cache_path, project)))){
      dir.create(
        file.path(
          gdc_cna_cache_path, project),
        recursive = T)
    }

    ## filter biospecimen XML data frame for current project
    project_cna_df <-
      gdc_cna_df |>
      dplyr::filter(
        project_id == project
      )


    ## set path to project GDC download cache
    GenomicDataCommons::gdc_set_cache(
      directory = file.path(
        gdc_cna_cache_path, project))

    ## download biospecimen XML files
    purrr::walk(
      project_cna_df$file_id,
      ~ GenomicDataCommons::gdcdata(
        .x,
        use_cached = T,
        progress = F
      )
    )

    system(
      paste0(
        "gzip ",
        file.path(
          gdc_cna_cache_path,
          project,
          "*",
          "*.tsv"
        )
      )
    )
  }

}

#' Function that downloads GDC Gene Level Copy Number Estimates
#' using ASCAT3 workflow (TCGA projects only)
#'
#' @param gdc_rna_cache_path character, path to GDC MAF download cache
#' @param tcga logical, whether to download TCGA MAF files only
#' @param tcga_projects character vector, TCGA project IDs
#'
#'
download_rna_tsv <- function(
    tcga = TRUE,
    gdc_rna_cache_path = NULL,
    tcga_projects = NULL){

  assertthat::assert_that(
    !is.null(gdc_rna_cache_path)
  )
  assertthat::assert_that(
    !is.null(tcga_projects)
  )
  assertthat::assert_that(
    dir.exists(gdc_cna_cache_path))

  rna_GDCfilesResults <-
    GenomicDataCommons::files() |>

    ## filter on Masked Somatic Mutation / MAF files
    GenomicDataCommons::filter(
      ~ data_format == "TSV" &
        data_type == "Gene Expression Quantification" &
        analysis.workflow_type ==
        'STAR - Counts') |>
    ## EXPAND to reach samples / aliquots
    GenomicDataCommons::expand(
      c("cases.samples")
    ) |>
    GenomicDataCommons::select(
      c(
        # file submitter ID
        "file_id",
        "submitter_id",
        "data_type",
        "data_format",
        "data_category",
        "analysis.workflow_type",

        ## case-level
        "cases.case_id",
        "cases.project.project_id",
        "cases.submitter_id"
      )) |>
    GenomicDataCommons::results_all()

  gdc_rna_sample_info <-
    rna_GDCfilesResults |>
    tibble::as_tibble() |>
    dplyr::mutate(
      tumor_sample_barcode =
        purrr::map(
          cases, ~ .x$samples |>
            dplyr::bind_rows() |>
            dplyr::pull(submitter_id))
    ) |>
    dplyr::select(-c("cases")) |>
    tidyr::unnest(tumor_sample_barcode) |>
    dplyr::mutate(
      workflow_type = analysis$workflow_type)

  gdc_rna_sample_info$analysis <- NULL

  gdc_rnaseq_metadata <-
    as.data.frame(
      dplyr::bind_rows(
        rna_GDCfilesResults$cases, .id = "id")) |>
    dplyr::mutate(project_id = project$project_id) |>
    dplyr::select(-c("project","samples")) |>
    dplyr::rename("bcr_patient_barcode" = "submitter_id") |>
    ## only keep TCGA projects
    dplyr::filter(
      project_id %in% tcga_projects
    ) |>
    dplyr::inner_join(
      gdc_rna_sample_info,
      by = "id"
    ) |>
    dplyr::select(-c("case_id","submitter_id")) |>
    dplyr::select(
      c("id",
        "bcr_patient_barcode",
        "tumor_sample_barcode",
        "project_id",
        "file_id",
        "data_type",
        "data_format",
        "data_category",
        "workflow_type")
    )

  readr::write_tsv(
    gdc_rnaseq_metadata,
    file = file.path(
      gdc_rna_cache_path,
      "rnaseq_gdc.metadata.tsv.gz"
    )
  )

  ## Download MAF files for each TCGA project
  for(project in tcga_projects){

    ## Create project directory in GDC cache path
    if(!(dir.exists(
      file.path(
        gdc_rna_cache_path, project)))){
      dir.create(
        file.path(
          gdc_rna_cache_path, project),
        recursive = T)
    }

    ## filter biospecimen XML data frame for current project
    project_rna_df <-
      gdc_rna_df |>
      dplyr::filter(
        project_id == project
      )


    ## set path to project GDC download cache
    GenomicDataCommons::gdc_set_cache(
      directory = file.path(
        gdc_rna_cache_path, project))

    ## download biospecimen XML files
    purrr::walk(
      project_rna_df$file_id,
      ~ GenomicDataCommons::gdcdata(
        .x,
        use_cached = T,
        progress = F
      )
    )

    system(
      paste0(
        "gzip ",
        file.path(
          gdc_rna_cache_path,
          project,
          "*",
          "*.tsv"
        )
      )
    )
  }

}


#' Function that downloads GDC Simple Somatic Mutation MAF files (TCGA projects only)
#'
#' @param gdc_maf_cache_path character, path to GDC MAF download cache
#' @param tcga logical, whether to download TCGA MAF files only
#' @param tcga_projects character vector, TCGA project IDs
#'
#'
download_ssm_maf <- function(
    tcga = TRUE,
    gdc_maf_cache_path = NULL,
    tcga_projects = c("TCGA-READ","TCGA-COAD",
                      "TCGA-STAD","TCGA-UCEC")){

  maf_GDCfilesResults <-
    GenomicDataCommons::files() |>

    ## filter on Masked Somatic Mutation / MAF files
    GenomicDataCommons::filter(
      ~ data_format == "maf" &
        data_type == "Masked Somatic Mutation" &
        analysis.workflow_type ==
        'Aliquot Ensemble Somatic Variant Merging and Masking') |>
    GenomicDataCommons::select(
      c("file_id",
        "data_type",
        "data_format",
        "data_category",
        "submitter_id",
        "analysis.workflow_type",
        "cases.project.project_id")) |>
    GenomicDataCommons::results_all()

  ## make data frame from complex nested list
  gdc_maf_df <-
    as.data.frame(
      dplyr::bind_rows(
        maf_GDCfilesResults$cases, .id = "id")) |>
    dplyr::mutate(project_id = project$project_id) |>
    dplyr::select(-c("project")) |>
    ## only keep TCGA projects
    dplyr::filter(
      project_id %in% tcga_projects
    ) |>
    dplyr::inner_join(
      data.frame(
        workflow_type = maf_GDCfilesResults$analysis$workflow_type,
        data_format = maf_GDCfilesResults$data_format,
        id = maf_GDCfilesResults$id,
        file_id = maf_GDCfilesResults$file_id,
        submitter_id = maf_GDCfilesResults$submitter_id,
        data_type = maf_GDCfilesResults$data_type
      ),
      by = "id"
    )

  ## Download MAF files for each TCGA project
  for(project in tcga_projects){

    ## Create project directory in GDC cache path
    if(!(dir.exists(
      file.path(
        gdc_maf_cache_path, project)))){
      dir.create(
        file.path(
          gdc_maf_cache_path, project),
        recursive = T)
    }

    ## filter biospecimen XML data frame for current project
    project_maf_df <-
      gdc_maf_df |>
      dplyr::filter(
        project_id == project
      )


    ## set path to project GDC download cache
    GenomicDataCommons::gdc_set_cache(
      directory = file.path(
        gdc_maf_cache_path, project))

    ## download biospecimen XML files
    purrr::walk(
      project_maf_df$file_id,
      ~ GenomicDataCommons::gdcdata(
        .x,
        use_cached = T,
        progress = F
      )
    )
  }

}

#' Function that loads GDC Simple Somatic Mutation MAF files
#' (TCGA projects only)
#'
#' @param gdc_projects character vector, GDC/TCGA project IDs to load
#' @param gdc_maf_cache_path character, path to GDC MAF download cache
#'
#' @export
load_ssm_maf <- function(
    gdc_maf_cache_path = NULL,
    gdc_projects = NULL){

  assertthat::assert_that(
    !is.null(gdc_maf_cache_path)
  )
  assertthat::assert_that(
    !is.null(gdc_projects)
  )
  assertthat::assert_that(
    dir.exists(gdc_maf_cache_path))
  maf_df_list <- list()

  for(project in gdc_projects){

    maf_fnames <-
      list.files(
        path = file.path(
          gdc_maf_cache_path, project),
        pattern = "\\.maf.gz$",
        full.names = T,
        recursive = T
      )

    project_maf_df <-
      purrr::map_dfr(
        maf_fnames,
        ~ readr::read_tsv(
          .x,
          comment = "#",
          col_types = readr::cols(.default = "c")) |>
          dplyr::mutate(
            gdc_maf_fname = .x,
            project_id = project
          )
      )

    maf_df_list[[project]] <- project_maf_df

  }

  maf_df <- dplyr::bind_rows(maf_df_list)

  return(maf_df)

}



#' Function that downloads GDC Biospecimen supplement XML files (TCGA projects only)
#'
#' @param gdc_biospecimen_cache_path character, path to GDC download cache
#' @param tcga logical, whether to download TCGA biospecimen XML files
#' @param tcga_projects character vector, TCGA project IDs to download biospecimen
#'
#'
download_biospecimen_xml <- function(
    tcga = TRUE,
    gdc_biospecimen_cache_path = NULL,
    tcga_projects = c("TCGA-READ","TCGA-COAD",
                      "TCGA-STAD","TCGA-UCEC")){

    biospecimen_GDCfilesResults <-
      GenomicDataCommons::files() |>

      ## filter on biospecimen supplement files
      GenomicDataCommons::filter(
        ~ data_category == "biospecimen" &
          data_type == "Biospecimen Supplement") |>
      GenomicDataCommons::select(
        c("file_id","data_type",
          "submitter_id","cases.project.project_id")) |>
      GenomicDataCommons::results_all()

    ## make data frame from complex nested list
    biospecimen_xml_df <-
      as.data.frame(
        dplyr::bind_rows(
          biospecimen_xml_GDCfilesResults$cases, .id = "id")) |>
      dplyr::mutate(project_id = project$project_id) |>
      dplyr::select(-c("project")) |>
      ## only keep TCGA projects
      dplyr::filter(
        project_id %in% tcga_projects
      ) |>
      dplyr::inner_join(
        data.frame(
          id =
            biospecimen_xml_GDCfilesResults$id,
          file_id =
            biospecimen_xml_GDCfilesResults$file_id,
          submitter_id =
            biospecimen_xml_GDCfilesResults$submitter_id,
          data_type =
            biospecimen_xml_GDCfilesResults$data_type
        ),
        by = "id"
      ) |>

      ## only keep biospecimen files
      dplyr::filter(
        stringr::str_detect(
          .data$submitter_id,"org_biospecimen.TCGA"))

    ## Download biospecimen XML files for each TCGA project
    for(project in tcga_projects){

      ## Create project directory in GDC cache path
      if(!(dir.exists(
        file.path(
          gdc_biospecimen_cache_path, project)))){
        dir.create(
          file.path(
            gdc_biospecimen_cache_path, project),
          recursive = T)
      }

      ## filter biospecimen XML data frame for current project
      project_xml_df <-
        biospecimen_xml_df |>
        dplyr::filter(
          project_id == project
        )


      ## set path to project GDC download cache
      GenomicDataCommons::gdc_set_cache(
        directory = file.path(
          gdc_biospecimen_cache_path, project))

      ## download biospecimen XML files
      purrr::walk(
        project_xml_df$file_id,
        ~ GenomicDataCommons::gdcdata(
          .x,
          use_cached = T,
          progress = F
        )
      )
    }

  }


#' Function that parses GDC Biospecimen supplement XML file
#'
#' Returns data frame with columns:
#' - bcr_patient_barcode
#' - tumor_sample_barcode
#' - project
#' - msi_status
#' - sample_type
#'
#' @param xml_fname character, path to GDC biospecimen XML file
#'
#' @export
#'
parse_biospecimen_xml <- function(xml_fname = NULL){

  df <- data.frame()

  if(file.exists(xml_fname) &
     endsWith(xml_fname,".xml")){

    #cat(xml_fname, "\n")
    bcr_patient_barcode_fname <-
      stringr::str_replace(
        basename(xml_fname),
        "nationwidechildrens.org_ssf\\.|\\.xml","")

    msi_status <- NA
    bcr_patient_barcode <- NA
    tumor_sample_barcode <- NA
    bio_sample_type_id <- NA
    bio_vial_number <- NA
    sample_type <- NA
    x <- xml2::read_xml(xml_fname)

    msi_status <- xml2::xml_find_first(
      x, ".//bio:msi_mono_di_nucleotide_assay_status") |>
      xml2::xml_text(trim = TRUE)

    bcr_patient_barcode <- xml2::xml_find_first(
      x, ".//shared:bcr_patient_barcode") |>
      xml2::xml_text(trim = TRUE)

    bio_sample_type_id <- xml2::xml_find_first(
      x, ".//bio:sample_type_id") |>
      xml2::xml_text(trim = TRUE)

    bio_vial_number <- xml2::xml_find_first(
      x, ".//bio:vial_number") |>
      xml2::xml_text(trim = TRUE)

    sample_type <- xml2::xml_find_first(
      x, ".//bio:sample_type") |>
      xml2::xml_text(trim = TRUE)

    if(!is.na(bcr_patient_barcode) &
       !is.na(bio_sample_type_id) &
       !is.na(bio_vial_number)){
      tumor_sample_barcode <- paste0(
        bcr_patient_barcode,"-",
          bio_sample_type_id,
          bio_vial_number)
    }

    if(is.na(msi_status)){
      cat("No MSI status found in ", xml_fname, "\n")
    }

    df <- data.frame(
      'bcr_patient_barcode' = bcr_patient_barcode,
      'tumor_sample_barcode' = tumor_sample_barcode,
      'msi_status' = msi_status,
      'sample_type' = sample_type,
      stringsAsFactors = F
    )
  }


  return(df)

}

#' Function that retrieves MSI status from GDC Biospecimen supplement XML files (TCGA projects only)
#'
#' @param gdc_biospecimen_cache_path character, path to GDC download cache
#' @param tcga_projects character vector, TCGA project IDs to download biospecimen
#'
#' @export
#'
#'
get_msi_status <- function(
    gdc_biospecimen_cache_path = NULL,
    tcga_projects = c("TCGA-READ","TCGA-COAD",
                      "TCGA-STAD","TCGA-UCEC")){

  assertthat::assert_that(
    !is.null(gdc_biospecimen_cache_path),
    !is.null(tcga_projects)
  )
  assertthat::assert_that(
    dir.exists(gdc_biospecimen_cache_path))

  msi_status <- list()
  for(project in tcga_projects){

    project_xml_fnames <-
      list.files(
        path = file.path(
          gdc_biospecimen_cache_path, project),
        pattern = "\\.xml$",
        full.names = T,
        recursive = T
      )

    project_msi_status_df <-
      purrr::map_df(
        project_xml_fnames,
        parse_biospecimen_xml
      )

    project_msi_status_df$project <- project

    msi_status[[project]] <- project_msi_status_df
  }

  return(dplyr::bind_rows(msi_status))

}

#' Function that classifies GDC Gene Level Copy Number Estimates
#' from ASCAT3 into categorical CNA calls
#'
#' @param cna_tsv_fname character, path to GDC CNA TSV file
#' @param gencode_xref data frame, gene cross-reference data frame
#' @param ignore_neutral logical, whether to ignore neutral calls
#'
#' @export
#'
get_gene_cna_calls <- function(
    cna_tsv_fname = NULL,
    gencode_xref = NULL,
    ignore_neutral_unknown = TRUE,
    protein_coding_only = TRUE){

  # Determine gene copy number state based on ASCAT3 copy number estimates
  #
  # Homozygous deletion: copy number = 0
  # Loss: copy number > 0, but < sample ploidy
  # Neutral: copy number = sample ploidy
  # Gain: copy number > sample ploidy, but < 2 × sample ploidy
  # Amplification: copy number ≥ 2 × sample ploidy

  cna_calls <- data.frame()
  if(file.exists(cna_tsv_fname)){
    cna_df <- readr::read_tsv(
      file = cna_tsv_fname,
      guess_max = 100000,
      show_col_types = FALSE
    )

    ## Use mode of copy number as sample ploidy
    sample_ploidy <- as.numeric(
      names(which.max(table(
        cna_df$copy_number, useNA = "no"))))

    cna_calls <- cna_df |>
      dplyr::mutate(
        ploidy = sample_ploidy,
        copy_number = as.numeric(copy_number)) |>
      dplyr::mutate(
        mut_status = dplyr::case_when(
          copy_number == 0 ~ "HOMDEL",
          copy_number > 0 & copy_number < sample_ploidy ~ "HEMDEL",
          copy_number == sample_ploidy ~ "NEUTRAL",
          copy_number > sample_ploidy &
            copy_number < (2 * sample_ploidy) ~ "GAIN",
          copy_number >= (2 * sample_ploidy) ~ "AMPL",
          TRUE ~ "UNKNOWN"
        )
      ) |>
      dplyr::mutate(ENSEMBL_GENE_ID = stringr::str_replace(
        .data$gene_id, "\\.[0-9]{1,}$","")) |>
      dplyr::select(
        c("ENSEMBL_GENE_ID",
          "ploidy",
          "copy_number",
          "mut_status")
      ) |>
      dplyr::left_join(
        gencode_xref$all,
        by = "ENSEMBL_GENE_ID"
      ) |>
      dplyr::filter(!is.na(SYMBOL)) |>
      dplyr::rename_with(tolower) |>
      dplyr::distinct() |>
      dplyr::select(
        c("ensembl_gene_id",
          "symbol",
          "entrezgene",
          "biotype",
          "genename",
          "ploidy",
          "copy_number",
          "mut_status")
      )

    if(ignore_neutral_unknown){
      cna_calls <- cna_calls |>
        dplyr::filter(
          mut_status != "NEUTRAL" &
            mut_status != "UNKNOWN"
        )
    }

    if(protein_coding_only){
      cna_calls <- cna_calls |>
        dplyr::filter(
          .data$biotype == "protein_coding"
        )
    }

    cna_calls <- cna_calls |>
      dplyr::select(-c("biotype","genename"))
  }

  return(cna_calls)

}


#' Read and process a single RNAseq sample file
#'
#' @param fname File name of the RNAseq sample file
#' @param gencode_xref Data frame mapping GENCODE gene identifiers
#' @return Data frame with processed RNAseq data
#' @export
#'
read_sample_rnaseq_data <- function(
    fname,
    gencode_xref = NULL){

  rnaseq_data <- readr::read_tsv(
    file = fname,
    show_col_types = F,
    comment = "#") |>
    dplyr::select(
      c("gene_id","tpm_unstranded")) |>
    dplyr::rename(TPM = "tpm_unstranded") |>
    dplyr::mutate(
      ENSEMBL_GENE_ID = stringr::str_replace(
        .data$gene_id,"\\.[0-9]{1,}$","")) |>
    dplyr::select(ENSEMBL_GENE_ID, TPM) |>
    dplyr::filter(!is.na(TPM)) |>
    dplyr::inner_join(
        gencode_xref$all,
      by = "ENSEMBL_GENE_ID"
    ) |>
    dplyr::distinct() |>
    as.data.frame()

  return(rnaseq_data)

}
