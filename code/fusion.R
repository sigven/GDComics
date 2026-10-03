#' Map a GDC TCGA project ID to its cBioPortal PanCancer Atlas study ID
#'
#' @param gdc_project Character - GDC project ID (e.g. "TCGA-OV")
#'
pancan_atlas_study_id <- function(gdc_project = NULL){

  code <- tolower(
    stringr::str_replace(gdc_project, "TCGA-", ""))
  ## COAD and READ are combined in one PanCancer Atlas study
  code <- dplyr::if_else(
    code %in% c("coad", "read"), "coadread", code)
  paste0(code, "_tcga_pan_can_atlas_2018")

}

#' Fetch RNA fusion calls (structural variants profile) for a cBioPortal
#' TCGA PanCancer Atlas study, cached as TSV in the raw data directory
#'
#' @param study_id cBioPortal study ID (e.g. "ov_tcga_pan_can_atlas_2018")
#' @param cache_dir Directory for cached raw fusion calls
#' @param overwrite Whether to re-fetch cached data
#'
#' @export
#'
fetch_cbioportal_fusions <- function(
    study_id = NULL,
    cache_dir = NULL,
    overwrite = FALSE){

  cache_fname <- file.path(
    cache_dir, paste0(study_id, ".structural_variants.tsv.gz"))

  if(file.exists(cache_fname) & overwrite == FALSE){
    return(readr::read_tsv(
      cache_fname, guess_max = 100000, show_col_types = F))
  }

  sv_calls <- httr2::request(
    "https://www.cbioportal.org/api/structural-variant/fetch") |>
    httr2::req_body_json(
      list(molecularProfileIds = list(
        paste0(study_id, "_structural_variants")))) |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = TRUE) |>
    tibble::as_tibble()

  if(nrow(sv_calls) > 0){
    readr::write_tsv(sv_calls, file = cache_fname)
  }
  return(sv_calls)

}

#' Lift single-base positions between assemblies
#'
#' @param chrom Character vector of chromosome names (without 'chr')
#' @param pos Numeric vector of positions
#' @param chain rtracklayer Chain object
#'
#' @return numeric vector of lifted positions (NA if unmapped, ambiguous or
#' mapped to a different chromosome)
#'
liftover_positions <- function(chrom = NULL, pos = NULL, chain = NULL){

  idx <- which(!is.na(chrom) & !is.na(pos) & pos > 0)
  lifted_pos <- rep(NA_real_, length(pos))
  if(length(idx) == 0){
    return(lifted_pos)
  }
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", chrom[idx]),
    ranges = IRanges::IRanges(start = pos[idx], width = 1))
  lifted <- rtracklayer::liftOver(gr, chain)
  ## keep unique, same-chromosome mappings only
  ok <- S4Vectors::elementNROWS(lifted) == 1
  lifted_df <- as.data.frame(lifted[ok])
  same_chrom <- as.character(lifted_df$seqnames) ==
    paste0("chr", chrom[idx][ok])
  lifted_pos[idx[ok][same_chrom]] <- lifted_df$start[same_chrom]
  return(lifted_pos)

}

#' Get or make TCGA RNA fusion calls from the PanCancer Atlas
#' (Gao et al., Cell Reports 2018), as distributed by cBioPortal
#'
#' Breakpoints are provided on GRCh38 (verified against GENCODE/Ensembl gene
#' coordinates: >99% of breakpoints fall within their partner gene on GRCh38,
#' vs. ~15% on GRCh37) and lifted to GRCh37. Gene partners are ordered
#' 5' (gene5) - 3' (gene3). Output is written to
#' 'fusion/tcga_fusions.rds' and 'fusion/tcga_fusions.tsv.gz'
#'
#' @param gdc_projects Character vector of GDC project IDs
#' @param tcga_clinical_info Data frame with TCGA clinical data
#' (project_id, bcr_patient_barcode), used to assign projects to patients
#' (COAD/READ share one PanCancer Atlas study)
#' @param data_raw_dir Directory containing raw data
#' @param output_dir Directory to save/load processed data
#' @param gdc_release TCGA data release version
#' @param overwrite Whether to overwrite existing processed data
#' (and re-fetch cached cBioPortal data)
#'
#' @export
#'
get_tcga_fusions <- function(
    gdc_projects = NULL,
    tcga_clinical_info = NULL,
    data_raw_dir = NULL,
    output_dir = NULL,
    gdc_release = "release45_20251204",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(gdc_projects),
    !is.null(tcga_clinical_info),
    !is.null(data_raw_dir),
    !is.null(output_dir)
  )

  fusion_output_dir <- file.path(output_dir, gdc_release, "fusion")
  output_fname <- file.path(fusion_output_dir, "tcga_fusions.rds")

  if(file.exists(output_fname) & overwrite == FALSE){
    return(readRDS(file = output_fname))
  }

  cache_dir <- file.path(data_raw_dir, "cbioportal", "fusion")
  for(d in c(fusion_output_dir, cache_dir, file.path(data_raw_dir, "liftover"))){
    if(!dir.exists(d)){
      dir.create(d, recursive = T)
    }
  }

  chain_fname <- file.path(
    data_raw_dir, "liftover", "hg38ToHg19.over.chain.gz")
  if(!file.exists(chain_fname)){
    download.file(
      url = paste0(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/",
        "liftOver/hg38ToHg19.over.chain.gz"),
      destfile = chain_fname)
  }
  chain_tmp <- tempfile(fileext = ".chain")
  R.utils::gunzip(chain_fname, destname = chain_tmp, remove = FALSE)
  chain <- rtracklayer::import.chain(chain_tmp)
  file.remove(chain_tmp)

  sv_calls <- purrr::map_dfr(
    unique(pancan_atlas_study_id(gdc_projects)),
    function(study_id){
      cat(study_id, sep = "\n")
      fetch_cbioportal_fusions(
        study_id = study_id,
        cache_dir = cache_dir,
        overwrite = overwrite) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    })

  patient_project <- tcga_clinical_info |>
    dplyr::select(c("bcr_patient_barcode", "project_id")) |>
    dplyr::distinct()

  na_values <- c("NA", "", "-1")
  fusions <- sv_calls |>
    dplyr::mutate(dplyr::across(
      dplyr::everything(),
      ~ dplyr::if_else(.x %in% na_values, NA_character_, .x))) |>
    dplyr::transmute(
      bcr_patient_barcode = patientId,
      tumor_sample_barcode = sampleId,
      cbioportal_study = studyId,
      gene5 = site1HugoSymbol,
      gene3 = site2HugoSymbol,
      entrezgene5 = as.integer(site1EntrezGeneId),
      entrezgene3 = as.integer(site2EntrezGeneId),
      chrom5 = site1Chromosome,
      pos5_grch38 = as.numeric(site1Position),
      chrom3 = site2Chromosome,
      pos3_grch38 = as.numeric(site2Position),
      split_reads = as.integer(tumorSplitReadCount),
      spanning_reads = as.integer(tumorPairedEndReadCount),
      frame_effect = site2EffectOnFrame,
      fusion_label = eventInfo,
      sv_status = svStatus) |>
    dplyr::mutate(
      pos5_grch37 = liftover_positions(chrom5, pos5_grch38, chain),
      pos3_grch37 = liftover_positions(chrom3, pos3_grch38, chain)) |>
    dplyr::left_join(
      patient_project, by = "bcr_patient_barcode") |>
    dplyr::filter(
      is.na(project_id) | project_id %in% gdc_projects) |>
    dplyr::select(
      c("project_id", "bcr_patient_barcode",
        "tumor_sample_barcode", "gene5", "gene3"),
      dplyr::everything()) |>
    dplyr::distinct()

  n_unmapped <- sum(is.na(fusions$project_id))
  if(n_unmapped > 0){
    cat(glue::glue(
      "{n_unmapped} fusion calls from patients without GDC clinical data",
      " (project_id = NA)"), sep = "\n")
  }

  saveRDS(fusions, file = output_fname)
  readr::write_tsv(
    fusions,
    file = file.path(fusion_output_dir, "tcga_fusions.tsv.gz"))

  return(fusions)

}
