
## GISTIC CNA ENCODING
## -2: homozygous deletion (HOMDEL)
## -1: hemizygous deletion (HEMDEL)
##  0: neutral/no change (NONE)
##  1: gain (GAIN)
##  2: high-level amplification (AMPL)

get_gistic_calls <- function(
  gdc_projects = NA,
  tcga_clinical_info = NA,
  tcga_release = NA,
  data_raw_dir = NA,
  gOncoX = NA,
  overwrite = F,
  clear_cache = T,
  output_dir = NA){

  options(timeout = 400000000)

  if(!dir.exists(
    file.path(output_dir, tcga_release))){
    dir.create(
      file.path(output_dir, tcga_release))
  }

  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "cna"))){
    dir.create(
      file.path(output_dir, tcga_release, "cna"))
  }

  k <- 0
  all_cna_calls <- data.frame()

  for(t in unique(gdc_projects$tumor)){
    cat(t,sep="\n")

    #rds_fname <- paste0('output/cna/gistic_',toupper(t),'_20190523.rds')
    rds_fname <- file.path(
      output_dir,
      tcga_release,
      "cna", paste0("tcga_cna_gistic2_", toupper(t),'.rds')
    )

    cna_calls_final <- NULL
    if(!file.exists(rds_fname) | overwrite == TRUE){

      ## OLD:
      #gistic_calls <- TCGAbiolinks::getGistic(t)

      lastAnalyseDate <- RTCGAToolbox::getFirehoseAnalyzeDates(1)

      if(!dir.exists(
        file.path(
          data_raw_dir, "TCGAFirehose"))){
        dir.create(
          file.path(data_raw_dir, "TCGAFirehose"))
      }

      gistic <- RTCGAToolbox::getFirehoseData(
        t, gistic2_Date = lastAnalyseDate,
        GISTIC = TRUE,
        destdir = file.path(
          data_raw_dir, "TCGAFirehose")
      )

      gistic_calls_raw <- list()

      # get GISTIC results
      gistic_calls_raw[['all']] <-
        RTCGAToolbox::getData(
          gistic, type = "GISTIC",
          platform = "AllByGene")

      gistic_calls_raw[['thresholded']] <-
        RTCGAToolbox::getData(
          gistic, type = "GISTIC",
          platform = "ThresholdedByGene")

      cna_df <- list()

      for(e in c('all','thresholded')){

        entrez2sym_gistic <-
          data.frame('entrezgene_gistic' = as.integer(gistic_calls_raw[[e]][,"Locus.ID"]),
                     'symbol_gistic' = as.character(gistic_calls_raw[[e]][,"Gene.Symbol"]),
                     stringsAsFactors = F) |>
          dplyr::left_join(
            dplyr::select(gOncoX$alias, value, entrezgene),
            by = c("symbol_gistic" = "value")
          ) |>
          dplyr::mutate(entrezgene = dplyr::if_else(
            is.na(entrezgene) & entrezgene_gistic > 0,
            as.integer(entrezgene_gistic),
            as.integer(entrezgene)
          )) |>
          dplyr::select(
            entrezgene_gistic,
            entrezgene,
          ) |>
          dplyr::filter(
            !is.na(entrezgene)) |>
          dplyr::left_join(
            dplyr::select(
              gOncoX$basic$records,
              entrezgene,
              gene_biotype,
              symbol),
            by = "entrezgene"
          ) |>
          dplyr::filter(
            !is.na(gene_biotype) &
              gene_biotype == "protein-coding"
          )

        rownames(gistic_calls_raw[[e]]) <- gistic_calls_raw[[e]][,"Locus.ID"]
        for(m in c("Gene.Symbol","Locus.ID","Cytoband")){
          gistic_calls_raw[[e]][,m] <- NULL
        }

        colnames(gistic_calls_raw[[e]]) <- stringr::str_replace_all(
          colnames(gistic_calls_raw[[e]]),"\\.","-")

        cna_calls <- dfrtopics::gather_matrix(
          as.matrix(gistic_calls_raw[[e]]),
          col_names = c('entrezgene_gistic','tumor_sample_barcode','cna_code'))
        cna_calls$n_samples <- length(unique(cna_calls$tumor_sample_barcode))
        cna_calls$entrezgene_gistic <- as.integer(cna_calls$entrezgene_gistic)
        cna_calls <- cna_calls |>
          dplyr::left_join(entrez2sym_gistic, by = "entrezgene_gistic") |>
          dplyr::select(-entrezgene_gistic) |>
          dplyr::distinct() |>
          dplyr::filter(stringr::str_detect(
            tumor_sample_barcode,"TCGA-([:alnum:]){2}-([:alnum:]){4}-0[0-9][A-Z]")) |>
          dplyr::mutate(tumor_sample_barcode = stringr::str_extract(
            tumor_sample_barcode,"TCGA-([:alnum:]){2}-([:alnum:]){4}-0[0-9][A-Z]")) |>
          dplyr::mutate(bcr_patient_barcode = stringr::str_extract(
            tumor_sample_barcode,"TCGA-([:alnum:]){2}-([:alnum:]){4}")) |>
          dplyr::mutate(sample_type = dplyr::if_else(
            stringr::str_detect(
              tumor_sample_barcode,
              "-0(1|5)[A-Z]$"),"Solid Tumor - Primary",
            as.character(NA))) |>
          dplyr::mutate(sample_type = dplyr::if_else(
            stringr::str_detect(
              tumor_sample_barcode,
              "-02[A-Z]$"),"Solid Tumor - Recurrent",
            as.character(sample_type))) |>
          dplyr::mutate(sample_type = dplyr::if_else(
            stringr::str_detect(
              tumor_sample_barcode,"-0(6|7)[A-Z]$"),
            "Metastatic",
            as.character(sample_type))) |>
          dplyr::mutate(sample_type = dplyr::if_else(
            stringr::str_detect(
              tumor_sample_barcode,
              "-0(3|4|9)[A-Z]$"),
            "Blood-Derived Cancer",
            as.character(sample_type)))

        if(e == 'thresholded'){
          cna_calls <- cna_calls |>
            dplyr::mutate(mut_status = dplyr::case_when(
              cna_code == "2" ~ "AMPL",
              cna_code == "1" ~ "GAIN",
              cna_code == "0" ~ "NA",
              cna_code == "-1" ~ "HEMDEL",
              cna_code == "-2" ~ "HOMDEL",
              TRUE ~ as.character(cna_code))
            )
          cna_calls <- cna_calls |>
            dplyr::filter(cna_code != "0") |>
            #dplyr::filter(mut_status == "HOMDEL" | mut_status == "AMPL") |>
            dplyr::mutate(tumor = t) |>
            dplyr::filter(!is.na(entrezgene)) |>
            dplyr::distinct()
        }
        cna_df[[e]] <- cna_calls
      }

      tmp <- as.data.frame(
        cna_df$all |>
          dplyr::select(tumor_sample_barcode,
                        cna_code, entrezgene,
                        bcr_patient_barcode) |>
          dplyr::rename(cna_signal_raw = cna_code) |>
          dplyr::group_by(
            dplyr::across(c(-cna_signal_raw))) |>
          dplyr::summarise(
            cna_signal_raw = max(as.numeric(cna_signal_raw)),
            .groups = "drop")
      )

      #cna_calls <- cna_df

      cna_calls_final <- cna_df$thresholded |>
        dplyr::inner_join(
          tmp,
          by = c("entrezgene",
                 "tumor_sample_barcode",
                 "bcr_patient_barcode"),
          relationship = "many-to-many") |>
        dplyr::left_join(
          dplyr::select(tcga_clinical_info$slim, bcr_patient_barcode,
                        primary_site, primary_diagnosis,
                        primary_diagnosis_very_simplified),
          by = "bcr_patient_barcode", relationship = "many-to-many") |>
        dplyr::select(
          bcr_patient_barcode,
          tumor_sample_barcode,
          tumor, sample_type,
          primary_site,
          primary_diagnosis,
          primary_diagnosis_very_simplified,
          cna_code,
          mut_status,
          cna_signal_raw,
          dplyr::everything()
        ) |>
        dplyr::mutate(
          lastAnalyseDateGistic = lastAnalyseDate
        )

      saveRDS(cna_calls_final, file = rds_fname)

      if(clear_cache == T){
        system(paste0(
          'rm -f ',
          file.path(
            data_raw_dir,
            "TCGAFirehose",
            paste0("*-",t,"-*"))))
      }
    }
    else{
      cna_calls_final <- readRDS(file = rds_fname)
    }
    all_cna_calls <- dplyr::bind_rows(
      all_cna_calls, cna_calls_final)

  }


  tcga_cna <- all_cna_calls |>
    dplyr::filter(
      !is.na(mut_status) &
        (mut_status == 'HOMDEL' | mut_status == 'AMPL')
    )

  saveRDS(
    tcga_cna,
    file = file.path(
      output_dir,
      tcga_release,
      "cna",
      "tcga_cna_gistic2.ampl_homdel.rds"
    )
  )

  readr::write_tsv(
    tcga_cna,
    file = file.path(
      output_dir,
      tcga_release,
      "cna",
      "tcga_cna_gistic2.ampl_homdel.tsv.gz"
    )
  )

}


#' Map GDC ASCAT3 gene-level copy number files to allele-specific segment files
#'
#' The GDC aliquot UUID in the gene-level file name is shared with the
#' corresponding allele-specific segment file
#' (<project>.<uuid>.ascat3.allelic_specific.seg.txt)
#'
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_projects Character vector of GDC project IDs
#' @param ascat_segment_dir Directory with GDC ASCAT3 allele-specific segment files
#'
#' @export
#'
get_ascat3_file_map <- function(
    data_raw_dir = NULL,
    gdc_projects = NULL,
    ascat_segment_dir = file.path(
      data_raw_dir, "ascat3", "202502", "sample_data")){

  cna_sample_metadata <-
    readr::read_tsv(
      file = file.path(
        data_raw_dir,
        "GDCdata",
        "cna",
        "cna_gdc.metadata.tsv.gz"),
      show_col_types = F)

  fname_info <- purrr::map_dfr(
    gdc_projects,
    ~ tibble::tibble(
      project_id = .x,
      fname = list.files(
        path = file.path(
          data_raw_dir, "GDCdata", "cna", .x
        ),
        pattern = "\\.tsv.gz$",
        full.names = TRUE,
        recursive = TRUE
      )
    )) |>
    dplyr::mutate(
      file_id = basename(dirname(fname)),
      gdc_aliquot_uuid = stringr::str_split_i(
        basename(fname), "\\.", 2),
      seg_fname = file.path(
        ascat_segment_dir,
        paste0(project_id, ".", gdc_aliquot_uuid,
               ".ascat3.allelic_specific.seg.txt"))
    ) |>
    dplyr::mutate(
      seg_fname = dplyr::if_else(
        file.exists(seg_fname), seg_fname, NA_character_))

  cna_sample_metadata |>
    dplyr::filter(
      project_id %in% gdc_projects) |>
    dplyr::inner_join(
      fname_info,
      by = c("project_id", "file_id")) |>
    dplyr::select(
      c("project_id", "file_id", "tumor_sample_barcode",
        "gdc_aliquot_uuid", "fname", "seg_fname"))

}

#' Read GDC ASCAT3 allele-specific copy number segments for a single sample
#'
#' @param seg_fname Path to ASCAT3 allele-specific segment file
#'
#' @export
#'
read_ascat3_segments <- function(seg_fname = NULL){

  readr::read_tsv(
    file = seg_fname,
    col_types = "ccnnnnn",
    show_col_types = F) |>
    dplyr::filter(
      !is.na(Major_Copy_Number) &
        !is.na(Minor_Copy_Number)) |>
    dplyr::select(
      chromosome = Chromosome,
      start = Start,
      end = End,
      copy_number = Copy_Number,
      nMajor = Major_Copy_Number,
      nMinor = Minor_Copy_Number)

}

#' Sample ploidy from ASCAT3 segments, defined (as in ASCAT) as the
#' length-weighted mean total copy number
#'
#' @param segments data frame with segments (from read_ascat3_segments())
#'
#' @export
#'
get_ascat3_sample_ploidy <- function(segments = NULL){

  if(is.null(segments) || nrow(segments) == 0){
    return(NA_real_)
  }
  seg_length <- segments$end - segments$start + 1
  round(
    sum(seg_length * (segments$nMajor + segments$nMinor)) /
      sum(seg_length), 4)

}

#' Get ASCAT purity/ploidy estimates for TCGA SNP6 samples, as released by
#' the ASCAT developers (one representative sample per case, QC 'Pass' or
#' 'Likely normal'):
#' https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData/TCGA_SNP6_hg38
#'
#' @param data_raw_dir Directory containing raw data
#'
#' @export
#'
get_ascat_tcga_release <- function(data_raw_dir = NULL){

  ascat_summary_fname <- file.path(
    data_raw_dir, "ascat3",
    "summary.ascatv3TCGA.penalty70.hg38.tsv")

  if(!file.exists(ascat_summary_fname)){
    download.file(
      url = paste0(
        "https://raw.githubusercontent.com/VanLoo-lab/ascat/master/",
        "ReleasedData/TCGA_SNP6_hg38/summary.ascatv3TCGA.penalty70.hg38.tsv"),
      destfile = ascat_summary_fname)
  }

  readr::read_tsv(
    file = ascat_summary_fname,
    show_col_types = F) |>
    dplyr::mutate(
      tumor_sample_barcode = substr(barcodeTumour, 1, 16)) |>
    dplyr::select(
      tumor_sample_barcode,
      ascat_tumor_aliquot = barcodeTumour,
      ascat_purity = purity,
      ascat_ploidy = ploidy,
      ascat_goodness_of_fit = goodness_of_fit,
      ascat_wgd = WGD,
      ascat_gi = GI,
      ascat_loh = LOH,
      ascat_qc = QC)

}

#' Get or make ASCAT3 allele-specific copy number segments and sample
#' purity/ploidy for a TCGA project
#'
#' Segments (all samples with GDC ASCAT3 data) are written to
#' 'tcga_cna_segments_ASCAT3_{project}.rds', sample-level estimates to
#' 'tcga_cna_purity_ploidy_ASCAT3_{project}.rds'. Sample ploidy is computed
#' from the GDC segments (ploidy), purity/ploidy/WGD/QC from the ASCAT
#' TCGA release are added where available (ascat_*).
#'
#' @param output_dir Directory to save/load processed CNA data
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_project Character - GDC project ID
#' @param gdc_release TCGA data release version
#' @param overwrite Whether to overwrite existing processed data
#'
#' @export
#'
gdc_tcga_cna_segments <- function(
    output_dir = NULL,
    data_raw_dir = NULL,
    gdc_project = NULL,
    gdc_release = "release45_20251204",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    !is.null(data_raw_dir),
    !is.null(gdc_project)
  )

  project <- stringr::str_replace_all(
    gdc_project, "TCGA-", "")

  segments_fname <- file.path(
    output_dir, gdc_release, "cna",
    glue::glue("tcga_cna_segments_ASCAT3_{project}.rds"))
  purity_ploidy_fname <- file.path(
    output_dir, gdc_release, "cna",
    glue::glue("tcga_cna_purity_ploidy_ASCAT3_{project}.rds"))

  if(file.exists(segments_fname) &
     file.exists(purity_ploidy_fname) &
     overwrite == FALSE){
    return(list(
      'segments' = readRDS(file = segments_fname),
      'purity_ploidy' = readRDS(file = purity_ploidy_fname)))
  }

  if(!dir.exists(file.path(output_dir, gdc_release, "cna"))){
    dir.create(
      file.path(output_dir, gdc_release, "cna"),
      recursive = T)
  }

  cna_files <- get_ascat3_file_map(
    data_raw_dir = data_raw_dir,
    gdc_projects = gdc_project) |>
    dplyr::filter(!is.na(seg_fname))

  segments <- purrr::pmap_dfr(
    cna_files,
    function(seg_fname, tumor_sample_barcode, gdc_aliquot_uuid, ...) {
      read_ascat3_segments(seg_fname = seg_fname) |>
        dplyr::mutate(
          tumor_sample_barcode = tumor_sample_barcode,
          gdc_aliquot_uuid = gdc_aliquot_uuid)
    }) |>
    dplyr::select(
      tumor_sample_barcode, gdc_aliquot_uuid,
      dplyr::everything())

  purity_ploidy <- segments |>
    dplyr::group_by(
      tumor_sample_barcode, gdc_aliquot_uuid) |>
    dplyr::group_modify(
      ~ tibble::tibble(
        ploidy = get_ascat3_sample_ploidy(.x),
        n_segments = nrow(.x))) |>
    dplyr::ungroup() |>
    dplyr::mutate(project_id = gdc_project) |>
    dplyr::left_join(
      get_ascat_tcga_release(data_raw_dir = data_raw_dir),
      by = "tumor_sample_barcode") |>
    dplyr::select(
      project_id, tumor_sample_barcode,
      dplyr::everything())

  saveRDS(segments, file = segments_fname)
  saveRDS(purity_ploidy, file = purity_ploidy_fname)

  return(list(
    'segments' = segments,
    'purity_ploidy' = purity_ploidy))

}


#' Get or make processed TCGA CNA data from GDC
#'
#' Gene-level copy number states are determined relative to the sample
#' ploidy derived from the ASCAT3 allele-specific segments
#' (see gdc_tcga_cna_segments())
#'
#' @param output_dir Directory to save/load processed CNA data
#' @param gencode_xref Data frame with GENCODE mapping gene identifiers
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_project Character - GDC project ID
#' @param gdc_release TCGA data release version
#' @param overwrite Whether to overwrite existing processed data
#'
#' @export
#'
gdc_tcga_cna <- function(
    output_dir = NULL,
    gencode_xref = NULL,
    data_raw_dir = NULL,
    gdc_project = NULL,
    gdc_release = "release45_20251204",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    !is.null(data_raw_dir),
    !is.null(gencode_xref),
    !is.null(gdc_project)
  )
  assertthat::assert_that(
    dir.exists(output_dir),
    dir.exists(data_raw_dir)
  )

  project <- stringr::str_replace_all(
    gdc_project,
    "TCGA-","")

  output_fname <- file.path(
    output_dir,
    gdc_release,
    "cna",
    glue::glue(
      "tcga_cna_ASCAT3_{project}.rds")
  )

  if(file.exists(output_fname) & overwrite == FALSE){
    all_calls <- readRDS(file = output_fname)
    return(all_calls)
  }

  purity_ploidy <- gdc_tcga_cna_segments(
    output_dir = output_dir,
    data_raw_dir = data_raw_dir,
    gdc_project = gdc_project,
    gdc_release = gdc_release,
    overwrite = overwrite)$purity_ploidy

  cna_files <- get_ascat3_file_map(
    data_raw_dir = data_raw_dir,
    gdc_projects = gdc_project) |>
    dplyr::left_join(
      dplyr::select(
        purity_ploidy, gdc_aliquot_uuid, ploidy),
      by = "gdc_aliquot_uuid")

  n_missing_ploidy <- sum(is.na(cna_files$ploidy))
  if(n_missing_ploidy > 0){
    cat(glue::glue(
      "{gdc_project}: no ASCAT3 segments for {n_missing_ploidy} sample(s)",
      " - using mode of gene copy number as ploidy"), sep = "\n")
  }

  # Process all CNA files efficiently
  all_calls <- purrr::pmap_dfr(
    cna_files,
    function(fname, tumor_sample_barcode, ploidy, ...) {

      cna_calls <- get_gene_cna_calls(
        cna_tsv_fname = fname,
        gencode_xref = gencode_xref,
        sample_ploidy = ploidy,
        ignore_neutral_unknown = TRUE,
        protein_coding_only = TRUE
      )

      if (nrow(cna_calls) == 0) return(NULL)

      cna_calls |>
        dplyr::mutate(
          tumor_sample_barcode = tumor_sample_barcode,
        )
    }
  )

  saveRDS(
    all_calls,
    file = output_fname
  )

  return(all_calls)


}
