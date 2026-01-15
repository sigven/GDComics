
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


#' Get or make processed TCGA CNA data from GDC
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

  gdc_projects <- c(gdc_project)
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
      file_id = stringr::str_split(
        fname, "/", simplify = TRUE)[,12]
    )

  # Join metadata with filenames
  cna_files <- cna_sample_metadata |>
    dplyr::filter(
      project_id %in% gdc_projects) |>
    dplyr::inner_join(
      fname_info,
      by = c("project_id", "file_id"))

  # Process all CNA files efficiently
  all_calls <- purrr::pmap_dfr(
    cna_files,
    function(fname, tumor_sample_barcode, ...) {

      cna_calls <- get_gene_cna_calls(
        cna_tsv_fname = fname,
        gencode_xref = gencode_xref,
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













