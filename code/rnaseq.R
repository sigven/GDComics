
# get_bulk_rnaseq <- function(
#     gdc_projects = NA,
#     tcga_clinical_info = NA,
#     tcga_release = NA,
#     data_raw_dir = NA,
#     overwrite = FALSE,
#     gOncoX = NA,
#     output_dir = NA){
#
#   if(!dir.exists(
#     file.path(output_dir, tcga_release))){
#     dir.create(
#       file.path(output_dir, tcga_release))
#   }
#
#   if(!dir.exists(
#     file.path(
#       output_dir, tcga_release, "rnaseq"))){
#     dir.create(
#       file.path(output_dir, tcga_release, "rnaseq"))
#   }
#
#   ## NOTE: As of GDC release v32, GENCODE v36 is used for annotation of RNA-seq genes/transcripts
#   ## https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-320
#   update_ensembl_gene <- F
#   gencode_xref_fname <-
#     file.path(data_raw_dir, "gencode_v36_xref.rds")
#   gencode_xref <- list()
#   if(!file.exists(gencode_xref_fname) | update_ensembl_gene == T){
#
#     gencode_gtf_fname <- file.path(
#       data_raw_dir,
#       "gencode",
#       "gencode.v36.annotation.gtf.gz")
#
#     gencode_xref_df <- as.data.frame(
#       read_gtf(path = gencode_gtf_fname) |>
#       dplyr::select(gene_id, hgnc_id) |>
#       dplyr::mutate(ensembl_gene_id = stringr::str_replace(
#         gene_id,"\\.[0-9]{1,}","")) |>
#       dplyr::select(ensembl_gene_id, hgnc_id) |>
#       dplyr::filter(!is.na(hgnc_id)) |>
#       dplyr::filter(!stringr::str_detect(
#         ensembl_gene_id,"_PAR_Y")) |>
#       dplyr::mutate(hgnc_id = stringr::str_replace(
#         hgnc_id, "HGNC:",""
#       )) |>
#       ## remove ensembl gene id's that are pointing
#       ## to the same gene symbol as other
#       ## ensembl identifiers
#       dplyr::filter(ensembl_gene_id != "ENSG00000145075") |>
#       dplyr::filter(ensembl_gene_id != "ENSG00000285437") |>
#       dplyr::filter(ensembl_gene_id != "ENSG00000286065") |>
#       dplyr::distinct() |>
#
#       ## join with geneOncoX basic records to get
#       ## gene symbols, entrezgene, gene biotype info
#       ## etc. that are up-to-date
#       dplyr::inner_join(
#         dplyr::select(gOncoX$basic$records,
#                       symbol, entrezgene,
#                       name,
#                       hgnc_id, gene_biotype),
#         by = "hgnc_id"
#       ) |>
#         dplyr::rename(genename = "name") |>
#       dplyr::select(
#         c("ensembl_gene_id",
#           "symbol",
#           "entrezgene",
#           "genename",
#           #"hgnc_id",
#           "gene_biotype")
#       )
#     )
#
#     dups <- as.data.frame(
#       dplyr::group_by(gencode_xref_df, entrezgene) |>
#       dplyr::summarise(n = dplyr::n()) |>
#       dplyr::filter(n > 1)
#     )
#
#     gencode_xref[['all']] <- gencode_xref_df |>
#       dplyr::anti_join(dups, by = "entrezgene") |>
#       dplyr::arrange(.data$entrezgene) |>
#       dplyr::distinct()
#
#     gencode_xref[['protein_coding']] <-
#       gencode_xref[['all']] |>
#       dplyr::filter(.data$gene_biotype == "protein_coding")
#
#     saveRDS(gencode_xref, file = gencode_xref_fname)
#   }else{
#     gencode_xref <- readRDS(file = gencode_xref_fname)
#   }
#
#   for(i in 1:nrow(gdc_projects)){
#
#     rnaseq_query <- NULL
#     rnaseq_data <- list()
#
#     project_code <- gdc_projects[i,]$project_id
#     tumor <- gdc_projects[i,]$tumor
#
#     cat(i, tumor, sep=" - ")
#     cat('\n')
#
#     downloads_raw_gdc <- file.path(
#       data_raw_dir, "GDCdata",
#       project_code, "harmonized",
#       "Transcriptome_Profiling",
#       "Gene_Expression_Quantification"
#     )
#
#     rds_fname <-
#       file.path(
#         output_dir, tcga_release, "rnaseq",
#         paste0('rnaseq_',tumor,'.rds')
#       )
#
#     if(file.exists(rds_fname)){
#       next
#     }
#
#     rnaseq_data <- list()
#     workflow_type <- "STAR - Counts"
#     rnaseq_query <-
#       TCGAbiolinks::GDCquery(
#         project = project_code,
#         data.category = "Transcriptome Profiling",
#         data.type = "Gene Expression Quantification",
#         workflow.type = workflow_type
#       )
#
#     TCGAbiolinks::GDCdownload(
#       rnaseq_query,
#       directory = file.path(
#         data_raw_dir, "GDCdata"),
#       files.per.chunk = 20)
#
#     rnaseq_data[['tcga_release']] <- tcga_release
#     rnaseq_data[['tcga_project']] <- tumor
#     rnaseq_data[['sample_metadata']] <- data.frame()
#     rnaseq_data[['tcga_biolinks']] <- list()
#     rnaseq_data[['tcga_biolinks']][['data.category']] <-
#       'Transcriptome Profiling'
#     rnaseq_data[['tcga_biolinks']][['data.type']] <-
#       'Gene Expression Quantification'
#     rnaseq_data[['tcga_biolinks']][['workflow.type']] <-
#       workflow_type
#     rnaseq_data[['tcga_biolinks']][['version']] <-
#       Biobase::package.version("TCGAbiolinks")
#     rnaseq_raw_SE <-
#       TCGAbiolinks::GDCprepare(
#         rnaseq_query, directory =
#           file.path(data_raw_dir, "GDCdata"))
#
#
#     #system(paste0('rm -rf ', downloads_raw_gdc, "/*"))
#     rnaseq_data[['matrix']] <- list()
#     rnaseq_data[['matrix']][['counts']] <- list()
#     rnaseq_data[['matrix']][['tpm']] <- list()
#
#     ## Get counts and TPM data matrices from SummarizedExperiment object
#     rnaseq_data[['matrix']][['tpm']][['all']] <-
#       SummarizedExperiment::assay(rnaseq_raw_SE, "tpm_unstrand")
#
#     rnaseq_data[['matrix']][['counts']][['all']] <-
#       SummarizedExperiment::assay(rnaseq_raw_SE)
#
#
#     ## Strip away barcode information related to Portion/Plate/Center
#     ## https://www.biostars.org/p/313063/
#     ## Example: TCGA-93-A4JP-01A-11H-A24S-13 --> TCGA-93-A4JP-01A
#     ##          "Barcode" --> "Sample"
#     ##
#     ## Strip away Ensembl gene version, e.g.:
#     ## ENSG00000000003.15 --> ENSG00000000003
#     for(c in c('tpm', 'counts')){
#       colnames(rnaseq_data[['matrix']][[c]][['all']]) <-
#         stringr::str_replace(colnames(rnaseq_data[['matrix']][[c]][['all']]),
#                              "-[0-9]{2}[A-Z]{1}-[A-Z|0-9]{4}-\\S+$",
#                              "")
#       rownames(rnaseq_data[['matrix']][[c]][['all']]) <-
#         stringr::str_replace(rownames(rnaseq_data[['matrix']][[c]][['all']]),
#                              "\\.[0-9]{1,}$",
#                              "")
#
#     }
#
#     rnaseq_data[['sample_metadata']] <-
#       as.data.frame(SummarizedExperiment::colData(rnaseq_raw_SE)) |>
#       dplyr::select(barcode, bcr_patient_barcode, sample_type,
#                     name, primary_diagnosis,
#                     site_of_resection_or_biopsy)
#
#     rownames(rnaseq_data[['sample_metadata']]) <- NULL
#
#     sample_collection <- list()
#     sample_collection[['tumor']] <- rnaseq_data$sample_metadata |>
#       dplyr::filter(
#         stringr::str_detect(.data$bcr_patient_barcode, "-(0[1-9][A-Z])$"))
#     sample_collection[['control']] <- rnaseq_data$sample_metadata |>
#       dplyr::filter(
#         stringr::str_detect(.data$bcr_patient_barcode, "-(1[0-4][A-Z])$"))
#
#     ## Assign data from tumor and control/normal samples in dedicated TPM/count matrices
#     for(exp_measure in c('tpm','counts')){
#       for(c in c('tumor', 'control')){
#
#         ## Filter complete sample set against tumor/control sample set
#         rnaseq_data[['matrix']][[exp_measure]][[c]] <- as.matrix(
#           rnaseq_data[['matrix']][[exp_measure]][['all']][, colnames(rnaseq_data[['matrix']][[exp_measure]][['all']]) %in%
#                                                             sample_collection[[c]]$bcr_patient_barcode]
#         )
#
#
#         if(!is.null(rnaseq_data[['matrix']][[exp_measure]][[c]])){
#           if(NCOL(rnaseq_data[['matrix']][[exp_measure]][[c]]) > 1){
#
#             ## Exclude duplicated sample names
#             rnaseq_data[['matrix']][[exp_measure]][[c]] <-
#               rnaseq_data[['matrix']][[exp_measure]][[c]][, !duplicated(
#                 colnames(rnaseq_data[['matrix']][[exp_measure]][[c]]), fromLast = TRUE)]
#
#             df <- data.frame(
#               'ensembl_gene_id' = rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]),
#               stringsAsFactors = F) |>
#               dplyr::left_join(gencode_xref, by = "ensembl_gene_id")
#
#
#             rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]) <- df$symbol
#             rnaseq_data[['matrix']][[exp_measure]][[c]] <-
#               as.matrix(rnaseq_data[['matrix']][[exp_measure]][[c]][!is.na(rownames(rnaseq_data[['matrix']][[exp_measure]][[c]])), ])
#             rnaseq_data[['matrix']][[exp_measure]][[c]] <-
#               as.matrix(rnaseq_data[['matrix']][[exp_measure]][[c]][rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]) != "NA", ])
#           }
#         }
#       }
#     }
#
#     ## remove matrices with all samples (tumor + normal created)
#     rnaseq_data$matrix$counts$all <- NULL
#     rnaseq_data$matrix$tpm$all <- NULL
#
#     saveRDS(rnaseq_data, rds_fname)
#
#   }
# }

#' Get or make processed (TPM) bulk TCGA RNAseq data
#' from GDC
#'
#' @param output_dir Directory to save/load processed RNAseq data
#' @param gencode_xref Data frame with GENCODE gene identifiers
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_project Character - GDC/TCGA project ID
#' @param gdc_release TCGA data release version
#' @param overwrite Whether to overwrite existing processed data
#'
#' @export
#'
gdc_tcga_rnaseq <- function(
    output_dir = NULL,
    gencode_xref = NULL,
    data_raw_dir = NULL,
    gdc_project = NULL,
    gdc_release = "release45_20251204",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    !is.null(data_raw_dir),
    !is.null(gene_xref),
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
    "rnaseq",
    glue::glue(
      "tcga_rnaseq_TPM_{project}.rds")
  )

  if(file.exists(output_fname) & overwrite == FALSE){
    expr_df <- readRDS(file = output_fname)
    return(expr_df)
  }

  fname_info <- purrr::map_dfr(
    gdc_projects,
    ~ tibble::tibble(
      project_id = .x,
      fname = list.files(
        path = file.path(
          data_raw_dir, "GDCdata",
          "rnaseq", .x
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

  rnaseq_sample_metadata <-
    readr::read_tsv(
      file = file.path(
        data_raw_dir,
        "GDCdata",
        "rnaseq",
        "rnaseq_gdc.metadata.tsv.gz"),
      show_col_types = F) |>

    ## Keep only relevant columns
    dplyr::arrange(
      bcr_patient_barcode,
      tumor_sample_barcode, id,
      file_id) |>
    ## Keep only one entry per sample/file
    dplyr::group_by(
      dplyr::across(-c("id","file_id"))) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
  dplyr::filter(project_id == gdc_project)

  # Join metadata with filenames
  rnaseq_files <- rnaseq_sample_metadata |>
    dplyr::filter(
      project_id %in% gdc_projects) |>
    dplyr::inner_join(
      fname_info,
      by = c("project_id", "file_id")) |>
    dplyr::distinct()

  # Process all RNAseq files efficiently
  all_calls <- purrr::pmap_dfr(
    rnaseq_files,
    function(fname, tumor_sample_barcode, ...) {

      rnaseq_calls <- read_sample_rnaseq_data(
        fname = fname,
        gencode_xref = gencode_xref
      )

      if (nrow(rnaseq_calls) == 0) return(NULL)

      rnaseq_calls |>
        dplyr::mutate(
          tumor_sample_barcode = tumor_sample_barcode
        )
    }
  )

  expr_df <- as.data.frame(all_calls |>
    tidyr::pivot_wider(
      names_from  = tumor_sample_barcode,
      values_from = TPM)) |>
    dplyr::arrange(.data$ENTREZGENE)

  saveRDS(
    expr_df,
    file = output_fname
  )
  return(expr_df)

}

#' Get or make global UMAP embedding of TCGA RNAseq data
#' from GDC
#'
#' @param gdc_release TCGA data release version
#' @param output_dir Directory to save/load processed RNAseq data
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_projects Character vector - GDC/TCGA project IDs
#' @param overwrite Whether to overwrite existing processed data
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param metric UMAP metric parameter
#' @export
#'
#'
global_umap_tcga <- function(
    gdc_release = "release45_20251204",
    output_dir = NULL,
    data_raw_dir = NULL,
    gdc_projects = NULL,
    overwrite = FALSE,
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "cosine"){

  assertthat::assert_that(
    !is.null(output_dir),
    !is.null(data_raw_dir),
    !is.null(gdc_projects)
  )
  assertthat::assert_that(
    dir.exists(output_dir),
    dir.exists(data_raw_dir)
  )

  output_fname <- file.path(
    output_dir,
    gdc_release,
    "rnaseq",
    "tcga_rnaseq_global_umap.rds"
  )

  if(file.exists(output_fname) & overwrite == FALSE){
    umap_result <- readRDS(file = output_fname)
    return(umap_result)
  }

  tcga_clinical <- gdc_tcga_clinical(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir)

  all_expr <- data.frame()

  for(gdc_project in gdc_projects){
    cat("Processing ", gdc_project, "\n")
    project <-
      stringr::str_replace_all(
        gdc_project, 'TCGA-','')
    fname <-
      file.path(
        output_dir,
        gdc_release,
        "rnaseq",
        glue::glue(
          "tcga_rnaseq_TPM_{project}.rds")
      )
    expr_df <- readRDS(file = fname)
    if(NROW(all_expr) == 0){
      all_expr <- expr_df
    }else{
      for(e in c("ENSEMBL_GENE_ID",
                 "SYMBOL",
                 "GENENAME",
                 "BIOTYPE")){
        expr_df[[e]] <- NULL
      }
      all_expr <- all_expr |>
        dplyr::left_join(
          expr_df, by = "ENTREZGENE"
        )
    }

  }

  gene_ids <- all_expr$ENSEMBL_GENE_ID

  all_expr$GENENAME <- NULL
  all_expr$BIOTYPE <- NULL
  all_expr$SYMBOL <- NULL
  all_expr$ENTREZGENE <- NULL
  all_expr$ENSEMBL_GENE_ID <- NULL

  rnaseq_mat <- as.matrix(all_expr)
  storage.mode(rnaseq_mat) <- "numeric"
  rownames(rnaseq_mat) <- gene_ids
  sample_ids <- colnames(rnaseq_mat)

  set.seed(1234)


  # Remove duplicated genes
  rnaseq_mat <- rnaseq_mat[
    !duplicated(rownames(rnaseq_mat)), ]

  # Remove duplicated samples (if needed)
  rnaseq_mat <- rnaseq_mat[
    , !duplicated(colnames(rnaseq_mat))]

  # Transpose: samples x genes
  expr <- t(rnaseq_mat)

  # Select top variable genes
  vars <- matrixStats::colVars(expr)
  n_top <- min(7000, ncol(expr))
  expr_var <- expr[, order(vars, decreasing = TRUE)[seq_len(n_top)]]

  # Log-transform and scale
  expr_var <- log2(expr_var + 0.001)
  expr_var <- scale(expr_var)

  ## store UMAP training settings
  ## for future reference/projection of new samples

  ## 1. store gene set used for UMAP training
  umap_train_settings <- list()
  umap_train_settings[['gene_set']] <-
    colnames(expr_var)

  ## 2. save scaling attributes used for UMAP training data
  umap_train_settings[['scale_center']] <-
    attr(expr_var, "scaled:center")
  umap_train_settings[['scale_scale']] <-
    attr(expr_var, "scaled:scale")


  umap_res <- uwot::umap(
    expr_var,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    scale = FALSE,
    n_components = 2,
    verbose = FALSE,
    ret_model = TRUE,
    n_threads = 2
  )

  umap_df <- as.data.frame(umap_res$embedding)
  colnames(umap_df) <-
    c("UMAP1", "UMAP2")
  umap_df$tumor_sample_barcode <-
    rownames(umap_df)
  umap_df <- umap_df |>
    dplyr::inner_join(
      dplyr::select(
        tcga_clinical,
        c("tumor_sample_barcode",
          "sample_type",
          "primary_site",
          "primary_diagnosis_simplified",
          "tumor_stage_TNM",
          "tumor",
          "project_id")
      ), by = "tumor_sample_barcode"
    ) |>
    dplyr::mutate(tumor_type = primary_site) |>

    ## ignore tumor types with few samples
    dplyr::filter(
      tumor_type != "Bone" &
        tumor_type != "Other/Unknown"
    ) |>
    dplyr::mutate(
      hover_text = glue::glue(
        "<b>Sample type:</b> {sample_type}<br>",
        "<b>Primary tumor site:</b> {primary_site}<br>",
        "<b>Tumor:</b> {project_id}<br>",
        "<b>Primary diagnosis:</b> {primary_diagnosis_simplified}<br>"
      )
    ) |>
    dplyr::distinct()


  base_cols <- ggsci::pal_jco("default")(10)
  cancer_types <- sort(unique(umap_df$tumor_type))
  extended_cols <- grDevices::colorRampPalette(
    base_cols)(length(cancer_types))


  umap_ggplot <- ggplot2::ggplot(
    umap_df, ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      color = tumor_type,
      text = hover_text)) +
    ggplot2::theme_classic(
      base_family = "Helvetica") +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size=10)
    ) +
    ggplot2::scale_color_manual(values = extended_cols)

  umap_plotly <- plotly::ggplotly(
    umap_ggplot,
    tooltip = "text") |>
    plotly::layout(
      legend = list(
        title = list(text = ""),
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.2
      ),
      xaxis = list(
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        showline = FALSE
      ),
      yaxis = list(
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        showline = FALSE
      )
    )


  umap_result <- list('model' = umap_res,
                      'training_settings' = umap_train_settings,
                      'plot_ggplot2' = umap_ggplot,
                      'plot_plotly' = umap_plotly,
                      'embedding' = umap_df)

  saveRDS(umap_result, file = output_fname)

  return(umap_result)

}

