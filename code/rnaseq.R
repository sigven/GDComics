
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
#' @param include_impress_test Whether to include IMPRESS test samples in UMAP
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param spread UMAP spread parameter
#' @param metric UMAP metric parameter
#' @export
#'
#'
#' global_umap_tcga <- function(
#'     gdc_release = "release45_20251204",
#'     output_dir = NULL,
#'     data_raw_dir = NULL,
#'     gdc_projects = NULL,
#'     overwrite = FALSE,
#'     include_impress_test = FALSE,
#'     n_neighbors = 15,
#'     min_dist = 0.3,
#'     spread = 0.7,
#'     metric = "cosine"){
#'
#'   assertthat::assert_that(
#'     !is.null(output_dir),
#'     !is.null(data_raw_dir),
#'     !is.null(gdc_projects)
#'   )
#'   assertthat::assert_that(
#'     dir.exists(output_dir),
#'     dir.exists(data_raw_dir)
#'   )
#'
#'   output_fname <- file.path(
#'     output_dir,
#'     gdc_release,
#'     "rnaseq",
#'     "tcga_rnaseq_umap.rds"
#'   )
#'
#'   output_model_fname <- file.path(
#'     output_dir,
#'     gdc_release,
#'     "rnaseq",
#'     "tcga_rnaseq_umap_model.uwot"
#'   )
#'
#'
#'   if(file.exists(output_fname) &
#'      file.exists(output_model_fname) &
#'      overwrite == FALSE){
#'
#'     umap_result <- readRDS(file = output_fname)
#'     umap_result$model <- uwot::load_uwot(file = output_model_fname)
#'     return(umap_result)
#'   }
#'
#'   tcga_clinical <- gdc_clinical(
#'     gdc_release = gdc_release,
#'     overwrite = F,
#'     output_dir = output_dir)
#'
#'   all_expr <- data.frame()
#'
#'   for(gdc_project in gdc_projects){
#'     cat("Processing ", gdc_project, "\n")
#'     project <-
#'       stringr::str_replace_all(
#'         gdc_project, 'TCGA-','')
#'     fname <-
#'       file.path(
#'         output_dir,
#'         gdc_release,
#'         "rnaseq",
#'         glue::glue(
#'           "tcga_rnaseq_TPM_{project}.rds")
#'       )
#'     expr_df <- readRDS(file = fname)
#'     if(NROW(all_expr) == 0){
#'       all_expr <- expr_df
#'     }else{
#'       for(e in c("ENSEMBL_GENE_ID",
#'                  "SYMBOL",
#'                  "GENENAME",
#'                  "BIOTYPE")){
#'         expr_df[[e]] <- NULL
#'       }
#'       all_expr <- all_expr |>
#'         dplyr::left_join(
#'           expr_df, by = "ENTREZGENE"
#'         )
#'     }
#'
#'   }
#'
#'   if(include_impress_test == TRUE){
#'     impress_sample_data <- readRDS(
#'       file = file.path(
#'         data_raw_dir,
#'         "impress",
#'         "impress_test_exprdata.rds"
#'       ))
#'
#'     exprdata_impress <- impress_sample_data$data |>
#'       dplyr::select(
#'         -c("GENENAME",
#'            "BIOTYPE",
#'            "SYMBOL",
#'            "ENSEMBL_GENE_ID")
#'       )
#'
#'     all_expr <- all_expr |>
#'       dplyr::left_join(
#'         exprdata_impress,
#'         by = "ENTREZGENE"
#'       )
#'
#'     impress_clinical <- impress_sample_data$metadata |>
#'       dplyr::select(
#'         -c("pcgr_tsite_code","primary_site")
#'       ) |>
#'       dplyr::mutate(
#'         project_id = "IMPRESS"
#'       ) |>
#'       dplyr::rename(
#'         tumor_sample_barcode = sample_id,
#'         primary_site = primary_site2
#'       ) |>
#'       dplyr::select(
#'         c("tumor_sample_barcode",
#'           "sample_type",
#'           "primary_site",
#'           "project_id")
#'       ) |>
#'       dplyr::distinct()
#'
#'
#'     tcga_clinical <- tcga_clinical |>
#'       dplyr::bind_rows(impress_clinical)
#'   }
#'
#'   gene_ids <- all_expr$ENSEMBL_GENE_ID
#'
#'   all_expr$GENENAME <- NULL
#'   all_expr$BIOTYPE <- NULL
#'   all_expr$SYMBOL <- NULL
#'   all_expr$ENTREZGENE <- NULL
#'   all_expr$ENSEMBL_GENE_ID <- NULL
#'
#'   rnaseq_mat <- as.matrix(all_expr)
#'   storage.mode(rnaseq_mat) <- "numeric"
#'   rownames(rnaseq_mat) <- gene_ids
#'   sample_ids <- colnames(rnaseq_mat)
#'
#'   set.seed(1234)
#'
#'
#'   # Remove duplicated genes
#'   rnaseq_mat <- rnaseq_mat[
#'     !duplicated(rownames(rnaseq_mat)), ]
#'
#'   # Remove duplicated samples (if needed)
#'   rnaseq_mat <- rnaseq_mat[
#'     , !duplicated(colnames(rnaseq_mat))]
#'
#'   # Transpose: samples x genes
#'   expr <- t(rnaseq_mat)
#'
#'   # Select top variable genes
#'   vars <- matrixStats::colVars(expr)
#'   n_top <- min(7000, ncol(expr))
#'   expr_var <- expr[, order(vars, decreasing = TRUE)[seq_len(n_top)]]
#'
#'   # Log-transform and scale
#'   expr_var <- log2(expr_var + 0.001)
#'   expr_var <- scale(expr_var)
#'
#'   ## store UMAP training settings
#'   ## for future reference/projection of new samples
#'
#'   ## 1. store gene set used for UMAP training
#'   umap_train_settings <- list()
#'   umap_train_settings[['gene_set']] <-
#'     colnames(expr_var)
#'
#'   ## 2. save scaling attributes used for UMAP training data
#'   umap_train_settings[['scale_center']] <-
#'     attr(expr_var, "scaled:center")
#'   umap_train_settings[['scale_scale']] <-
#'     attr(expr_var, "scaled:scale")
#'
#'
#'   umap_model <- uwot::umap(
#'     expr_var,
#'     n_neighbors = n_neighbors,
#'     min_dist = min_dist,
#'     spread = spread,
#'     metric = metric,
#'     scale = FALSE,
#'     n_components = 2,
#'     verbose = FALSE,
#'     ret_model = TRUE,
#'     n_threads = 2
#'   )
#'
#'   umap_df <- as.data.frame(umap_model$embedding)
#'   colnames(umap_df) <-
#'     c("UMAP1", "UMAP2")
#'   umap_df$tumor_sample_barcode <-
#'     rownames(umap_df)
#'   umap_df <- umap_df |>
#'     dplyr::inner_join(
#'       dplyr::select(
#'         tcga_clinical,
#'         c("tumor_sample_barcode",
#'           "sample_type",
#'           "primary_site",
#'           "primary_diagnosis_simplified",
#'           "tumor_stage_TNM",
#'           "tumor",
#'           "project_id")
#'       ), by = "tumor_sample_barcode"
#'     ) |>
#'     dplyr::mutate(tumor_type = primary_site) |>
#'
#'     ## ignore tumor types with few samples
#'     dplyr::filter(
#'       tumor_type != "Bone" &
#'         tumor_type != "Other/Unknown"
#'     ) |>
#'     dplyr::mutate(
#'       hover_text = glue::glue(
#'         "<b>Sample type:</b> {sample_type}<br>",
#'         "<b>Primary tumor site:</b> {primary_site}<br>",
#'         "<b>Tumor:</b> {project_id}<br>",
#'         "<b>Primary diagnosis:</b> {primary_diagnosis_simplified}<br>"
#'       )
#'     ) |>
#'     dplyr::distinct()
#'
#'
#'   base_cols <- ggsci::pal_jco("default")(10)
#'   cancer_types <- sort(unique(umap_df$tumor_type))
#'   extended_cols <- grDevices::colorRampPalette(
#'     base_cols)(length(cancer_types))
#'
#'   umap_plot <- list()
#'   umap_plot[['static']] <- ggplot2::ggplot(
#'     umap_df, ggplot2::aes(
#'       x = UMAP1,
#'       y = UMAP2,
#'       color = tumor_type,
#'       text = hover_text)) +
#'     ggplot2::theme_classic(
#'       base_family = "Helvetica") +
#'     ggplot2::geom_point(size = 0.8, alpha = 0.6) +
#'     ggplot2::theme(
#'       axis.title.x = ggplot2::element_blank(),
#'       axis.title.y = ggplot2::element_blank(),
#'       axis.ticks = ggplot2::element_blank(),
#'       legend.position = "bottom",
#'       legend.title = ggplot2::element_blank(),
#'       legend.text = ggplot2::element_text(size=10)
#'     ) +
#'     ggplot2::scale_color_manual(values = extended_cols)
#'
#'   umap_plot[['interactive']] <- plotly::ggplotly(
#'     umap_plot[['static']],
#'     tooltip = "text") |>
#'     plotly::layout(
#'       legend = list(
#'         title = list(text = ""),
#'         orientation = "h",
#'         x = 0.5,
#'         xanchor = "center",
#'         y = -0.2
#'       ),
#'       xaxis = list(
#'         title = "",
#'         showticklabels = FALSE,
#'         showgrid = FALSE,
#'         zeroline = FALSE,
#'         showline = FALSE
#'       ),
#'       yaxis = list(
#'         title = "",
#'         showticklabels = FALSE,
#'         showgrid = FALSE,
#'         zeroline = FALSE,
#'         showline = FALSE
#'       )
#'     )
#'
#'
#'   umap_result <- list('training_settings' = umap_train_settings,
#'                       'plot' = umap_plot,
#'                       'embedding' = umap_df)
#'
#'   saveRDS(umap_result, file = output_fname)
#'   uwot_model <-
#'     uwot::save_uwot(umap_model, file = output_model_fname)
#'
#' }
#'
#'
#' #' Get or make global UMAP embedding of TCGA RNAseq data
#' #' from GDC
#' #'
#' #' @param gdc_release TCGA data release version
#' #' @param output_dir Directory to save/load processed RNAseq data
#' #' @param data_raw_dir Directory containing raw GDC data
#' #' @param gdc_projects Character vector - GDC/TCGA project IDs
#' #' @param overwrite Whether to overwrite existing processed data
#' #' @param include_impress_test Whether to include IMPRESS test samples in UMAP
#' #' @param n_neighbors UMAP n_neighbors parameter
#' #' @param min_dist UMAP min_dist parameter
#' #' @param spread UMAP spread parameter
#' #' @param metric UMAP metric parameter
#' #' @param batch_correct Whether to apply ComBat batch correction when including IMPRESS samples
#' #' @export
#' #'
#' #'
#' global_umap_tcga2 <- function(
#'     gdc_release = "release45_20251204",
#'     output_dir = NULL,
#'     data_raw_dir = NULL,
#'     gdc_projects = NULL,
#'     overwrite = FALSE,
#'     include_impress_test = FALSE,
#'     n_neighbors = 15,
#'     min_dist = 0.3,
#'     spread = 0.7,
#'     metric = "cosine",
#'     batch_correct = TRUE){
#'
#'   assertthat::assert_that(
#'     !is.null(output_dir),
#'     !is.null(data_raw_dir),
#'     !is.null(gdc_projects)
#'   )
#'   assertthat::assert_that(
#'     dir.exists(output_dir),
#'     dir.exists(data_raw_dir)
#'   )
#'
#'   output_fname <- file.path(
#'     output_dir,
#'     gdc_release,
#'     "rnaseq",
#'     "tcga_rnaseq_umap.rds"
#'   )
#'
#'   output_model_fname <- file.path(
#'     output_dir,
#'     gdc_release,
#'     "rnaseq",
#'     "tcga_rnaseq_umap_model.uwot"
#'   )
#'
#'
#'   if(file.exists(output_fname) &
#'      file.exists(output_model_fname) &
#'      overwrite == FALSE){
#'
#'     umap_result <- readRDS(file = output_fname)
#'     umap_result$model <- uwot::load_uwot(file = output_model_fname)
#'     return(umap_result)
#'   }
#'
#'   tcga_clinical <- gdc_clinical(
#'     gdc_release = gdc_release,
#'     overwrite = F,
#'     tumor_samples_only = TRUE,
#'     output_dir = output_dir)
#'
#'   all_expr <- data.frame()
#'
#'   for(gdc_project in gdc_projects){
#'     cat("Processing ", gdc_project, "\n")
#'     project <-
#'       stringr::str_replace_all(
#'         gdc_project, 'TCGA-','')
#'     fname <-
#'       file.path(
#'         output_dir,
#'         gdc_release,
#'         "rnaseq",
#'         glue::glue(
#'           "tcga_rnaseq_TPM_{project}.rds")
#'       )
#'     expr_df <- readRDS(file = fname)
#'     if(NROW(all_expr) == 0){
#'       all_expr <- expr_df
#'     }else{
#'       for(e in c("ENSEMBL_GENE_ID",
#'                  "SYMBOL",
#'                  "GENENAME",
#'                  "BIOTYPE")){
#'         expr_df[[e]] <- NULL
#'       }
#'       all_expr <- all_expr |>
#'         dplyr::left_join(
#'           expr_df, by = "ENTREZGENE"
#'         )
#'     }
#'
#'   }
#'
#'   if(include_impress_test == TRUE){
#'     impress_sample_data <- readRDS(
#'       file = file.path(
#'         data_raw_dir,
#'         "impress",
#'         "impress_test_exprdata.rds"
#'       ))
#'
#'     exprdata_impress <- impress_sample_data$data |>
#'       dplyr::select(
#'         -c("GENENAME",
#'            "BIOTYPE",
#'            "SYMBOL",
#'            "ENSEMBL_GENE_ID")
#'       )
#'
#'     all_expr <- all_expr |>
#'       dplyr::left_join(
#'         exprdata_impress,
#'         by = "ENTREZGENE"
#'       )
#'
#'     impress_clinical <- impress_sample_data$metadata |>
#'       dplyr::select(
#'         -c("pcgr_tsite_code")
#'       ) |>
#'       dplyr::mutate(
#'         project_id = "IMPRESS"
#'       ) |>
#'       dplyr::rename(
#'         tumor_sample_barcode = sample_id,
#'         #primary_site = primary_site2
#'       ) |>
#'       dplyr::select(
#'         c("tumor_sample_barcode",
#'           "sample_type",
#'           "primary_site",
#'           "project_id")
#'       ) |>
#'       dplyr::distinct()
#'
#'
#'     tcga_clinical <- tcga_clinical |>
#'       dplyr::bind_rows(impress_clinical)
#'   }
#'
#'   gene_ids <- all_expr$ENSEMBL_GENE_ID
#'
#'   all_expr$GENENAME <- NULL
#'   all_expr$BIOTYPE <- NULL
#'   all_expr$SYMBOL <- NULL
#'   all_expr$ENTREZGENE <- NULL
#'   all_expr$ENSEMBL_GENE_ID <- NULL
#'
#'   rnaseq_mat <- as.matrix(all_expr)
#'   storage.mode(rnaseq_mat) <- "numeric"
#'   rownames(rnaseq_mat) <- gene_ids
#'   sample_ids <- colnames(rnaseq_mat)
#'
#'   set.seed(1234)
#'
#'
#'   # Remove duplicated genes
#'   rnaseq_mat <- rnaseq_mat[
#'     !duplicated(rownames(rnaseq_mat)), ]
#'
#'   # Remove duplicated samples (if needed)
#'   rnaseq_mat <- rnaseq_mat[
#'     , !duplicated(colnames(rnaseq_mat))]
#'
#'   # Transpose: samples x genes
#'   expr <- t(rnaseq_mat)
#'
#'   # Select top variable genes (based on TCGA only if including IMPRESS)
#'   if(include_impress_test == TRUE & batch_correct == TRUE){
#'     tcga_samples <- grep("^TCGA", rownames(expr), value = TRUE)
#'     vars <- matrixStats::colVars(expr[tcga_samples, ])
#'   } else {
#'     vars <- matrixStats::colVars(expr)
#'   }
#'
#'   n_top <- min(7000, ncol(expr))
#'   expr_var <- expr[, order(vars, decreasing = TRUE)[seq_len(n_top)]]
#'
#'   # Log-transform (before batch correction)
#'   expr_var <- log2(expr_var + 1)
#'
#'   # Apply limma batch correction if IMPRESS samples are included
#'   if(include_impress_test == TRUE & batch_correct == TRUE){
#'
#'     cat("Applying batch correction with limma::removeBatchEffect...\n")
#'
#'     # Create batch vector
#'     batch <- ifelse(grepl("^TCGA", rownames(expr_var)), "TCGA", "IMPRESS")
#'
#'     if(!requireNamespace("limma", quietly = TRUE)){
#'       stop("Package 'limma' is required for batch correction. Please install it with:\n",
#'            "  BiocManager::install('limma')")
#'     }
#'
#'     # Apply removeBatchEffect (expects genes x samples)
#'     expr_corrected <- limma::removeBatchEffect(
#'       x = t(expr_var),
#'       batch = batch
#'     )
#'
#'     expr_var <- t(expr_corrected)  # Back to samples x genes
#'
#'     # Remove genes with NA values after batch correction
#'     na_genes <- colSums(is.na(expr_var)) > 0
#'     if(any(na_genes)){
#'       n_na <- sum(na_genes)
#'       cat(sprintf("Removing %d genes with NA values after batch correction\n", n_na))
#'       expr_var <- expr_var[, !na_genes]
#'     }
#'
#'     cat("Batch correction completed.\n")
#'   }
#'
#'   # Scale the data
#'   expr_var <- scale(expr_var)
#'
#'   ## store UMAP training settings
#'   ## for future reference/projection of new samples
#'
#'   ## 1. store gene set used for UMAP training
#'   umap_train_settings <- list()
#'   umap_train_settings[['gene_set']] <-
#'     colnames(expr_var)
#'
#'   ## 2. save scaling attributes used for UMAP training data
#'   umap_train_settings[['scale_center']] <-
#'     attr(expr_var, "scaled:center")
#'   umap_train_settings[['scale_scale']] <-
#'     attr(expr_var, "scaled:scale")
#'
#'   ## 3. save batch correction flag
#'   umap_train_settings[['batch_corrected']] <-
#'     (include_impress_test == TRUE & batch_correct == TRUE)
#'
#'
#'   umap_model <- uwot::umap(
#'     expr_var,
#'     n_neighbors = n_neighbors,
#'     min_dist = min_dist,
#'     spread = spread,
#'     metric = metric,
#'     scale = FALSE,
#'     n_components = 2,
#'     verbose = FALSE,
#'     ret_model = TRUE,
#'     n_threads = 2
#'   )
#'
#'   umap_df <- as.data.frame(umap_model$embedding)
#'   colnames(umap_df) <-
#'     c("UMAP1", "UMAP2")
#'   umap_df$tumor_sample_barcode <-
#'     rownames(umap_df)
#'   umap_df <- umap_df |>
#'     dplyr::inner_join(
#'       dplyr::select(
#'         tcga_clinical,
#'         c("tumor_sample_barcode",
#'           "sample_type",
#'           "primary_site",
#'           "primary_diagnosis_simplified",
#'           "tumor_stage_TNM",
#'           "tumor",
#'           "project_id")
#'       ), by = "tumor_sample_barcode"
#'     ) |>
#'     dplyr::mutate(tumor_type = primary_site) |>
#'
#'     ## ignore tumor types with few samples
#'     dplyr::filter(
#'       tumor_type != "Bone" &
#'         tumor_type != "Other/Unknown"
#'     ) |>
#'     dplyr::mutate(
#'       hover_text = glue::glue(
#'         "<b>Sample type:</b> {sample_type}<br>",
#'         "<b>Primary tumor site:</b> {primary_site}<br>",
#'         "<b>Tumor:</b> {project_id}<br>",
#'         "<b>Primary diagnosis:</b> {primary_diagnosis_simplified}<br>"
#'       )
#'     ) |>
#'     dplyr::distinct()
#'
#'
#'   base_cols <- ggsci::pal_jco("default")(10)
#'   cancer_types <- sort(unique(umap_df$tumor_type))
#'   extended_cols <- grDevices::colorRampPalette(
#'     base_cols)(length(cancer_types))
#'
#'   umap_plot <- list()
#'   umap_plot[['static']] <- ggplot2::ggplot(
#'     umap_df, ggplot2::aes(
#'       x = UMAP1,
#'       y = UMAP2,
#'       color = tumor_type,
#'       text = hover_text)) +
#'     ggplot2::theme_classic(
#'       base_family = "Helvetica") +
#'     ggplot2::geom_point(size = 0.8, alpha = 0.6) +
#'     ggplot2::theme(
#'       axis.title.x = ggplot2::element_blank(),
#'       axis.title.y = ggplot2::element_blank(),
#'       axis.ticks = ggplot2::element_blank(),
#'       legend.position = "bottom",
#'       legend.title = ggplot2::element_blank(),
#'       legend.text = ggplot2::element_text(size=10)
#'     ) +
#'     ggplot2::scale_color_manual(values = extended_cols)
#'
#'   umap_plot[['interactive']] <- plotly::ggplotly(
#'     umap_plot[['static']],
#'     tooltip = "text") |>
#'     plotly::layout(
#'       legend = list(
#'         title = list(text = ""),
#'         orientation = "h",
#'         x = 0.5,
#'         xanchor = "center",
#'         y = -0.2
#'       ),
#'       xaxis = list(
#'         title = "",
#'         showticklabels = FALSE,
#'         showgrid = FALSE,
#'         zeroline = FALSE,
#'         showline = FALSE
#'       ),
#'       yaxis = list(
#'         title = "",
#'         showticklabels = FALSE,
#'         showgrid = FALSE,
#'         zeroline = FALSE,
#'         showline = FALSE
#'       )
#'     )
#'
#'
#'   umap_result <- list('training_settings' = umap_train_settings,
#'                       'plot' = umap_plot,
#'                       'embedding' = umap_df)
#'
#'   saveRDS(umap_result, file = output_fname)
#'   uwot_model <-
#'     uwot::save_uwot(umap_model, file = output_model_fname)
#'
#' }


#' Get or train global UMAP embedding of TCGA RNAseq data
#'
#' @param gdc_release TCGA data release version
#' @param output_dir Directory to save/load processed RNAseq data
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_projects Character vector - GDC/TCGA project IDs
#' @param overwrite Whether to overwrite existing processed data
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param spread UMAP spread parameter
#' @param metric UMAP metric parameter
#' @param n_reference Total number of reference samples (default: 2000)
#' @param min_per_primary Minimum samples per primary site (default: 20)
#' @param min_per_subtype Minimum samples per subtype (default: 3)
#' @param max_per_subtype Maximum samples per subtype (default: 50)
#' @export
train_tcga_umap <- function(
    gdc_release = "release45_20251204",
    output_dir = NULL,
    data_raw_dir = NULL,
    gdc_projects = NULL,
    overwrite = FALSE,
    n_neighbors = 15,
    min_dist = 0.3,
    spread = 0.7,
    metric = "cosine",
    n_top_variable_genes = 6000,
    n_reference = 5000,
    min_per_primary = 60,
    min_per_subtype = 20,
    max_per_subtype = 90){

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
    "tcga_rnaseq_umap.rds"
  )

  output_model_fname <- file.path(
    output_dir,
    gdc_release,
    "rnaseq",
    "tcga_rnaseq_umap_model.uwot"
  )

  output_reference_fname <- file.path(
    output_dir,
    gdc_release,
    "rnaseq",
    "tcga_reference_samples.rds"
  )

  output_centroids_fname <- file.path(
    output_dir,
    gdc_release,
    "rnaseq",
    "tcga_centroids.rds"
  )


  if(file.exists(output_fname) &
     file.exists(output_model_fname) &
     file.exists(output_reference_fname) &
     file.exists(output_centroids_fname) &
     overwrite == FALSE){

    umap_result <- readRDS(file = output_fname)
    umap_result$model <- uwot::load_uwot(file = output_model_fname)
    umap_result$reference_samples <- readRDS(file = output_reference_fname)
    umap_result$centroids <- readRDS(file = output_centroids_fname)
    return(umap_result)
  }

  tcga_clinical <- gdc_clinical(
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

  # Select top variable genes (TCGA only)
  vars <- matrixStats::colVars(expr)
  n_top <- min(n_top_variable_genes, ncol(expr))
  expr_var <- expr[, order(vars, decreasing = TRUE)[seq_len(n_top)]]

  # Log-transform
  expr_var <- log2(expr_var + 1)

  expr_var_unscaled <- expr_var

  # Scale the data
  expr_var <- scale(expr_var)

  ## store UMAP training settings
  umap_train_settings <- list()
  umap_train_settings[['gene_set']] <-
    colnames(expr_var)
  umap_train_settings[['scale_center']] <-
    attr(expr_var, "scaled:center")
  umap_train_settings[['scale_scale']] <-
    attr(expr_var, "scaled:scale")
  umap_train_settings[['batch_corrected']] <- FALSE

  # Train UMAP model
  umap_model <- uwot::umap(
    expr_var,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    metric = metric,
    scale = FALSE,
    n_components = 2,
    verbose = FALSE,
    ret_model = TRUE,
    n_threads = 2
  )

  umap_df <- as.data.frame(umap_model$embedding)
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
          "project_id",
          "subtype_selected")
      ), by = "tumor_sample_barcode"
    ) |>
    dplyr::mutate(tumor_type = primary_site) |>
    dplyr::filter(
      tumor_type != "Bone" &
        tumor_type != "Other/Unknown"
    ) |>
    dplyr::mutate(
      hover_text = glue::glue(
        "<b>Sample type:</b> {sample_type}<br>",
        "<b>Primary tumor site:</b> {primary_site}<br>",
        "<b>Tumor:</b> {project_id}<br>",
        "<b>Primary diagnosis:</b> {primary_diagnosis_simplified}<br>",
        "<b>Subtype:</b> {ifelse(is.na(subtype_selected), 'N/A', subtype_selected)}<br>"
      )
    ) |>
    dplyr::distinct()

  # SELECT REFERENCE SAMPLES
  cat("\nSelecting reference samples for batch correction...\n")

  reference_result <- select_subtype_aware_reference_samples(
    tcga_clinical = umap_df,  # Use filtered clinical data
    n_total = n_reference,
    min_per_primary = min_per_primary,
    min_per_subtype = min_per_subtype,
    max_per_subtype = max_per_subtype
  )

  reference_sample_ids <- reference_result$sample_ids

  # Extract reference expression data (raw TPM before log/scale)
  reference_expr_raw <- expr[reference_sample_ids, ]

  reference_samples <- list(
    sample_ids = reference_sample_ids,
    expression_raw = reference_expr_raw,
    clinical = umap_df[umap_df$tumor_sample_barcode %in% reference_sample_ids, ],
    summary = reference_result$summary
  )

  cat(sprintf("Selected %d reference samples\n", length(reference_sample_ids)))

  # CREATE CENTROIDS
  cat("\nCreating expression centroids for primary sites and subtypes...\n")

  centroids <- create_expression_centroids(
    expression_matrix = expr_var_unscaled,  # Use scaled/unscaled data
    clinical_data = umap_df,
    training_genes = umap_train_settings[['gene_set']]
  )

  cat(sprintf("Created %d centroids\n", nrow(centroids$metadata)))

  # CREATE PLOTS
  base_cols <- ggsci::pal_jco("default")(10)
  cancer_types <- sort(unique(umap_df$tumor_type))
  extended_cols <- grDevices::colorRampPalette(
    base_cols)(length(cancer_types))

  umap_plot <- list()
  umap_plot[['static']] <- ggplot2::ggplot(
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

  umap_plot[['interactive']] <- plotly::ggplotly(
    umap_plot[['static']],
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

  umap_result <- list(
    'training_settings' = umap_train_settings,
    'plot' = umap_plot,
    'embedding' = umap_df
    #'reference_samples' = reference_samples,
    #'centroids' = centroids
  )

  saveRDS(umap_result, file = output_fname)
  saveRDS(reference_samples, file = output_reference_fname)
  saveRDS(centroids, file = output_centroids_fname)
  uwot::save_uwot(umap_model, file = output_model_fname)

  cat("\nTraining complete!\n")
  cat(sprintf("  UMAP embedding: %d samples\n", nrow(umap_df)))
  cat(sprintf("  Reference samples: %d samples\n", length(reference_sample_ids)))
  cat(sprintf("  Centroids: %d groups\n", nrow(centroids$metadata)))

  invisible(umap_result)
}

#' Select subtype-aware reference samples
#'
select_subtype_aware_reference_samples <- function(
    tcga_clinical,
    n_total = 2000,
    min_per_primary = 50,
    min_per_subtype = 20,
    max_per_subtype = 90) {

  set.seed(123)

  primary_subtype_counts <- tcga_clinical |>
    dplyr::group_by(primary_site, subtype_selected) |>
    dplyr::summarise(n_available = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(primary_site, dplyr::desc(n_available))

  reference_samples <- c()

  for(site in unique(primary_subtype_counts$primary_site)) {

    site_data <- primary_subtype_counts |>
      dplyr::filter(primary_site == site)

    site_total_samples <- sum(site_data$n_available)
    has_subtypes <- any(!is.na(site_data$subtype_selected))

    if(!has_subtypes || all(is.na(site_data$subtype_selected))) {

      # No subtypes - proportional sampling
      n_site <- max(
        min_per_primary,
        round(n_total * site_total_samples / nrow(tcga_clinical))
      )

      site_samples <- tcga_clinical |>
        dplyr::filter(primary_site == site) |>
        dplyr::pull(tumor_sample_barcode)

      selected <- sample(site_samples, min(n_site, length(site_samples)))
      reference_samples <- c(reference_samples, selected)

    } else {

      # Has subtypes - ensure representation
      for(i in 1:nrow(site_data)) {

        subtype <- site_data$subtype_selected[i]
        n_available <- site_data$n_available[i]

        if(is.na(subtype)) {
          n_select <- min(
            max(min_per_primary, round(n_available * 0.3)),
            n_available
          )
        } else {
          if(n_available < min_per_subtype) {
            n_select <- n_available
          } else if(n_available > max_per_subtype) {
            n_select <- max_per_subtype
          } else {
            n_select <- max(
              min_per_subtype,
              round(n_available * 0.4)
            )
          }
        }

        subtype_samples <- tcga_clinical |>
          dplyr::filter(
            primary_site == site,
            (is.na(subtype_selected) & is.na(subtype)) |
              (!is.na(subtype_selected) & subtype_selected == subtype)
          ) |>
          dplyr::pull(tumor_sample_barcode)

        selected <- sample(subtype_samples, n_select)
        reference_samples <- c(reference_samples, selected)
      }
    }
  }

  ref_summary <- tcga_clinical |>
    dplyr::filter(tumor_sample_barcode %in% reference_samples) |>
    dplyr::group_by(primary_site, subtype_selected) |>
    dplyr::summarise(n_ref = dplyr::n(), .groups = "drop") |>
    dplyr::left_join(primary_subtype_counts, by = c("primary_site", "subtype_selected")) |>
    dplyr::mutate(pct_sampled = round(100 * n_ref / n_available, 1))

  cat(sprintf("Selected %d reference samples across %d primary sites\n",
              length(reference_samples),
              length(unique(ref_summary$primary_site))))

  return(list(
    sample_ids = reference_samples,
    summary = ref_summary
  ))
}


#' Create expression centroids for primary sites and subtypes
#'
#' @param expression_matrix Scaled expression matrix (samples x genes)
#' @param clinical_data Clinical data with primary_site and subtype_selected
#' @param training_genes Gene IDs used in training
#' @param min_samples_for_centroids Minimum samples required to create centroid
#'
#' @return List with centroid_matrix and metadata
create_expression_centroids <- function(
    expression_matrix = NULL,
    clinical_data = NULL,
    training_genes = NULL,
    min_samples_for_centroids = 25) {

  # Ensure clinical data matches expression matrix
  clinical_data <- clinical_data |>
    dplyr::filter(tumor_sample_barcode %in% rownames(expression_matrix))

  # Create grouping variable
  clinical_data <- clinical_data |>
    dplyr::mutate(
      centroid_group = dplyr::case_when(
        !is.na(subtype_selected) ~
          paste0(primary_site, "::", subtype_selected),
        TRUE ~ primary_site
      )
    )

  # Calculate centroids for each group
  centroid_list <- list()
  centroid_metadata <- data.frame()

  for(group in unique(clinical_data$centroid_group)) {

    group_samples <- clinical_data |>
      dplyr::filter(centroid_group == group) |>
      dplyr::pull(tumor_sample_barcode)

    # Skip if too few samples
    if(length(group_samples) < min_samples_for_centroids) {
      cat(sprintf("  Skipping %s (only %d samples)\n", group, length(group_samples)))
      next
    }

    # Calculate median expression across samples (more robust than mean)
    group_expr <- expression_matrix[group_samples, , drop = FALSE]
    centroid <- apply(group_expr, 2, median, na.rm = TRUE)

    centroid_list[[group]] <- centroid

    # Extract primary site and subtype
    parts <- strsplit(group, "::")[[1]]

    centroid_metadata <- rbind(
      centroid_metadata,
      data.frame(
        centroid_id = group,
        primary_site = parts[1],
        subtype = ifelse(length(parts) > 1, parts[2], NA),
        n_samples = length(group_samples),
        stringsAsFactors = FALSE
      )
    )
  }

  # Convert to matrix
  centroid_matrix <- do.call(rbind, centroid_list)
  colnames(centroid_matrix) <- training_genes

  cat(sprintf("Created %d centroids:\n", nrow(centroid_matrix)))
  cat(sprintf("  - %d primary site-only centroids\n",
              sum(is.na(centroid_metadata$subtype))))
  cat(sprintf("  - %d subtype-specific centroids\n",
              sum(!is.na(centroid_metadata$subtype))))

  return(list(
    centroids = centroid_matrix,
    metadata = centroid_metadata
  ))
}


#' Calculate similarities between sample and centroids
#'
#' @param sample_vector Scaled expression vector for single sample
#' @param centroid_matrix Matrix of centroids (rows = centroids, cols = genes)
#' @param centroid_metadata Metadata for centroids
#' @param method Similarity method: "spearman", "pearson", or "cosine"
#' @param min_samples_for_comparison Minimum samples required for comparison
#'
#' @return Data frame with similarities ranked by score
calculate_centroid_similarities <- function(
    sample_vector = NULL,
    centroid_matrix = NULL,
    centroid_metadata = NULL,
    method = "spearman",
    min_samples_for_comparison = 25) {

  # Filter out small centroids
  valid_centroids <- centroid_metadata$n_samples >= min_samples_for_comparison

  if(sum(!valid_centroids) > 0) {
    cat(sprintf("  Excluding %d centroids with < %d samples\n",
                sum(!valid_centroids), min_samples_for_comparison))
  }

  centroid_matrix_filtered <- centroid_matrix[valid_centroids, , drop = FALSE]
  centroid_metadata_filtered <- centroid_metadata[valid_centroids, ]

  similarities <- numeric(nrow(centroid_matrix_filtered))

  for(i in 1:nrow(centroid_matrix_filtered)) {

    if(method == "spearman") {
      similarities[i] <- cor(
        sample_vector,
        centroid_matrix_filtered[i, ],
        method = "spearman",
        use = "complete.obs"
      )
    } else if(method == "pearson") {
      similarities[i] <- cor(
        sample_vector,
        centroid_matrix_filtered[i, ],
        method = "pearson",
        use = "complete.obs"
      )
    } else if(method == "cosine") {
      # Cosine similarity
      similarities[i] <- sum(sample_vector * centroid_matrix_filtered[i, ]) /
        (sqrt(sum(sample_vector^2)) * sqrt(sum(centroid_matrix_filtered[i, ]^2)))
    }
  }

  # CRITICAL FIX: Use filtered metadata
  result <- centroid_metadata_filtered  #
  result$similarity <- similarities
  result$method <- "Centroid"

  # Rank by similarity
  result <- result |>
    dplyr::arrange(dplyr::desc(similarity))

  return(result)
}


#' Classify single sample using centroid-based similarity
#'
#' @param new_sample Named numeric vector of TPM values
#' @param sample_id Sample identifier
#' @param umap_train_settings Training settings (for gene set and scaling parameters)
#' @param reference_samples Reference sample data for batch correction
#' @param centroids Centroid data for similarity-based classification
#' @param n_reference_for_batch Number of reference samples for batch correction (default: 1000)
#' @param use_batch_correction Whether to apply batch correction (default: TRUE)
#' @param sample_metadata Optional metadata for the sample
#'
#' @return List with centroid classification results
#' @export
classify_sample_by_centroid <- function(
    new_sample = NULL,
    sample_id = "QUERY_SAMPLE",
    umap_train_settings = NULL,
    reference_samples = NULL,
    centroids = NULL,
    n_reference_for_batch = 1000,
    use_batch_correction = TRUE,
    sample_metadata = NULL) {

  # Convert to matrix if vector
  if(is.vector(new_sample)) {
    genenames <- names(new_sample)
    new_sample_mat <- matrix(new_sample, nrow = 1)
    colnames(new_sample_mat) <- genenames
    rownames(new_sample_mat) <- sample_id
  } else {
    new_sample_mat <- new_sample
  }

  training_genes <- umap_train_settings[['gene_set']]

  # ===== HANDLE MISSING GENES =====

  available_genes <- intersect(training_genes, colnames(new_sample_mat))
  missing_genes <- setdiff(training_genes, colnames(new_sample_mat))

  if(length(missing_genes) > 0) {
    cat(sprintf("Warning: %d/%d training genes missing in sample %s\n",
                length(missing_genes), length(training_genes), sample_id))
    cat("Missing genes will be imputed as zero TPM\n")

    # Create complete matrix with all training genes
    new_sample_complete <- matrix(0, nrow = 1, ncol = length(training_genes))
    colnames(new_sample_complete) <- training_genes
    rownames(new_sample_complete) <- sample_id

    # Fill in available genes
    new_sample_complete[, available_genes] <- new_sample_mat[, available_genes]
    new_sample_genes <- new_sample_complete

  } else {
    # All genes present - reorder to match training genes
    new_sample_genes <- new_sample_mat[, training_genes, drop = FALSE]
  }

  cat(sprintf("Sample prepared: %d genes (%d imputed as zero)\n",
              ncol(new_sample_genes), length(missing_genes)))

  # ===== PROCESSING PIPELINE =====

  if(use_batch_correction) {

    cat("\nApplying batch correction...\n")

    # Get reference expression data
    reference_expr_raw <- reference_samples$expression_raw
    reference_clinical <- reference_samples$clinical

    cat(sprintf("  Available reference samples: %d\n", nrow(reference_expr_raw)))

    # Subsample reference samples for batch correction
    if(nrow(reference_expr_raw) > n_reference_for_batch) {

      cat(sprintf("  Subsampling to %d samples...\n", n_reference_for_batch))

      # Stratified sampling by primary site
      sites_needed <- unique(reference_clinical$primary_site)
      n_per_site <- ceiling(n_reference_for_batch / length(sites_needed))

      site_counts <- reference_clinical |>
        dplyr::group_by(primary_site) |>
        dplyr::summarise(n_available = dplyr::n(), .groups = "drop")

      reference_subset_ids <- c()

      for(site in sites_needed) {
        site_n_available <- site_counts$n_available[site_counts$primary_site == site]
        n_to_sample <- min(n_per_site, site_n_available)

        site_samples <- reference_clinical |>
          dplyr::filter(primary_site == site) |>
          dplyr::pull(tumor_sample_barcode)

        selected <- sample(site_samples, n_to_sample)
        reference_subset_ids <- c(reference_subset_ids, selected)
      }

      if(length(reference_subset_ids) > n_reference_for_batch) {
        reference_subset_ids <- sample(reference_subset_ids, n_reference_for_batch)
      }

      reference_subset <- reference_expr_raw[reference_subset_ids, training_genes]

      cat(sprintf("  Sampled %d samples from %d sites\n",
                  length(reference_subset_ids), length(sites_needed)))

    } else {
      cat(sprintf("  Using all %d reference samples\n", nrow(reference_expr_raw)))
      reference_subset <- reference_expr_raw[, training_genes]
      reference_subset_ids <- rownames(reference_expr_raw)
    }

    # Combine and log-transform
    combined <- rbind(reference_subset, new_sample_genes)
    combined_log <- log2(combined + 1)

    # Batch correction
    batch <- c(rep("TCGA", nrow(reference_subset)), "NEW")

    combined_corrected <- limma::removeBatchEffect(
      x = t(combined_log),
      batch = batch
    )

    cat(sprintf("  Batch correction complete\n"))

    # Transpose back: samples x genes
    combined_corrected <- t(combined_corrected)

    # Extract corrected samples
    reference_corrected <- combined_corrected[1:(nrow(combined_corrected)-1), , drop = FALSE]
    new_sample_corrected <- combined_corrected[nrow(combined_corrected), , drop = FALSE]

    # Remove NA genes if any
    na_genes <- is.na(new_sample_corrected[1, ])
    if(any(na_genes)) {
      cat(sprintf("  Removing %d NA genes\n", sum(na_genes)))
      new_sample_corrected <- new_sample_corrected[, !na_genes, drop = FALSE]
      reference_corrected <- reference_corrected[, !na_genes, drop = FALSE]
    }

    # Recompute scaling parameters on batch-corrected reference
    corrected_center <- colMeans(reference_corrected, na.rm = TRUE)
    corrected_scale <- apply(reference_corrected, 2, sd, na.rm = TRUE)
    corrected_scale[corrected_scale == 0] <- 1

    # Scale using batch-corrected parameters
    new_sample_scaled <- scale(
      new_sample_corrected,
      center = corrected_center,
      scale = corrected_scale
    )

    scale_genes <- colnames(new_sample_scaled)

  } else {

    cat("\nProcessing WITHOUT batch correction...\n")

    # Log-transform
    new_sample_log <- log2(new_sample_genes + 1)

    # Scale using TCGA training parameters
    new_sample_scaled <- scale(
      new_sample_log,
      center = umap_train_settings[['scale_center']],
      scale = umap_train_settings[['scale_scale']]
    )

    scale_genes <- training_genes
  }

  cat(sprintf("Sample scaled: %d genes\n", ncol(new_sample_scaled)))

  # ===== CENTROID CLASSIFICATION =====

  cat("\nCalculating centroid similarities...\n")

  centroid_similarities <- calculate_centroid_similarities(
    sample_vector = new_sample_scaled[1, ],
    centroid_matrix = centroids$centroids[, scale_genes],
    centroid_metadata = centroids$metadata,
    method = "spearman"
  )

  # ===== FORMAT RESULTS =====

  # Create sample info data frame
  sample_info <- data.frame(
    sample_id = sample_id,
    n_genes = ncol(new_sample_scaled),
    n_missing_genes = length(missing_genes),
    pct_missing = round(100 * length(missing_genes) / length(training_genes), 2),
    batch_corrected = use_batch_correction,
    stringsAsFactors = FALSE
  )

  # Add metadata if provided
  if(!is.null(sample_metadata)) {
    for(col in names(sample_metadata)) {
      sample_info[[col]] <- sample_metadata[[col]]
    }
  }

  # ===== DISPLAY RESULTS =====

  cat("\n=== CLASSIFICATION RESULTS ===\n")
  cat(sprintf("Sample: %s\n", sample_id))
  cat(sprintf("Batch correction: %s\n", ifelse(use_batch_correction, "ENABLED", "DISABLED")))
  if(!is.null(sample_metadata) && "known_site" %in% names(sample_metadata)) {
    cat(sprintf("Known site: %s\n", sample_metadata$known_site))
  }
  cat("\nTop 5 predictions:\n")
  print(head(centroid_similarities, 5))

  # ===== RETURN RESULTS =====

  result <- list(
    sample_info = sample_info,
    classification = centroid_similarities,
    top_prediction = list(
      primary_site = centroid_similarities$primary_site[1],
      subtype = centroid_similarities$subtype[1],
      similarity = centroid_similarities$similarity[1],
      n_samples = centroid_similarities$n_samples[1]
    )
  )

  return(result)
}
