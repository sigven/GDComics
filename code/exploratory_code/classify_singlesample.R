source('code/utils.R')
source('code/rnaseq.R')

#' Plot TCGA UMAP with projected new sample(s) highlighted
#'
#' @param umap_result Result object from global_umap_tcga()
#' @param new_projections Single projection result or list of projection results from project_sample_to_tcga_umap()
#' @param highlight_color Color for highlighted new samples (default = "yellow")
#' @param highlight_size Size of highlighted points (default = 3)
#' @param highlight_shape Shape for highlighted points (default = 23, filled diamond)
#' @param show_interactive Whether to return plotly object (TRUE) or ggplot (FALSE)
#'
#' @return ggplot2 or plotly object with highlighted new sample(s)
#' @export
#'
# plot_umap_with_projection <- function(
#     umap_result,
#     new_projections = NULL,
#     highlight_color = "yellow",
#     highlight_size = 3,
#     highlight_shape = 23,
#     show_interactive = TRUE) {
#
#   # Get original embedding data
#   umap_df <- umap_result$embedding
#   umap_df$is_new_sample <- FALSE
#
#   # Process new projections
#   if (!is.null(new_projections)) {
#     # Handle single projection or list of projections
#     if ("projection_data" %in% names(new_projections)) {
#       # Single projection
#       new_data <- new_projections$projection_data
#     } else {
#       # List of projections
#       new_data <- dplyr::bind_rows(
#         lapply(new_projections, function(x) x$projection_data)
#       )
#     }
#
#     # Combine with original data
#     umap_df <- dplyr::bind_rows(umap_df, new_data)
#   }
#
#   # Get color palette from original plot
#   base_cols <- ggsci::pal_jco("default")(10)
#   cancer_types <- sort(unique(umap_df$tumor_type[!umap_df$is_new_sample]))
#   extended_cols <- grDevices::colorRampPalette(base_cols)(length(cancer_types))
#   names(extended_cols) <- cancer_types
#
#   # Create base plot matching your original style
#   p <- ggplot2::ggplot(
#     umap_df,
#     ggplot2::aes(
#       x = UMAP1,
#       y = UMAP2,
#       color = tumor_type,
#       text = hover_text
#     )) +
#     ggplot2::theme_classic(base_family = "Helvetica") +
#     ggplot2::geom_point(
#       data = umap_df[!umap_df$is_new_sample, ],
#       size = 0.8,
#       alpha = 0.6
#     ) +
#     ggplot2::theme(
#       axis.title.x = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_blank(),
#       axis.ticks = ggplot2::element_blank(),
#       legend.position = "bottom",
#       legend.title = ggplot2::element_blank(),
#       legend.text = ggplot2::element_text(size = 10)
#     ) +
#     ggplot2::scale_color_manual(values = extended_cols)
#
#   # Add highlighted new samples if present
#   if (!is.null(new_projections)) {
#     p <- p +
#       ggplot2::geom_point(
#         data = umap_df[umap_df$is_new_sample, ],
#         ggplot2::aes(x = UMAP1, y = UMAP2),
#         color = "black",
#         fill = highlight_color,
#         size = highlight_size,
#         shape = highlight_shape,
#         stroke = 1.5,
#         inherit.aes = FALSE
#       ) +
#       ggrepel::geom_text_repel(
#         data = umap_df[umap_df$is_new_sample, ],
#         ggplot2::aes(x = UMAP1, y = UMAP2, label = tumor_sample_barcode),
#         size = 3.5,
#         fontface = "bold",
#         box.padding = 0.5,
#         point.padding = 0.5,
#         inherit.aes = FALSE
#       )
#   }
#
#   # Return interactive or static plot
#   if (show_interactive) {
#     p_interactive <- plotly::ggplotly(p, tooltip = "text") |>
#       plotly::layout(
#         legend = list(
#           title = list(text = ""),
#           orientation = "h",
#           x = 0.5,
#           xanchor = "center",
#           y = -0.2
#         ),
#         xaxis = list(
#           title = "",
#           showticklabels = FALSE,
#           showgrid = FALSE,
#           zeroline = FALSE,
#           showline = FALSE
#         ),
#         yaxis = list(
#           title = "",
#           showticklabels = FALSE,
#           showgrid = FALSE,
#           zeroline = FALSE,
#           showline = FALSE
#         )
#       )
#     return(p_interactive)
#   } else {
#     return(p)
#   }
# }

#
gdc_release <- "release45_20251204"

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

# Load your existing UMAP result
umap_result <- train_tcga_umap(
  data_raw_dir = file.path(here::here(), "data-raw"),
  output_dir = file.path(here::here(), "output"),
  gdc_release = gdc_release,
  gdc_projects = gdc_projects,
  overwrite = FALSE
)

data_raw_dir <- file.path(
  here::here(), "data-raw"
)

# Get Gencode Xref -----
gencode_xref <- get_gencode_xref(
  data_raw_dir = data_raw_dir,
  update_ensembl_gene = FALSE
)

# Load new sample TPM data (assuming gene names as rownames)
impress_path <-
  "/Users/sigven/project_data/data/data__impress/impress/test_samples/"
expression_files <-
  list.files(path = impress_path, pattern = "\\.gene_expression.grch38.tsv")



i <- 1
impress_expr_data <- data.frame()
while(i <= 9){
  expression_fname <-
    file.path(impress_path, "SECOND_SET",
              paste0("WGS-T00", i,".gene_expression.tsv"))
  new_sample_tpm <-
    readr::read_tsv(expression_fname, show_col_types = F) |>
    dplyr::select(c("TargetID","TPM")) |>
    dplyr::rename(ensembl_transcript_id = TargetID) |>
    dplyr::mutate(ensembl_transcript_id = stringr::str_remove(
      ensembl_transcript_id, "\\.*$")
    ) |>
    dplyr::mutate(
      SAMPLE_ID = paste0("IMPRESS-WGS-T00", i)
    ) |>
    dplyr::inner_join(
      gencode_xref[['transcript']],
      by = c("ensembl_transcript_id" = "ENSEMBL_TRANSCRIPT_ID")
    ) |>
    dplyr::group_by(
      ENSEMBL_GENE_ID,
      SYMBOL,
      ENTREZGENE,
      GENENAME,
      BIOTYPE,
      SAMPLE_ID) |>
    dplyr::summarise(TPM = sum(TPM), .groups = "drop") |>
    dplyr::ungroup()

  i <- i + 1
  impress_expr_data <- dplyr::bind_rows(
    impress_expr_data, new_sample_tpm)

}

impress_expr_data_wide <- as.data.frame(
  impress_expr_data |>
  dplyr::arrange(ENTREZGENE) |>
  tidyr::pivot_wider(
    names_from = "SAMPLE_ID",
    values_from = "TPM"))

#new_sample_tpm <- readRDS("new_sample_tpm.rds")  # or read from file

# Optional: provide metadata for the new sample
sample_metadata <- readr::read_tsv(
  file = file.path(impress_path, "RAW_DATA", "clinical.txt"),
  show_col_types = F) |>
  dplyr::mutate(primary_site2 = paste0(
    "IMPRESS-", primary_site
  ))

impress_samples <- list()
impress_samples[['data']] <-
  impress_expr_data_wide
impress_samples[['metadata']] <-
  sample_metadata

saveRDS(
  impress_samples,
  file = "data-raw/impress_test_exprdata.rds"
)


source('code/rnaseq.R')

gdc_release <- "release45_20251204"
data_raw_dir <- file.path(
  here::here(), "data-raw"
)
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

tcga_umap_model_fname <-
  file.path(here::here(), "output", gdc_release,
            "rnaseq", "tcga_rnaseq_umap_model.uwot")

tcga_umap_data_fname <-
  file.path(here::here(), "output", gdc_release,
            "rnaseq", "tcga_rnaseq_umap.rds")

tcga_centroids_fname <-
  file.path(here::here(), "output", gdc_release,
            "rnaseq", "tcga_centroids.rds")

tcga_reference_samples_fname <-
  file.path(here::here(), "output", gdc_release,
            "rnaseq", "tcga_reference_samples.rds")



# Load for projection
umap_model <- uwot::load_uwot(tcga_umap_model_fname)
umap_data <- readRDS(file = tcga_umap_data_fname)
umap_data[['centroids']] <-
  readRDS(file = tcga_centroids_fname)
umap_data[['reference_samples']] <-
  readRDS(file = tcga_reference_samples_fname)

impress_sample_data <- readRDS(
  file = file.path(
    data_raw_dir,
    "impress",
    "impress_test_exprdata.rds"
  ))

# Extract individual samples as named vectors

samples <- list()

samples[['Uterus']] <- impress_sample_data$data$`WGS-T001`
names(samples[['Uterus']]) <- impress_sample_data$data$ENSEMBL_GENE_ID
samples[['Ovarian']] <- impress_sample_data$data$`WGS-T002`
names(samples[['Ovarian']]) <- impress_sample_data$data$ENSEMBL_GENE_ID
samples[['Bladder']] <- impress_sample_data$data$`WGS-T004`
names(samples[['Bladder']]) <- impress_sample_data$data$ENSEMBL_GENE_ID
samples[['Thyroid']] <- impress_sample_data$data$`WGS-T005`
names(samples[['Thyroid']]) <- impress_sample_data$data$ENSEMBL_GENE_ID
samples[['CRC']] <- impress_sample_data$data$`WGS-T007`
names(samples[['CRC']]) <- impress_sample_data$data$ENSEMBL_GENE_ID


# Classify a sample
for(sample_type in c('Uterus',
                     'Ovarian',
                     'Bladder',
                     'Thyroid',
                     'CRC')){

  result <- classify_sample_by_centroid(
    new_sample = samples[[sample_type]],
    sample_id = "IMPR-WGS-XXX",
    umap_train_settings = umap_data$training_settings,
    reference_samples = umap_data$reference_samples,
    centroids = umap_data$centroids,
    n_reference_for_batch = 1000,
    use_batch_correction = TRUE,
    sample_metadata = list(known_site = sample_type, stage = "III")
  )

  # Top prediction
  cat("\nTop prediction:\n")
  cat(sprintf("  Primary site: %s\n", result$top_prediction$primary_site))
  cat(sprintf("  Subtype: %s\n", result$top_prediction$subtype))
  cat(sprintf("  Similarity: %.3f\n", result$top_prediction$similarity))
  cat(sprintf("  Based on %d samples\n", result$top_prediction$n_samples))


}

# View results
#print(result$sample_info)
#print(result$classification)

