# ## Core dependencies
# library(TCGAbiolinks)
# library(valr)
# library(SummarizedExperiment)
# library(dplyr)
# library(stringr)
# library(tidyr)
# library(GenomeInfoDb)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(reshape2)
# library(ggplot2)
# library(ggpubr)
# library(caret)
# library(RTCGAToolbox)
#
# #remotes::install_github('sigven/geneOncoX')
# #remotes::install_github('sigven/cancerHotspots')
# #remotes::install_github('sigven/vcfHelpR')
#
# ## GDC omics data retrieval scripts
# source('code/utils.R')
# source('code/clinical.R')
# source('code/mutation.R')
# source('code/msi.R')
# source('code/cna.R')
# source('code/rnaseq.R')
#
# msi_report_template_rmarkdown <-
#   file.path(here::here(),"code","msi_classifier.Rmd")
#
# tcga_release <- "release45_20251204"
# gdc_release <- "GDC Data Release 45.0, Dec 4th, 2025"
#
# gdc_projects <- get_gdc_projects()
#
# output_dir <- file.path(
#   here::here(), "output"
# )
#
# data_raw_dir <- file.path(
#   here::here(), "data-raw"
# )
#
# gOncoX <- list()
#
# gOncoX[['basic']] <- geneOncoX::get_basic(
#   cache_dir = data_raw_dir
# )
#
# gOncoX[['basic']]$records <-
#   gOncoX[['basic']]$records |>
#   dplyr::mutate(entrezgene = as.integer(
#     entrezgene
#   ))
#
#
# gOncoX[['gencode']] <- geneOncoX::get_gencode(
#   cache_dir = data_raw_dir
# )
#
# for(b in c('grch37','grch38')){
#   gOncoX[['gencode']]$records[[b]] <-
#     gOncoX[['gencode']]$records[[b]] |>
#     dplyr::mutate(entrezgene = as.integer(
#       entrezgene
#     ))
# }
#
# gOncoX[['alias']] <- geneOncoX::get_alias(
#   cache_dir = data_raw_dir)$records |>
#   dplyr::filter((n_primary_map == 1 &
#                   alias != symbol) |
#                   (alias == symbol & n_primary_map == 1))|>
#   dplyr::select(
#     alias, entrezgene
#   ) |>
#   dplyr::arrange(entrezgene) |>
#   dplyr::mutate(property = "alias") |>
#   dplyr::rename(value = alias) |>
#   dplyr::mutate(entrezgene = as.integer(
#     entrezgene
#   ))
#
#
# #####--- Clinical -----#####
# tcga_clinical_info <- get_tcga_clinical(
#   gdc_projects = gdc_projects,
#   tcga_release = tcga_release,
#   overwrite = F,
#   output_dir = output_dir)
#
# #####--- SNVs/InDels -----#####
# tcga_calls <- get_tcga_snv(
#   gdc_projects = gdc_projects,
#   tcga_clinical_info = tcga_clinical_info,
#   tcga_release = tcga_release,
#   data_raw_dir = data_raw_dir,
#   gdc_files_per_chunk = 20,
#   gOncoX = gOncoX,
#   overwrite = F,
#   output_dir = output_dir)
#
# #####--- VCF output -----#####
# write_tcga_vcf(
#   tcga_calls = tcga_calls,
#   output_dir = output_dir,
#   tcga_release = tcga_release)
#
# #####--- TMB -----#####
# tmb_stats <- calculate_sample_tmb(
#   tcga_calls = tcga_calls,
#   t_depth_min = 30,
#   t_vaf_min = 0.05,
#   tcga_release = tcga_release,
#   tcga_clinical_info = tcga_clinical_info,
#   overwrite = T,
#   output_dir = output_dir)
#
#
# #####---- CNA -------#####
# get_cna_calls(
#     gdc_projects = gdc_projects,
#     tcga_clinical_info = tcga_clinical_info,
#     tcga_release = tcga_release,
#     data_raw_dir = data_raw_dir,
#     gOncoX = gOncoX,
#     overwrite = T,
#     clear_cache = F,
#     output_dir = output_dir)
#
# #####--- Gene aberration rates -----#####
# gene_mutation_stats <- calculate_gene_mutation_rate(
#   tcga_calls = tcga_calls,
#   tcga_release = tcga_release,
#   overwrite = T,
#   tcga_clinical_info = tcga_clinical_info,
#   output_dir = output_dir)
#
#
# #####--- MSI classifier -----#####
#
# ## NOTE!! For some reason (not yet clear), the output object
# ## for the msi classifier ('msi_classifier.rds') gets immensely huge
# ## when running the function from here. When sourcing the internal code
# ## of the function, the output object size gets as expected
# generate_msi_classifier(
#     msi_report_template_rmarkdown = msi_report_template_rmarkdown,
#     tcga_calls = tcga_calls,
#     tcga_clinical_info = tcga_clinical_info,
#     tcga_release = tcga_release,
#     gdc_release = gdc_release,
#     data_raw_dir = data_raw_dir,
#     overwrite = FALSE,
#     output_dir = output_dir)
#
# #####--- RNAseq -----#####
#
# get_bulk_rnaseq(
#   gdc_projects = tail(gdc_projects,6),
#   tcga_clinical_info = tcga_clinical_info,
#   tcga_release = tcga_release,
#   data_raw_dir = data_raw_dir,
#   overwrite = F,
#   gOncoX = gOncoX,
#   output_dir = output_dir)
#
#
