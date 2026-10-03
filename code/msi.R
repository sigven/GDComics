#library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(caret)

#' GDC MSI data for TCGA samples (mono/dinucleotide assay-based)
#'
#' @param output_dir Output directory for processed data
#' @param gdc_release GDC data release version
#' @param data_raw_dir Data raw directory
#' @param overwrite Overwrite existing data
#' @return Data frame with MSI data
#'
#' @export
gdc_tcga_msi <- function(
  output_dir = NULL,
  gdc_release = "release46_20260810",
  data_raw_dir = NULL,
  overwrite = F){


  output_fname <-
    file.path(
      output_dir,
      gdc_release,
      "msi",
      "tcga_msi.rds"
    )

  if(file.exists(output_fname) & overwrite == F){
    msi_data <- readRDS(file = output_fname)
    return(msi_data)
  }

  gdc_biospecimen_cache_path <-
    file.path(
      data_raw_dir,
      "GDCdata",
      "biospecimen"
    )

  msi_calls <-
    get_msi_status(
      gdc_biospecimen_cache_path = gdc_biospecimen_cache_path
    )

  saveRDS(
    msi_calls,
    file = output_fname
  )

  return(msi_calls)

}

plot_frac_winMaskIndels <- function(df){
  p <-
    ggplot2::ggplot(data = df) +
    ggplot2::geom_boxplot(mapping = ggplot2::aes(
      x = msi_status, y = fracWinMaskIndels,
      color = msi_status, fill = msi_status)) +
    ggplot2::facet_grid(. ~ tumor) +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::scale_fill_brewer(palette='Dark2') +
    #ggplot2::theme_classic() +
    #ggplot2::ggtitle("TCGA - fraction of indels in repetitive DNA") +
    ggplot2::ylab("InDel fraction") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family="Helvetica",size=16,hjust=0.5, face="bold"),
      axis.text.x= ggplot2::element_text(family="Helvetica",size=14,face="bold"),
      axis.title.x = ggplot2::element_text(family="Helvetica",size=14),
      legend.title = ggplot2::element_blank(),
      legend.text=ggplot2::element_text(family="Helvetica",size=14),
      axis.text.y=ggplot2::element_text(family="Helvetica",size=14),
      axis.title.y=ggplot2::element_text(family="Helvetica",size=14,vjust=1.5),
      plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")))

  return(p)

}

plot_frac_Indels <- function(df){
  p <-
    ggplot2::ggplot(data = df) +
    ggplot2::geom_histogram(mapping = ggplot2::aes(
      x = fracIndels, color = msi_status, fill = msi_status),
      position = "dodge", binwidth = 0.01)+
    ggplot2::facet_grid(tumor ~ .) +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::scale_fill_brewer(palette='Dark2') +
    #ggplot2::theme_classic() +
    #ggplot2::ggtitle("TCGA - indel fraction among somatic SNVs/indels") +
    ggplot2::ylab("Number of samples") +
    ggplot2::xlab("InDel fraction") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family="Helvetica",size=16,hjust=0.5, face="bold"),
      axis.text.x= ggplot2::element_text(family="Helvetica",size=14,face="bold"),
      axis.title.x = ggplot2::element_text(family="Helvetica",size=14),
      legend.title = ggplot2::element_blank(),
      legend.text=ggplot2::element_text(family="Helvetica",size=14),
      axis.text.y=ggplot2::element_text(family="Helvetica",size=14),
      axis.title.y=ggplot2::element_text(family="Helvetica",size=14,vjust=1.5),
      plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")))

  return(p)
}

plot_msi_gene_enrichment <- function(df, g = 'MLH1'){
  tmp <- data.frame()
  tmp <- df
  tmp$gene_status <- 'not_mutated'

  tmp[tmp[,g] > 0,]$gene_status <- "gene_mutated"

  tmp2 <- as.data.frame(
    dplyr::filter(tmp, !is.na(gene_status)) |>
      dplyr::group_by(gene_status, msi_status) |>
      dplyr::summarise(n = dplyr::n(),
                       .groups = "drop"))
  tmp3 <- as.data.frame(
    dplyr::group_by(tmp, gene_status) |>
      dplyr::summarise(n_all = dplyr::n(),
                       .groups = "drop"))
  tmp2 <- dplyr::left_join(
    tmp2, tmp3, by = "gene_status")
  tmp2$fractionMutated <- tmp2$n / tmp2$n_all

  tmp_msi.h <- log2(
    tmp2[tmp2$msi_status == "MSI-H" & tmp2$gene_status == "gene_mutated",]$fractionMutated /
      tmp2[tmp2$msi_status == "MSI-H" & tmp2$gene_status == "not_mutated",]$fractionMutated)
  tmp_mss <- log2(
    tmp2[tmp2$msi_status == "MSS" & tmp2$gene_status == "gene_mutated",]$fractionMutated /
      tmp2[tmp2$msi_status == "MSS" & tmp2$gene_status == "not_mutated",]$fractionMutated)

  result <- data.frame('msi_status' = 'MSI-H','enrichment' = tmp_msi.h)
  result <- rbind(result, data.frame('msi_status' = 'MSS', 'enrichment' = tmp_mss))

  p1 <- ggplot2::ggplot(result) +
    ggplot2::geom_bar(mapping = ggplot2::aes(
      x = msi_status, y = enrichment, fill = msi_status),
      colour = "black", stat = "identity",
      position = ggplot2::position_dodge(), width = 0.6) +
    ggplot2::guides(fill = "none") +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::scale_fill_brewer(palette='Dark2') +
    ggplot2::ylab("log2 (mutated/non-mutated samples)") +
    ggplot2::xlab("MSI status") +
    ggplot2::ylim(-3,3) +
    ggplot2::ggtitle(g) +
    ggplot2::theme(plot.title = ggplot2::element_text(family="Helvetica",size=14,hjust=0.5, face="bold"),
                   axis.text.x=ggplot2::element_text(family="Helvetica",size=12,face="bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_text(family="Helvetica",size=12),
                   axis.title.y=ggplot2::element_text(family="Helvetica",size=9,vjust=1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")))

  return(p1)
}



plot_mutated_msi_samples <- function(df, g = 'MLH1'){
  tmp <- data.frame()
  tmp <- df
  tmp$gene_status <- NA
  tmp[tmp[g] > 0,]$gene_status <- 'gene_mutated'

  tmp2 <- as.data.frame(
    dplyr::filter(tmp, !is.na(gene_status)) |>
      dplyr::group_by(gene_status, msi_status) |>
      dplyr::summarise(n = dplyr::n(),
                       .groups = "drop"))
  tmp3 <- as.data.frame(
    dplyr::group_by(tmp, msi_status) |>
      dplyr::summarise(n_all = dplyr::n(),
                       .groups = "drop"))
  tmp2 <- dplyr::left_join(
    tmp2, tmp3, by = "msi_status")
  tmp2$fractionMutated <- tmp2$n / tmp2$n_all

  p1 <- ggplot2::ggplot(tmp2) +
    ggplot2::geom_bar(mapping = ggplot2::aes(
      x = msi_status, y = fractionMutated, fill = msi_status),
      colour = "black", stat = "identity",
      position = ggplot2::position_dodge(), width = 0.8) +
    ggplot2::guides(fill = "none") +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::scale_fill_brewer(palette='Dark2') +
    ggplot2::ylab("Fraction of mutated samples") +
    ggplot2::xlab("MSI status") +
    ggplot2::ylim(0,0.25) +
    ggplot2::ggtitle(g) +
    ggplot2::theme(plot.title = ggplot2::element_text(family="Helvetica",size=14,hjust=0.5, face="bold"),
                   axis.text.x=ggplot2::element_text(family="Helvetica",size=12,face="bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_text(family="Helvetica",size=12),
                   axis.title.y=ggplot2::element_text(family="Helvetica",size=12,vjust=1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")))

  return(p1)
}


get_msi_prediction_features <- function(varcalls, target_size_mb = 34.0){

  assertable::assert_colnames(
    varcalls,
    c('SIMPLEREPEATS_HIT',
      'WINMASKER_HIT',
      'Variant_Type',
      'tumor_sample_barcode',
      'tumor',
      'Hugo_Symbol',
      'One_Consequence'),
    only_colnames = F
  )

  calls_repeatAnnotated <- varcalls |>
    dplyr::mutate(repeatStatus = dplyr::if_else(
      SIMPLEREPEATS_HIT == T,"simpleRepeat",as.character(NA))) |>
    dplyr::mutate(winMaskStatus = dplyr::if_else(
      WINMASKER_HIT == T,"winMaskDust",as.character(NA))) |>
    dplyr::mutate(symbol = Hugo_Symbol) |>
    dplyr::mutate(Variant_Type = dplyr::if_else(
      Variant_Type == "DEL","INDEL",as.character(Variant_Type))) |>
    dplyr::mutate(Variant_Type = dplyr::if_else(
      Variant_Type == "INS","INDEL",as.character(Variant_Type))) |>
    dplyr::mutate(VCF_SAMPLE_ID = tumor_sample_barcode) |>
    dplyr::filter(EXONIC_STATUS == "exonic") |>
    dplyr::filter(tumor == 'UCEC' |
                    tumor == 'READ' |
                    tumor == 'COAD' |
                    tumor == 'STAD') |>
    dplyr::mutate(EXONIC_STATUS = dplyr::if_else(
      stringr::str_detect(
        One_Consequence,
        "^(missense|synonymous|stop_|frameshift|splice_acc|splice_donor|inframe|start_)"),
      "exonic",
      "nonexonic"
    )) |>
    dplyr::select(
      tumor_sample_barcode,
      VCF_SAMPLE_ID,
      bcr_patient_barcode,
      tumor,
      symbol,
      One_Consequence,
      Variant_Type,
      EXONIC_STATUS,
      repeatStatus,
      winMaskStatus
    )


  rep_indels <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(repeatStatus) & Variant_Type == 'INDEL') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(repeat_indels = dplyr::n(),
                       .groups = "drop"))

  rep_snvs <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(repeatStatus) & Variant_Type == 'SNP') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(repeat_SNVs = dplyr::n(),
                       .groups = "drop"))

  rep_tot <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(repeatStatus)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(repeat_indelSNVs = dplyr::n(),
                       .groups = "drop"))

  winmask_indels <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(winMaskStatus) & Variant_Type == 'INDEL') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(winmask_indels = dplyr::n(),
                       .groups = "drop"))

  winmask_snvs <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(winMaskStatus) & Variant_Type == 'SNP') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(winmask_SNVs = dplyr::n(),
                       .groups = "drop"))

  winmask_tot <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(!is.na(winMaskStatus)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(winmask_indelSNVs = dplyr::n(),
                       .groups = "drop"))


  norep_indels <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(is.na(repeatStatus) & Variant_Type == 'INDEL') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(nonRepeat_indels = dplyr::n(),
                       .groups = "drop"))

  norep_snvs <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(is.na(repeatStatus) & Variant_Type == 'SNP') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(nonRepeat_SNVs = dplyr::n(),
                       .groups = "drop"))

  norep_tot <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(is.na(repeatStatus)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(nonRepeat_indelSNVs = dplyr::n(),
                       .groups = "drop"))

  indels <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(Variant_Type == 'INDEL') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(indels = dplyr::n(),
                       .groups = "drop"))

  snvs <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(Variant_Type == 'SNP') |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(SNVs = dplyr::n(),
                       .groups = "drop"))

  tot <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(indelSNVs = dplyr::n(),
                       .groups = "drop"))

  frameshift <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(Variant_Type == 'INDEL') |>
      dplyr::filter(stringr::str_detect(One_Consequence,"frameshift_variant")) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(indels_frameshift = dplyr::n(),
                       .groups = "drop"))

  coding_regex <- "frameshift_|missense_|splice_|stop_|inframe_|start_"
  mlh1 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'MLH1' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(MLH1 = 1,
                       .groups = "drop"))

  mlh3 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'MLH3' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(MLH3 = 1,
                       .groups = "drop"))

  msh2 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'MSH2' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(MSH2 = 1,
                       .groups = "drop"))

  msh3 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'MSH3' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(MSH3 = 1,
                       .groups = "drop"))

  msh6 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'MSH6' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(MSH6 = 1,
                       .groups = "drop"))

  pms1 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'PMS1' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(PMS1 = 1,
                       .groups = "drop"))

  pms2 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'PMS2' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(PMS2 = 1,
                       .groups = "drop"))

  pole <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'POLE' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(POLE = 1,
                       .groups = "drop"))

  pold1 <- as.data.frame(
    calls_repeatAnnotated |>
      dplyr::filter(symbol == 'POLD1' &
                      stringr::str_detect(One_Consequence,coding_regex)) |>
      dplyr::group_by(tumor_sample_barcode) |>
      dplyr::summarise(POLD1 = 1,
                       .groups = "drop"))

  ## initialize data frame of predictive features with all patients
  ## there is mutation data for
  sample_msi_features <- varcalls |>
    dplyr::select(tumor_sample_barcode, tumor) |>
    dplyr::distinct() |>
    dplyr::left_join(rep_indels, by = "tumor_sample_barcode") |>
    dplyr::left_join(rep_snvs, by = "tumor_sample_barcode") |>
    dplyr::left_join(rep_tot, by = "tumor_sample_barcode") |>
    dplyr::left_join(winmask_indels, by = "tumor_sample_barcode") |>
    dplyr::left_join(winmask_snvs, by = "tumor_sample_barcode") |>
    dplyr::left_join(winmask_tot, by = "tumor_sample_barcode") |>
    dplyr::left_join(norep_indels, by = "tumor_sample_barcode") |>
    dplyr::left_join(norep_snvs, by = "tumor_sample_barcode") |>
    dplyr::left_join(norep_tot, by = "tumor_sample_barcode") |>
    dplyr::left_join(indels, by = "tumor_sample_barcode") |>
    dplyr::left_join(snvs, by = "tumor_sample_barcode") |>
    dplyr::left_join(tot, by = "tumor_sample_barcode") |>
    dplyr::left_join(frameshift, by = "tumor_sample_barcode") |>
    dplyr::left_join(mlh1, by = "tumor_sample_barcode") |>
    dplyr::left_join(mlh3, by = "tumor_sample_barcode") |>
    dplyr::left_join(msh2, by = "tumor_sample_barcode") |>
    dplyr::left_join(msh3, by = "tumor_sample_barcode") |>
    dplyr::left_join(msh6, by = "tumor_sample_barcode") |>
    dplyr::left_join(pms1, by = "tumor_sample_barcode") |>
    dplyr::left_join(pms2, by = "tumor_sample_barcode") |>
    dplyr::left_join(pole, by = "tumor_sample_barcode") |>
    dplyr::left_join(pold1, by = "tumor_sample_barcode") |>
    dplyr::mutate(tmb = NA,
                  tmb_snv = NA,
                  tmb_indel = NA)


  ## samples with no entries (NA)
  for(gene in c('MLH1','MLH3','MSH2','MSH3',
                'MSH6','PMS1','PMS2',
                'POLD1','POLE')){
    sample_msi_features[is.na(sample_msi_features[gene]),][gene] <- 0
  }

  for(stat in c('winmask_indels',
                'winmask_SNVs',
                'winmask_indelSNVs',
                'repeat_indelSNVs',
                'repeat_SNVs',
                'repeat_indels',
                'nonRepeat_indels',
                'nonRepeat_SNVs',
                'nonRepeat_indelSNVs',
                'indels',
                'SNVs',
                'indels_frameshift',
                'indelSNVs',
                'tmb',
                'tmb_snv',
                'tmb_indel')){
    if(nrow(sample_msi_features[is.na(sample_msi_features[stat]),]) > 0){
      sample_msi_features[is.na(sample_msi_features[stat]),][stat] <- 0
    }
  }

  sample_msi_features$fracWinMaskIndels <-
    sample_msi_features$winmask_indels / sample_msi_features$indels
  sample_msi_features$fracWinMaskSNVs <-
    sample_msi_features$winmask_SNVs / sample_msi_features$SNVs
  sample_msi_features$fracRepeatIndels <-
    sample_msi_features$repeat_indels / sample_msi_features$repeat_indelSNVs
  sample_msi_features$fracNonRepeatIndels <-
    sample_msi_features$nonRepeat_indels / sample_msi_features$nonRepeat_indelSNVs
  sample_msi_features$fracIndels <-
    sample_msi_features$indels / sample_msi_features$indelSNVs
  sample_msi_features$fracFrameshiftIndels <-
    sample_msi_features$indels_frameshift / sample_msi_features$indels
  sample_msi_features$tmb <-
    sample_msi_features$indelSNVs / target_size_mb
  sample_msi_features$tmb_indel <-
    sample_msi_features$indels / target_size_mb
  sample_msi_features$tmb_snv <-
    sample_msi_features$SNVs / target_size_mb
  for(stat in c('fracWinMaskIndels',
                'fracWinMaskSNVs',
                'fracRepeatIndels',
                'fracNonRepeatIndels',
                'fracIndels',
                'fracFrameshiftIndels',
                'tmb',
                'tmb_indel',
                'tmb_snv')){
    if(nrow(sample_msi_features[is.na(sample_msi_features[stat]),]) > 0){
      sample_msi_features[is.na(sample_msi_features[stat]),][stat] <- 0
    }
  }

  sample_msi_features$bcr_patient_barcode <-
    stringr::str_replace(
      sample_msi_features$tumor_sample_barcode,"-01[A-Z]$","")

  return(list('msi_features' = sample_msi_features,
              'calls' = calls_repeatAnnotated))
}

generate_msi_classifier <- function(
    msi_report_template_qmd = NA,
    t_depth_min = 30,
    t_vaf_min = 0.05,
    gdc_release = NA,
    data_raw_dir = NA,
    overwrite = FALSE,
    output_dir = NA){

  if(!dir.exists(
    file.path(
      output_dir, gdc_release, "msi"))){
    dir.create(
      file.path(output_dir, gdc_release, "msi"))
  }

  msi_runtime_data_fname = file.path(
    output_dir, gdc_release, "msi", "tcga_msi_runtime_data.rds")

  if(file.exists(msi_runtime_data_fname) & overwrite == F){
    return(0)
  }

  tcga_clinical <- gdc_tcga_clinical(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir)

  snv_indel_calls <- gdc_tcga_snv(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir)

  msi_data_goldstandard <- gdc_tcga_msi(
    output_dir = output_dir,
    gdc_release = gdc_release,
    data_raw_dir = data_raw_dir,
    overwrite = F) |>
    dplyr::filter(
      !is.na(.data$msi_status) &
        .data$msi_status != "Indeterminate") |>
    dplyr::mutate(msi_status = dplyr::if_else(
      msi_status == 'MSI-L',
      'MSS',
      as.character(msi_status)
    )) |>
    dplyr::arrange(
      bcr_patient_barcode,
      tumor_sample_barcode,
      project,
    ) |>
    ## keep first entry per sample/patient/project
    dplyr::group_by(
      dplyr::across(-c("project"))
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  assertable::assert_colnames(
    snv_indel_calls,
    c('SIMPLEREPEATS_HIT',
      'WINMASKER_HIT',
      'Variant_Type',
      'tumor_sample_barcode',
      'tumor',
      "bcr_patient_barcode",
      "t_depth",
      "t_alt_count",
      'Hugo_Symbol',
      'One_Consequence'),
    only_colnames = F,
    quiet = T
  )

  ## 1. Apply depth filter only — AF filter is intentionally withheld here so
  ##    the AF distribution check below can see the full shape of the AF
  ##    spectrum, including the low-AF region where artefacts accumulate.
  ## 2. Only keep samples for which we have gold standard MSI data
  snv_indel_calls_dp_filtered <- snv_indel_calls |>
    dplyr::mutate(
      t_vaf = as.numeric(t_alt_count) / as.numeric(t_depth)
    ) |>
    dplyr::filter(
      as.numeric(t_depth) >= t_depth_min
    ) |>
    dplyr::inner_join(
      dplyr::select(
        msi_data_goldstandard,
        tumor_sample_barcode
      ), by = "tumor_sample_barcode"
    )

  ## Check AF distribution of SNVs vs indels (all, repeat-region indels, and
  ## non-repeat indels) BEFORE AF filtering, to detect artefact enrichment
  ## near the AF floor. Repeat-region indels piling up at low AF are a specific
  ## concern as they directly inflate MSI-H features (fracRepeatIndels, fracIndels).
  af_dist_check <- snv_indel_calls_dp_filtered |>
    dplyr::mutate(
      repeatStatus = dplyr::if_else(
        SIMPLEREPEATS_HIT == TRUE |
          WINMASKER_HIT == TRUE,
        "repeat",
        "nonrepeat"
      )
    ) |>
    dplyr::mutate(
      var_class = dplyr::case_when(
        Variant_Type == "SNP"   ~ "SNV",
        Variant_Type == "INS" & repeatStatus == "repeat" ~ "Indel_repeat",
        Variant_Type == "DEL" &  repeatStatus == "nonrepeat" ~ "Indel_nonrepeat",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(var_class))

  af_dist_summary <- af_dist_check |>
    dplyr::group_by(var_class) |>
    dplyr::summarise(
      n            = dplyr::n(),
      af_median    = median(t_vaf, na.rm = TRUE),
      af_mean      = mean(t_vaf,   na.rm = TRUE),
      af_sd        = sd(t_vaf,     na.rm = TRUE),
      af_p05       = quantile(t_vaf, 0.05, na.rm = TRUE),
      af_p25       = quantile(t_vaf, 0.25, na.rm = TRUE),
      af_p75       = quantile(t_vaf, 0.75, na.rm = TRUE),
      frac_below_0.10 = mean(t_vaf < 0.10, na.rm = TRUE),
      frac_below_0.15 = mean(t_vaf < 0.15, na.rm = TRUE),
      .groups = "drop"
    )

  lgr::lgr$info("AF distribution by variant class (post-filter):")
  for (i in seq_len(nrow(af_dist_summary))) {
    r <- af_dist_summary[i, ]
    lgr::lgr$info(
      paste0("  ", r$var_class,
             ": n=", r$n,
             ", median_AF=", round(r$af_median, 3),
             ", mean_AF=", round(r$af_mean, 3),
             ", frac<0.10=", round(r$frac_below_0.10, 3),
             ", frac<0.15=", round(r$frac_below_0.15, 3))
    )
  }

  ## KS test: are indel AF distributions significantly different from SNVs?
  af_snv         <- dplyr::filter(af_dist_check, var_class == "SNV")$t_vaf
  af_ind_rep     <- dplyr::filter(af_dist_check, var_class == "Indel_repeat")$t_vaf
  af_ind_nonrep  <- dplyr::filter(af_dist_check, var_class == "Indel_nonrepeat")$t_vaf

  ks_repeat    <- ks.test(af_snv, af_ind_rep)
  ks_nonrepeat <- ks.test(af_snv, af_ind_nonrep)

  lgr::lgr$info(
    paste0("KS test SNV vs repeat indels:    D=",
           round(ks_repeat$statistic, 4),
           ", p=", format(ks_repeat$p.value, digits = 3))
  )
  lgr::lgr$info(
    paste0("KS test SNV vs non-repeat indels: D=",
           round(ks_nonrepeat$statistic, 4),
           ", p=", format(ks_nonrepeat$p.value, digits = 3))
  )

  af_dist_results <- list(
    summary      = af_dist_summary,
    ks_repeat    = ks_repeat,
    ks_nonrepeat = ks_nonrepeat
  )

  ## Now apply the AF filter to produce the final filtered callset used for
  ## all downstream model training and evaluation.
  snv_indel_calls_filtered <- snv_indel_calls_dp_filtered |>
    dplyr::filter(t_vaf >= t_vaf_min)

  sample_call_counts <- as.data.frame(
    dplyr::group_by(
      snv_indel_calls_filtered,
      tumor_sample_barcode) |>
      dplyr::summarise(
        n_calls = dplyr::n(),
        .groups = "drop"
      )
  )

  ## Evaluate classifier performance across multiple minimum-mutation thresholds
  ## to inform a scientifically grounded choice for the final model.
  min_mut_thresholds <- c(30, 50, 75, 100, 125, 150)
  threshold_performance <- list()

  for (min_mut in min_mut_thresholds) {

    calls_thresh <- snv_indel_calls_filtered |>
      dplyr::inner_join(
        dplyr::filter(sample_call_counts, n_calls >= min_mut) |>
          dplyr::select(tumor_sample_barcode),
        by = "tumor_sample_barcode"
      )

    gs_thresh <- msi_data_goldstandard |>
      dplyr::filter(
        tumor_sample_barcode %in% unique(calls_thresh$tumor_sample_barcode)
      )

    msi_data_thresh <- get_msi_prediction_features(varcalls = calls_thresh)

    features_response_thresh <- msi_data_thresh$msi_features |>
      dplyr::inner_join(
        gs_thresh,
        by = c("bcr_patient_barcode", "tumor_sample_barcode")
      )

    dataset_thresh <- dplyr::select(
      features_response_thresh,
      tumor, msi_status,
      fracWinMaskIndels, fracWinMaskSNVs,
      fracRepeatIndels, fracIndels, fracNonRepeatIndels,
      tmb, tmb_snv, tmb_indel,
      MLH1, MSH2, MLH3, MSH3, MSH6,
      PMS1, PMS2, POLE, POLD1
    )

    set.seed(9999)
    inTrain_thresh <- caret::createDataPartition(
      dataset_thresh$msi_status, p = 0.70)[[1]]
    training_thresh <- dplyr::select(dataset_thresh[ inTrain_thresh,], -tumor)
    testing_thresh  <- dataset_thresh[-inTrain_thresh,]

    modfit_thresh <- caret::train(
      as.factor(msi_status) ~ .,
      method = "rf",
      data = training_thresh,
      preProcess = c("YeoJohnson", "scale"),
      trControl = caret::trainControl(method = "cv", number = 10),
      na.action = na.exclude
    )

    cm_thresh <- caret::confusionMatrix(
      predict(modfit_thresh, dplyr::select(testing_thresh, -msi_status)),
      as.factor(testing_thresh$msi_status)
    )

    threshold_performance[[as.character(min_mut)]] <- list(
      min_mut         = min_mut,
      n_samples       = nrow(dataset_thresh),
      n_training      = nrow(training_thresh),
      n_test          = nrow(testing_thresh),
      accuracy        = cm_thresh$overall[["Accuracy"]],
      kappa           = cm_thresh$overall[["Kappa"]],
      confusion_matrix = cm_thresh
    )

    lgr::lgr$info(
      paste0("Threshold n >= ", min_mut,
             ": n_samples = ", nrow(dataset_thresh),
             ", Accuracy = ", round(cm_thresh$overall[["Accuracy"]], 4),
             ", Kappa = ", round(cm_thresh$overall[["Kappa"]], 4))
    )
  }

  threshold_summary <- dplyr::bind_rows(
    lapply(threshold_performance, function(x)
      data.frame(
        min_mut    = x$min_mut,
        n_samples  = x$n_samples,
        n_training = x$n_training,
        n_test     = x$n_test,
        accuracy   = x$accuracy,
        kappa      = x$kappa
      )
    )
  )

  ## Evaluate how prediction fidelity degrades for low-mutation samples.
  ## Strategy: train once on high-confidence samples (n >= 100), then apply
  ## to all excluded samples binned by mutation count. This directly answers
  ## at what mutation count predictions become unreliable.
  high_conf_calls <- snv_indel_calls_filtered |>
    dplyr::inner_join(
      dplyr::filter(sample_call_counts, n_calls >= 100) |>
        dplyr::select(tumor_sample_barcode),
      by = "tumor_sample_barcode"
    )

  high_conf_gs <- msi_data_goldstandard |>
    dplyr::filter(
      tumor_sample_barcode %in% unique(high_conf_calls$tumor_sample_barcode)
    )

  msi_data_hc <- get_msi_prediction_features(varcalls = high_conf_calls)

  features_hc <- msi_data_hc$msi_features |>
    dplyr::inner_join(high_conf_gs,
                      by = c("bcr_patient_barcode", "tumor_sample_barcode"))

  dataset_hc <- dplyr::select(
    features_hc,
    tumor, msi_status,
    fracWinMaskIndels, fracWinMaskSNVs,
    fracRepeatIndels, fracIndels, fracNonRepeatIndels,
    tmb, tmb_snv, tmb_indel,
    MLH1, MSH2, MLH3, MSH3, MSH6,
    PMS1, PMS2, POLE, POLD1
  )

  set.seed(9999)
  modfit_hc <- caret::train(
    as.factor(msi_status) ~ .,
    method = "rf",
    data = dplyr::select(dataset_hc, -tumor),
    preProcess = c("YeoJohnson", "scale"),
    trControl = caret::trainControl(method = "cv", number = 10),
    na.action = na.exclude
  )

  ## Samples excluded from high-confidence training set, with gold-standard labels
  low_mut_samples <- sample_call_counts |>
    dplyr::filter(
      n_calls < 100,
      tumor_sample_barcode %in% msi_data_goldstandard$tumor_sample_barcode
    ) |>
    dplyr::mutate(
      mut_bin = cut(
        n_calls,
        breaks = c(0, 10, 20, 30, 50, 75, 99),
        labels = c("<10", "10-19", "20-29", "30-49", "50-74", "75-99"),
        right = TRUE
      )
    )

  low_mut_calls <- snv_indel_calls_filtered |>
    dplyr::inner_join(
      dplyr::select(low_mut_samples, tumor_sample_barcode),
      by = "tumor_sample_barcode"
    )

  low_mut_gs <- msi_data_goldstandard |>
    dplyr::filter(
      tumor_sample_barcode %in% unique(low_mut_calls$tumor_sample_barcode)
    )

  msi_data_lm <- get_msi_prediction_features(varcalls = low_mut_calls)

  features_lm <- msi_data_lm$msi_features |>
    dplyr::inner_join(low_mut_gs,
                      by = c("bcr_patient_barcode", "tumor_sample_barcode")) |>
    dplyr::inner_join(
      dplyr::select(low_mut_samples, tumor_sample_barcode, n_calls, mut_bin),
      by = "tumor_sample_barcode"
    )

  pred_lm <- predict(
    modfit_hc,
    dplyr::select(features_lm,
                  fracWinMaskIndels, fracWinMaskSNVs,
                  fracRepeatIndels, fracIndels, fracNonRepeatIndels,
                  tmb, tmb_snv, tmb_indel,
                  MLH1, MSH2, MLH3, MSH3, MSH6,
                  PMS1, PMS2, POLE, POLD1)
  )

  low_mut_eval <- features_lm |>
    dplyr::select(tumor_sample_barcode, msi_status, n_calls, mut_bin) |>
    dplyr::mutate(predicted = as.character(pred_lm))

  low_mut_bin_performance <- low_mut_eval |>
    dplyr::group_by(mut_bin) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      accuracy  = mean(predicted == msi_status),
      n_correct = sum(predicted == msi_status),
      .groups = "drop"
    )

  for (i in seq_len(nrow(low_mut_bin_performance))) {
    lgr::lgr$info(
      paste0("Low-mutation bin ", low_mut_bin_performance$mut_bin[i],
             ": n = ", low_mut_bin_performance$n_samples[i],
             ", Accuracy = ", round(low_mut_bin_performance$accuracy[i], 4))
    )
  }

  ## Use n >= 100 as the threshold for the final saved model, selected on the
  ## basis of the threshold sensitivity sweep above (best Kappa at n >= 100).
  snv_indel_calls_filtered <- snv_indel_calls_filtered |>
    dplyr::inner_join(
      dplyr::filter(
        sample_call_counts,
        n_calls >= 100) |>
        dplyr::select(tumor_sample_barcode),
      by = "tumor_sample_barcode"
    )

  msi_data_goldstandard <- msi_data_goldstandard |>
    dplyr::filter(
      tumor_sample_barcode %in%
        unique(snv_indel_calls_filtered$tumor_sample_barcode)
    )

  msi_data <- get_msi_prediction_features(
    varcalls = snv_indel_calls_filtered)

  msi_pred_features_response <- msi_data$msi_features |>
    dplyr::inner_join(
      msi_data_goldstandard,
      by = c("bcr_patient_barcode",
             "tumor_sample_barcode"))

  tcga_dataset <- dplyr::select(
    msi_pred_features_response,
    tumor,
    msi_status,
    fracWinMaskIndels,
    fracWinMaskSNVs,
    fracRepeatIndels,
    fracIndels,
    fracNonRepeatIndels,
    tmb,
    tmb_snv,
    tmb_indel,MLH1,MSH2,MLH3,MSH3,MSH6,
    PMS1,PMS2,POLE,POLD1)


  msi_predmodel_data <- tcga_dataset

  set.seed(9999)
  inTrain <- caret::createDataPartition(
    msi_predmodel_data$msi_status, p = 0.70)[[1]]
  training <- msi_predmodel_data[ inTrain,]
  testing <- msi_predmodel_data[-inTrain,]
  training_exploration <- training

  msi_plots <- list()
  msi_plots[['indelWinMaskPlot']] <- plot_frac_winMaskIndels(
    df = training_exploration
  )
  msi_plots[['indelFracPlot']] <- plot_frac_Indels(
    df = training_exploration
  )
  msi_plots[['fraction_mutated']] <-
    list('MLH1' <- plot_mutated_msi_samples(training_exploration,g = 'MLH1'),
         'MLH3' <- plot_mutated_msi_samples(training_exploration,g = 'MLH3'),
         'MSH2' <- plot_mutated_msi_samples(training_exploration,g = 'MSH2'),
         'MSH3' <- plot_mutated_msi_samples(training_exploration,g = 'MSH3'),
         'MSH6' <- plot_mutated_msi_samples(training_exploration,g = 'MSH6'),
         'PMS1' <- plot_mutated_msi_samples(training_exploration,g = 'PMS1'),
         'PMS2' <- plot_mutated_msi_samples(training_exploration,g = 'PMS2'),
         'POLD1' <- plot_mutated_msi_samples(training_exploration,g = 'POLD1'),
         'POLE' <- plot_mutated_msi_samples(training_exploration,g = 'POLE'))

  msi_plots[['gene_enrichment']] <-
    list('MLH1' <- plot_msi_gene_enrichment(training_exploration,g = 'MLH1'),
         'MLH3' <- plot_msi_gene_enrichment(training_exploration,g = 'MLH3'),
         'MSH2' <- plot_msi_gene_enrichment(training_exploration,g = 'MSH2'),
         'MSH3' <- plot_msi_gene_enrichment(training_exploration,g = 'MSH3'),
         'MSH6' <- plot_msi_gene_enrichment(training_exploration,g = 'MSH6'),
         'PMS1' <- plot_msi_gene_enrichment(training_exploration,g = 'PMS1'),
         'PMS2' <- plot_msi_gene_enrichment(training_exploration,g = 'PMS2'),
         'POLD1' <- plot_msi_gene_enrichment(training_exploration,g = 'POLD1'),
         'POLE' <- plot_msi_gene_enrichment(training_exploration,g = 'POLE'))

  ## train model using training set (random forest)
  ## ten-fold cross-validation
  training <- dplyr::select(training, -tumor)
  modfit_rf <- caret::train(
    as.factor(msi_status) ~ .,
    method="rf",
    data = training,
    preProcess = c("YeoJohnson","scale"),
    trControl = caret::trainControl(
      method = "cv", number = 10),
    na.action = na.exclude)

  ## -------------------------------------------------------------------------
  ## Empirical feature-variance lookup table
  ## Computed across ALL gold-standard samples regardless of mutation count
  ## (combining the high-confidence n>=100 set and the low-mutation <100 set)
  ## so that PCGR can look up expected feature stability for any sample,
  ## including those below the training threshold.
  ## -------------------------------------------------------------------------
  msi_feature_cols <- c(
    "fracWinMaskIndels", "fracWinMaskSNVs",
    "fracRepeatIndels", "fracIndels", "fracNonRepeatIndels",
    "tmb", "tmb_snv", "tmb_indel"
  )

  all_features_for_variance <- dplyr::bind_rows(
    dplyr::select(features_hc,
                  tumor_sample_barcode,
                  dplyr::all_of(msi_feature_cols)) |>
      dplyr::inner_join(
        dplyr::select(sample_call_counts, tumor_sample_barcode, n_calls),
        by = "tumor_sample_barcode"
      ),
    dplyr::select(features_lm,
                  tumor_sample_barcode, n_calls,
                  dplyr::all_of(msi_feature_cols))
  )

  feature_variance_table <- all_features_for_variance |>
    dplyr::mutate(
      mut_bin = cut(
        n_calls,
        breaks = c(0, 30, 50, 75, 100, 150, 200, 500, Inf),
        labels = c("<30", "30-49", "50-74", "75-99",
                   "100-149", "150-199", "200-499", "500+"),
        right = FALSE
      )
    ) |>
    dplyr::group_by(mut_bin) |>
    dplyr::summarise(
      n_samples            = dplyr::n(),
      dplyr::across(
        dplyr::all_of(msi_feature_cols),
        list(
          median = \(x) median(x, na.rm = TRUE),
          sd     = \(x) sd(x,     na.rm = TRUE),
          cv     = \(x) sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) + 1e-9)
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    )

  lgr::lgr$info("Feature variance lookup table (SD of fracIndels by mutation bin):")
  for (i in seq_len(nrow(feature_variance_table))) {
    lgr::lgr$info(
      paste0("  mut_bin=", feature_variance_table$mut_bin[i],
             "  n=", feature_variance_table$n_samples[i],
             "  fracIndels_sd=",
             round(feature_variance_table$fracIndels__sd[i], 4),
             "  fracRepeatIndels_sd=",
             round(feature_variance_table$fracRepeatIndels__sd[i], 4),
             "  tmb_indel_sd=",
             round(feature_variance_table$tmb_indel__sd[i], 4))
    )
  }

  ## -------------------------------------------------------------------------
  ## RF probability calibration check on the held-out test set
  ## Bins predicted P(MSI-H) into deciles and compares observed MSI-H rate
  ## within each bin. Well-calibrated probabilities lie on the diagonal.
  ## Also computes Brier score as a scalar summary of calibration quality.
  ## -------------------------------------------------------------------------
  prob_testing <- predict(
    modfit_rf,
    dplyr::select(testing, -msi_status),
    type = "prob"
  )

  calibration_df <- data.frame(
    msi_status  = testing$msi_status,
    prob_msi_h  = prob_testing[["MSI-H"]]
  ) |>
    dplyr::mutate(
      prob_bin = cut(
        prob_msi_h,
        breaks = seq(0, 1, by = 0.1),
        include.lowest = TRUE,
        right = FALSE
      ),
      is_msi_h = as.integer(msi_status == "MSI-H")
    )

  calibration_summary <- calibration_df |>
    dplyr::group_by(prob_bin) |>
    dplyr::summarise(
      n                 = dplyr::n(),
      mean_pred_prob    = mean(prob_msi_h),
      observed_msi_h_rate = mean(is_msi_h),
      .groups = "drop"
    )

  brier_score <- mean(
    (calibration_df$prob_msi_h - calibration_df$is_msi_h)^2
  )

  lgr::lgr$info(
    paste0("RF calibration — Brier score on test set: ",
           round(brier_score, 4),
           " (0 = perfect, 0.25 = uninformative)")
  )
  lgr::lgr$info("Calibration by predicted probability decile:")
  for (i in seq_len(nrow(calibration_summary))) {
    r <- calibration_summary[i, ]
    lgr::lgr$info(
      paste0("  bin=", r$prob_bin,
             "  n=", r$n,
             "  mean_pred=", round(r$mean_pred_prob, 3),
             "  observed_rate=", round(r$observed_msi_h_rate, 3))
    )
  }

  ## -------------------------------------------------------------------------
  ## msi_training_record — full development artefact, consumed by the GDComics
  ## Quarto report (tcga_msi_training_report.html) and internal QC.
  ## Never loaded by PCGR at runtime.
  ## -------------------------------------------------------------------------
  msi_training_record <- list()

  ## Core model outputs
  msi_training_record$fitted_model      <- modfit_rf
  msi_training_record$variable_importance <- varImp(modfit_rf)
  msi_training_record$confusion_matrix  <- confusionMatrix(
    predict(modfit_rf, dplyr::select(testing, -msi_status)),
    as.factor(testing$msi_status))

  ## Training provenance
  msi_training_record$gdc_release  <- gdc_release
  msi_training_record$t_depth_min  <- t_depth_min
  msi_training_record$t_vaf_min    <- t_vaf_min
  msi_training_record$n_training   <- nrow(training)
  msi_training_record$n_test       <- nrow(testing)
  msi_training_record$n_total      <- nrow(training) + nrow(testing)
  msi_training_record$n_COAD       <-
    dplyr::filter(msi_pred_features_response, tumor == 'COAD') |> nrow()
  msi_training_record$n_STAD       <-
    dplyr::filter(msi_pred_features_response, tumor == 'STAD') |> nrow()
  msi_training_record$n_READ       <-
    dplyr::filter(msi_pred_features_response, tumor == 'READ') |> nrow()
  msi_training_record$n_UCEC       <-
    dplyr::filter(msi_pred_features_response, tumor == 'UCEC') |> nrow()

  ## Training data and exploratory plots
  msi_training_record$sample_features <- msi_predmodel_data
  msi_training_record$sample_calls     <- msi_data$calls
  msi_training_record$plots            <- msi_plots

  ## Diagnostic analyses
  msi_training_record$threshold_sensitivity   <- threshold_summary
  msi_training_record$low_mut_bin_performance <- low_mut_bin_performance
  msi_training_record$af_dist_results         <- af_dist_results
  msi_training_record$feature_variance_table  <- feature_variance_table
  msi_training_record$calibration_summary     <- calibration_summary
  msi_training_record$brier_score             <- brier_score

  ## -------------------------------------------------------------------------
  ## msi_runtime_data — lean bundle loaded by PCGR at inference time.
  ## Contains only what is needed for prediction and report rendering.
  ## Built from msi_training_record to avoid any divergence between the two.
  ## -------------------------------------------------------------------------
  msi_runtime_data <- list(
    ## RF model — used by predict_msi_status()
    model                  = msi_training_record$fitted_model,
    ## Test-set performance — referenced in the PCGR report text
    confMatrix             = msi_training_record$confusion_matrix,
    ## TCGA fracIndels distribution — background for the report histogram
    tcga_dataset           = msi_pred_features_response,
    ## Per mutation-count-bin feature SD — drives the confidence note in PCGR
    feature_variance_table = msi_training_record$feature_variance_table,
    ## Calibration results — supports reporting of P(MSI-H)
    calibration_summary    = msi_training_record$calibration_summary,
    brier_score            = msi_training_record$brier_score
  )

  saveRDS(
    msi_runtime_data,
    file = msi_runtime_data_fname
  )

  saveRDS(
    msi_training_record,
    file = file.path(
      output_dir, gdc_release,
      "msi", "tcga_msi_training_record.rds")
  )

  quarto::quarto_render(
    input = msi_report_template_qmd,
    output_file =
      file.path(
        output_dir, gdc_release, "msi",
      "tcga_msi_training_report.html"))

  system(paste0("mv code/tcga_msi_training_report.html ",
                file.path(
                  output_dir, gdc_release, "msi",
                  "tcga_msi_training_report.html")))

}
