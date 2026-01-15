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
  gdc_release = "release45_20251204",
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
    msi_report_template_rmarkdown = NA,
    t_depth_min = 30,
    t_vaf_min = 0.05,
    gdc_release = NA,
    data_raw_dir = NA,
    overwrite = FALSE,
    output_dir = NA){

  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "msi"))){
    dir.create(
      file.path(output_dir, tcga_release, "msi"))
  }

  msi_classifier_fname = file.path(
    output_dir, tcga_release, "msi", "tcga_msi_classifier2.rds")

  if(file.exists(msi_classifier_fname) & overwrite == F){
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

  ## 1. Filter calls based on DP and VAF
  ## - t_depth_min: minimum tumor depth
  ## - t_vaf_min: minimum tumor variant allelic fraction
  ## 2. Only keep samples for which we have gold standard MSI data
  ## 3. Limit samples to those with a minimum of n = 50 SNV/indel
  ##    calls after DP/AF filtering
  snv_indel_calls_filtered <- snv_indel_calls |>
    dplyr::mutate(
      t_vaf = as.numeric(t_alt_count) / as.numeric(t_depth)
    ) |>
    dplyr::filter(
      as.numeric(t_depth) >= t_depth_min &
        as.numeric(t_vaf) >= t_vaf_min
    ) |>
    dplyr::inner_join(
      dplyr::select(
        msi_data_goldstandard,
        tumor_sample_barcode
      ), by = "tumor_sample_barcode"
    )

  sample_call_counts <- as.data.frame(
    dplyr::group_by(
      snv_indel_calls_filtered,
      tumor_sample_barcode) |>
      dplyr::summarise(
        n_calls = dplyr::n(),
        .groups = "drop"
      )
  )

  snv_indel_calls_filtered <- snv_indel_calls_filtered |>
    dplyr::inner_join(
      dplyr::filter(
        sample_call_counts,
        n_calls >= 50) |>
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

  msi_model <- list()
  msi_model$fitted_model <- modfit_rf
  msi_model$variable_importance <- varImp(modfit_rf)
  msi_model$confusion_matrix <- confusionMatrix(
    predict(
      modfit_rf, dplyr::select(testing,-msi_status)),
    as.factor(testing$msi_status))
  msi_model$sample_features <- msi_predmodel_data
  msi_model$sample_calls <- msi_data$calls
  msi_model$plots <- msi_plots
  msi_model$gdc_release <- gdc_release
  msi_model$t_depth_min <- t_depth_min
  msi_model$t_vaf_min <- t_vaf_min
  msi_model$n_test <- nrow(testing)
  msi_model$n_training <- nrow(training)
  msi_model$n_total <- nrow(training) + nrow(testing)
  msi_model$n_COAD <-
    msi_pred_features_response |>
    dplyr::filter(tumor == 'COAD') |> nrow()
  msi_model$n_STAD <-
    msi_pred_features_response |>
    dplyr::filter(tumor == 'STAD') |> nrow()
  msi_model$n_READ <-
    msi_pred_features_response |>
    dplyr::filter(tumor == 'READ') |> nrow()
  msi_model$n_UCEC <-
    msi_pred_features_response |>
    dplyr::filter(tumor == 'UCEC') |> nrow()

  msi_classifier <- list()
  msi_classifier$model <- modfit_rf
  msi_classifier$confMatrix <- msi_model$confusion_matrix
  msi_classifier$training_dataset <- msi_pred_features_response

  saveRDS(
    msi_classifier,
    file = file.path(
      output_dir, gdc_release, "msi", "tcga_msi_classifier.rds")
  )

  saveRDS(
    msi_model,
    file = file.path(
      output_dir, gdc_release, "msi", "tcga_msi_model.rds")
  )

  quarto::quarto_render(
    input = msi_report_template_rmarkdown,
    output_file =
      "tcga_msi_classifier.html")

  system(paste0("mv code/tcga_msi_classifier.html ",
                file.path(
                  output_dir, gdc_release, "msi",
                  "tcga_msi_classifier.html")))

}
