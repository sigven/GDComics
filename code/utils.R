library(reshape2)


read_gtf <- function(path, zero_based = TRUE) {
  gtf <- rtracklayer::import(path)
  gtf <- as.data.frame(gtf)
  gtf <- dplyr::mutate_if(gtf, is.factor, as.character)
  res <- dplyr::rename(gtf, chrom = seqnames)

  if (zero_based) {
    res <- dplyr::mutate(res, start = start - 1L)
  }

  tibble::as_tibble(res)
}


get_gencode_xref <- function(
    data_raw_dir = NULL,
    update_ensembl_gene = FALSE){

  gencode_xref_fname <-
    file.path(data_raw_dir, "gencode_v36_xref.rds")

  gencode_xref <- list()
  if(!file.exists(gencode_xref_fname) |
     update_ensembl_gene == T){

    gOncoX <- list()

    gOncoX[['basic']] <- geneOncoX::get_basic(
      cache_dir = data_raw_dir
    )

    gOncoX[['basic']]$records <-
      gOncoX[['basic']]$records |>
      dplyr::mutate(entrezgene = as.integer(
        entrezgene
      ))

    gencode_gtf_fname <- file.path(
      data_raw_dir,
      "gencode",
      "gencode.v36.annotation.gtf.gz")

    gencode_xref_df <- as.data.frame(
      read_gtf(path = gencode_gtf_fname) |>
        dplyr::select(gene_id, hgnc_id) |>
        dplyr::mutate(ensembl_gene_id = stringr::str_replace(
          gene_id,"\\.[0-9]{1,}","")) |>
        dplyr::select(ensembl_gene_id, hgnc_id) |>
        dplyr::filter(!is.na(hgnc_id)) |>
        dplyr::filter(!stringr::str_detect(
          ensembl_gene_id,"_PAR_Y")) |>
        dplyr::mutate(hgnc_id = stringr::str_replace(
          hgnc_id, "HGNC:",""
        )) |>
        ## remove ensembl gene id's that are pointing
        ## to the same gene symbol as other
        ## ensembl identifiers
        dplyr::filter(ensembl_gene_id != "ENSG00000145075") |>
        dplyr::filter(ensembl_gene_id != "ENSG00000285437") |>
        dplyr::filter(ensembl_gene_id != "ENSG00000286065") |>
        dplyr::distinct() |>

        ## join with geneOncoX basic records to get
        ## gene symbols, entrezgene, gene biotype info
        ## etc. that are up-to-date
        dplyr::inner_join(
          dplyr::select(gOncoX$basic$records,
                        symbol, entrezgene,
                        name,
                        hgnc_id, gene_biotype),
          by = "hgnc_id"
        ) |>
        dplyr::rename(
          genename = "name",
          biotype = "gene_biotype") |>
        dplyr::select(
          c("ensembl_gene_id",
            "symbol",
            "entrezgene",
            "genename",
            #"hgnc_id",
            "biotype")
        )
    )

    names(gencode_xref_df) <-
      toupper(names(gencode_xref_df))

    ## identify duplicates based on entrezgene
    duplicates <- as.data.frame(
      dplyr::group_by(gencode_xref_df, ENTREZGENE) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::filter(n > 1)
    )

    gencode_xref[['all']] <- gencode_xref_df |>
      dplyr::anti_join(duplicates, by = "ENTREZGENE") |>
      dplyr::arrange(.data$ENTREZGENE) |>
      dplyr::distinct()

    gencode_xref[['protein_coding']] <-
      gencode_xref[['all']] |>
      dplyr::filter(
        .data$BIOTYPE == "protein_coding")

    saveRDS(gencode_xref, file = gencode_xref_fname)
  }else{
    gencode_xref <- readRDS(file = gencode_xref_fname)
  }

  return(gencode_xref)


}


get_gdc_projects <- function(tcga = TRUE){
  gdc_projects <- as.data.frame(
    TCGAbiolinks::getGDCprojects())

  if(tcga == TRUE){
    gdc_projects <- gdc_projects |>
      dplyr::filter(is.na(dbgap_accession_number) &
                      startsWith(id,"TCGA"))
  }
  return(gdc_projects)

}


get_proper_maf_alleles <- function(maf_df, genome_seq, seqinfo){

  maf_SNV <- dplyr::filter(maf_df, Variant_Type == 'SNP' |
                             Variant_Type == 'SNV')
  maf_SNV$REF <- as.character(maf_SNV$Reference_Allele)
  maf_SNV$ALT <- as.character(maf_SNV$Tumor_Seq_Allele2)
  maf_SNV$POS <- maf_SNV$Start_Position
  maf_SNV$CHROM <- stringr::str_replace(maf_SNV$Chromosome,'chr','')
  maf_SNV$end <- maf_SNV$End_Position
  maf_ALL <- maf_SNV


  maf_MNV <- dplyr::filter(
    maf_df, (Variant_Type == "ONP" |
               Variant_Type == "TNP" |
               Variant_Type == "DNP") &
      nchar(Reference_Allele) ==
      nchar(Tumor_Seq_Allele2)
  )
  if(NROW(maf_MNV) > 0){
    maf_MNV <- maf_MNV |>
      dplyr::mutate(
        REF = as.character(Reference_Allele),
        ALT = as.character(Tumor_Seq_Allele2),
        CHROM = stringr::str_replace(
          Chromosome, "chr", ""),
        POS = Start_Position,
        end = End_Position)

    maf_ALL <- dplyr::bind_rows(
      maf_ALL, maf_MNV)
  }

  maf_INS <- dplyr::filter(
    maf_df, Variant_Type == 'INS' |
      Variant_Type == 'insertion')
  maf_DEL <- dplyr::filter(
    maf_df, Variant_Type == 'DEL' |
      Variant_Type == 'deletion')

  if(nrow(maf_DEL) > 0){
    ## get appropriate alleles (VCF-like) of reference and alternate (DELETIONS)
    maf_del_gr <- GenomicRanges::makeGRangesFromDataFrame(
      maf_DEL, keep.extra.columns = T, seqinfo = seqinfo,
      seqnames.field = 'Chromosome',start.field = 'Start_Position',
      end.field = 'End_Position',
      ignore.strand = T, starts.in.df.are.0based = F)

    maf_del_flank_gr <- GenomicRanges::flank(
      maf_del_gr, width = 1, start=T)
    maf_del_flank_seq <- Biostrings::getSeq(genome_seq, maf_del_flank_gr)
    maf_del_seq <- Biostrings::getSeq(genome_seq, maf_del_gr)
    vcf_alleles_alt <- data.frame('ALT'=as.character(
      toupper(unlist(strsplit(toString(maf_del_flank_seq),", ")))),
      stringsAsFactors=F)
    vcf_alleles_ref <- data.frame('REF'=as.character(
      toupper(unlist(strsplit(toString(maf_del_seq),", ")))),
      stringsAsFactors=F)
    vcf_alleles <- cbind(vcf_alleles_ref, vcf_alleles_alt)
    vcf_alleles$REF <- paste0(vcf_alleles$ALT,vcf_alleles$REF)
    maf_DEL <- cbind(maf_DEL,vcf_alleles)
    maf_DEL$POS <- maf_DEL$Start_Position - 1
    maf_DEL$CHROM <- stringr::str_replace(maf_DEL$Chromosome,'chr','')
    maf_DEL$end <- maf_DEL$End_Position

    maf_ALL <- dplyr::bind_rows(maf_ALL, maf_DEL)
  }

  if(nrow(maf_INS) > 0){
    ## get appropriate alleles (VCF-like) of reference and alternate (INSERTIONS)
    maf_ins_gr <- GenomicRanges::makeGRangesFromDataFrame(
      maf_INS, keep.extra.columns = T, seqinfo = seqinfo,
      seqnames.field = 'Chromosome',start.field = 'Start_Position',
      end.field = 'Start_Position', ignore.strand = T,
      starts.in.df.are.0based = F)
    maf_ins_seq <- Biostrings::getSeq(genome_seq, maf_ins_gr)
    vcf_alleles_alt <- data.frame('REF'=as.character(
      toupper(unlist(strsplit(toString(maf_ins_seq),", ")))),
      stringsAsFactors=F)
    maf_INS <- cbind(maf_INS,vcf_alleles_alt)
    maf_INS$ALT <- as.character(
      paste0(maf_INS$REF,maf_INS$Tumor_Seq_Allele2))
    maf_INS$POS <- maf_INS$Start_Position
    maf_INS$CHROM <- stringr::str_replace(maf_INS$Chromosome,'chr','')
    maf_INS$end <- maf_INS$End_Position

    maf_ALL <- dplyr::bind_rows(maf_ALL, maf_INS)
  }

  maf_ALL$GENOMIC_CHANGE <- paste(
    paste(
      paste(
        paste0("g.chr",
               maf_ALL$CHROM),
        maf_ALL$POS,sep=":"),
      maf_ALL$REF,sep=":"),
    maf_ALL$ALT,sep=">")
  maf_ALL <- maf_ALL |> dplyr::select(-end)

  return(maf_ALL)

}


cohort_mutation_stats <- function(calls_df,
                                  clinical_df,
                                  vartype = 'snv_indel',
                                  genomic_strata = "gene",
                                  clinical_strata = "site",
                                  num_digits = 2){

  mutation_stats <- NULL
  if(vartype == 'snv_indel'){
    mutation_status_sample <- calls_df |>
      dplyr::filter(CODING_STATUS == "coding") |>
      dplyr::select(bcr_patient_barcode, symbol, primary_site,
                    primary_diagnosis_very_simplified) |>
      dplyr::filter(!is.na(primary_site)) |>
      dplyr::distinct() |>
      dplyr::mutate(mut_status = "MUT") |>
      dplyr::mutate(genomic_strata = "gene",
                    variant_type = vartype,
                    consensus_calls = F)

    mutation_status_pancancer <- calls_df |>
      dplyr::filter(CODING_STATUS == "coding") |>
      dplyr::select(bcr_patient_barcode, symbol) |>
      dplyr::mutate(primary_site = "Pancancer",
                    primary_diagnosis_very_simplified = "Pancancer") |>
      dplyr::distinct() |>
      dplyr::mutate(mut_status = "MUT") |>
      dplyr::mutate(genomic_strata = "gene",
                    variant_type = vartype,
                    consensus_calls = F)

    mutation_status_sample <- dplyr::bind_rows(
      mutation_status_sample, mutation_status_pancancer)

    sample_numbers_cohort <- as.data.frame(
      dplyr::select(clinical_df, bcr_patient_barcode,
                    primary_site, primary_diagnosis_very_simplified) |>
        dplyr::distinct() |>
        dplyr::select(-bcr_patient_barcode) |>
        dplyr::filter(!is.na(primary_site)) |>
        dplyr::group_by(primary_site) |>
        dplyr::summarise(tot_samples = dplyr::n(),
                         .groups = "drop")
      )

    pancancer_sample_numbers <- as.data.frame(
      dplyr::select(clinical_df, bcr_patient_barcode) |>
      dplyr::distinct() |>
      dplyr::mutate(primary_site = "Pancancer") |>
      dplyr::select(-bcr_patient_barcode) |>
      dplyr::group_by(primary_site) |>
      dplyr::summarise(tot_samples = dplyr::n(), .groups = "drop")
    )

    sample_numbers_cohort <- dplyr::bind_rows(sample_numbers_cohort,
                                              pancancer_sample_numbers)
    sample_numbers_cohort$clinical_strata <- 'site'


    mutation_stats <- as.data.frame(
      dplyr::group_by(mutation_status_sample, symbol,
                      primary_site, genomic_strata,
                      variant_type, consensus_calls) |>
        dplyr::summarise(samples_mutated = dplyr::n(), .groups = "drop") |>
        #dplyr::left_join(sample_numbers_cohort) |>
        dplyr::left_join(sample_numbers_cohort, by = "primary_site") |>
        dplyr::mutate(percent_mutated =
                        round((samples_mutated / tot_samples) * 100,
                              digits = num_digits))
    )

    mutation_stats <- as.data.frame(
      mutation_stats |>
        dplyr::group_by(primary_site, variant_type) |>
        dplyr::mutate(percentile = dplyr::ntile(percent_mutated, 100),
                      decile = dplyr::ntile(percent_mutated,10))
    )

    if(clinical_strata == "site_diagnosis"){
      sample_numbers_cohort <- as.data.frame(
        dplyr::select(clinical_df, bcr_patient_barcode, primary_site,
                      primary_diagnosis_very_simplified) |>
        dplyr::distinct() |>
        dplyr::select(-bcr_patient_barcode) |>
        dplyr::filter(!is.na(primary_site)) |>
        dplyr::group_by(primary_site, primary_diagnosis_very_simplified) |>
        dplyr::summarise(tot_samples = dplyr::n(), .groups = "drop")
      )

      pancancer_sample_numbers <- as.data.frame(
        dplyr::select(clinical_df, bcr_patient_barcode) |>
        dplyr::distinct() |>
        dplyr::mutate(primary_site = "Pancancer",
                      primary_diagnosis_very_simplified = "Pancancer") |>
        dplyr::select(-bcr_patient_barcode) |>
        dplyr::group_by(primary_site, primary_diagnosis_very_simplified) |>
        dplyr::summarise(tot_samples = dplyr::n(), .groups = "drop")
      )

      sample_numbers_cohort <- dplyr::bind_rows(
        sample_numbers_cohort,pancancer_sample_numbers)
      sample_numbers_cohort$clinical_strata <- 'site_diagnosis'

      mutation_stats <- as.data.frame(
        dplyr::group_by(mutation_status_sample,
                        symbol, primary_site,
                        primary_diagnosis_very_simplified,
                        genomic_strata,
                        variant_type, consensus_calls) |>
          dplyr::summarise(samples_mutated = dplyr::n(), .groups = "drop") |>
          dplyr::left_join(sample_numbers_cohort,
                           by = c("primary_site",
                                  "primary_diagnosis_very_simplified")) |>
          dplyr::mutate(percent_mutated = round(
            (samples_mutated / tot_samples) * 100, digits = num_digits))
      )

      mutation_stats <- as.data.frame(
        mutation_stats |>
          dplyr::group_by(
            primary_site, primary_diagnosis_very_simplified,
            variant_type) |>
          dplyr::mutate(percentile = dplyr::ntile(
            percent_mutated, 100),
            decile = dplyr::ntile(percent_mutated,10))
      )
    }
  }else{
    if(vartype == "cna_homdel" | vartype == "cna_ampl"){
      mutation_status_sample <- calls_df |>
        dplyr::filter(mut_status == "HOMDEL") |>
        dplyr::select(
          bcr_patient_barcode, symbol,
          primary_site, primary_diagnosis_very_simplified) |>
        dplyr::filter(!is.na(primary_site)) |>
        dplyr::distinct() |>
        dplyr::mutate(mut_status = "MUT") |>
        dplyr::mutate(genomic_strata = "gene",
                      variant_type = vartype, consensus_calls = TRUE)

      mutation_status_pancancer <- calls_df |>
        dplyr::filter(mut_status == "HOMDEL") |>
        dplyr::select(bcr_patient_barcode, symbol) |>
        dplyr::mutate(primary_site = "Pancancer",
                      primary_diagnosis_very_simplified = "Pancancer") |>
        dplyr::distinct() |>
        dplyr::mutate(mut_status = "MUT") |>
        dplyr::mutate(genomic_strata = "gene",
                      variant_type = vartype,
                      consensus_calls = TRUE)

      mutation_status_sample <- dplyr::bind_rows(
        mutation_status_sample, mutation_status_pancancer)

      if(vartype == 'cna_ampl'){
        mutation_status_sample <- calls_df |>
          dplyr::filter(mut_status == "AMPL") |>
          dplyr::select(bcr_patient_barcode, symbol,
                        primary_site,
                        primary_diagnosis_very_simplified) |>
          dplyr::filter(!is.na(primary_site)) |>
          dplyr::distinct() |>
          dplyr::mutate(mut_status = "MUT") |>
          dplyr::mutate(genomic_strata = "gene",
                        variant_type = vartype,
                        consensus_calls = TRUE)

        mutation_status_pancancer <- calls_df |>
          dplyr::filter(mut_status == "AMPL") |>
          dplyr::select(bcr_patient_barcode, symbol) |>
          dplyr::mutate(primary_site = "Pancancer",
                        primary_diagnosis_very_simplified = "Pancancer") |>
          dplyr::distinct() |>
          dplyr::mutate(mut_status = "MUT") |>
          dplyr::mutate(genomic_strata = "gene",
                        variant_type = vartype,
                        consensus_calls = TRUE)

        mutation_status_sample <- dplyr::bind_rows(
          mutation_status_sample, mutation_status_pancancer)
      }

      sample_numbers_cohort <- as.data.frame(
        dplyr::select(clinical_df, bcr_patient_barcode,
                      primary_site,
                      primary_diagnosis_very_simplified) |>
          dplyr::distinct() |>
          dplyr::select(-bcr_patient_barcode) |>
          dplyr::filter(!is.na(primary_site)) |>
          dplyr::group_by(primary_site) |>
          dplyr::summarise(tot_samples = dplyr::n(),
                           .groups = "drop")
      )
      pancancer_sample_numbers <- as.data.frame(
        dplyr::select(clinical_df, bcr_patient_barcode) |>
          dplyr::distinct() |>
          dplyr::mutate(primary_site = "Pancancer") |>
          dplyr::select(-bcr_patient_barcode) |>
          dplyr::group_by(primary_site) |>
          dplyr::summarise(tot_samples = dplyr::n(),
                           .groups = "drop")
      )

      sample_numbers_cohort <- dplyr::bind_rows(
        sample_numbers_cohort, pancancer_sample_numbers)
      sample_numbers_cohort$clinical_strata <- 'site'

      mutation_stats <- as.data.frame(
        dplyr::group_by(mutation_status_sample, symbol,
                        primary_site, genomic_strata,
                        variant_type, consensus_calls) |>
          dplyr::summarise(samples_mutated = dplyr::n(), .groups = "drop") |>
          dplyr::left_join(sample_numbers_cohort, by ="primary_site") |>
          dplyr::mutate(percent_mutated = round(
            (samples_mutated / tot_samples) * 100, digits = num_digits))
      )
      mutation_stats <- as.data.frame(
        mutation_stats |>
          dplyr::group_by(primary_site, variant_type) |>
          dplyr::mutate(percentile = dplyr::ntile(
            percent_mutated, 100),
            decile = dplyr::ntile(percent_mutated,10))
      )


      if(clinical_strata == "site_diagnosis"){
        sample_numbers_cohort <- as.data.frame(
          dplyr::select(clinical_df, bcr_patient_barcode,
                        primary_site,
                        primary_diagnosis_very_simplified) |>
            dplyr::distinct() |>
            dplyr::select(-bcr_patient_barcode) |>
            dplyr::filter(!is.na(primary_site)) |>
            dplyr::group_by(primary_site, primary_diagnosis_very_simplified) |>
            dplyr::summarise(tot_samples = dplyr::n(),
                             .groups = "drop")
        )
        pancancer_sample_numbers <- as.data.frame(
          dplyr::select(clinical_df, bcr_patient_barcode) |>
            dplyr::distinct() |>
            dplyr::mutate(primary_site = "Pancancer",
                          primary_diagnosis_very_simplified = "Pancancer") |>
            dplyr::select(-bcr_patient_barcode) |>
            dplyr::group_by(primary_site, primary_diagnosis_very_simplified) |>
            dplyr::summarise(tot_samples = dplyr::n(),
                             .groups = "drop")
        )

        sample_numbers_cohort <- dplyr::bind_rows(
          sample_numbers_cohort, pancancer_sample_numbers)
        sample_numbers_cohort$clinical_strata <- 'site_diagnosis'

        mutation_stats <- as.data.frame(
          dplyr::group_by(
            mutation_status_sample, symbol, primary_site,
            primary_diagnosis_very_simplified,
            genomic_strata, variant_type, consensus_calls) |>
            dplyr::summarise(samples_mutated = dplyr::n(), .groups = "drop") |>
            dplyr::left_join(sample_numbers_cohort,
                             by = c("primary_site",
                                    "primary_diagnosis_very_simplified")) |>
            dplyr::mutate(percent_mutated = round(
              (samples_mutated / tot_samples) * 100, digits = num_digits))
        )
        mutation_stats <- as.data.frame(
          mutation_stats |>
          dplyr::group_by(
            primary_site, primary_diagnosis_very_simplified, variant_type) |>
          dplyr::mutate(percentile = dplyr::ntile(percent_mutated, 100),
                        decile = dplyr::ntile(percent_mutated,10))
        )

      }
    }
  }
  return(mutation_stats)


}

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }


coexpression_mat <- function(
    tumor = "BRCA",
    method = "pearson",
    min_tpm_expressed = 1,
    min_frac_expressed = 0.25,
    min_nonzero_overlap = 0.25,
    var_quantile = 0.25,
    fdr_cutoff = 0.05,
    correlation_threshold = 0.7,
    output_dir = "",
    input_dir = "",
    write_rds = T,
    release = "release45_20251204"){

  rds_fname <- glue::glue(
    "output/{release}/rnaseq/rnaseq_{tumor}.rds")
  corr_df <- NULL
  if(file.exists(rds_fname)){
    expression_data <- readRDS(file=rds_fname)

    ## 0) get TPM values
    mat <- expression_data$matrix$tpm$tumor

    ## transpose -> samples x genes
    t_mat <- t(mat)

    rm(mat, expression_data)
    n_samples <- nrow(t_mat)

    # 1) Filter genes expressed in >= 25% of samples
    #    (SCOPE: gene - biological )
    keep_expr <- colSums(
      t_mat > min_tpm_expressed) >=
      (min_frac_expressed * n_samples)
    t_mat <- t_mat[, keep_expr]

    # 2) Log-transform
    t_mat_log <- log2(t_mat + 1)

    # 3) Remove lowest 25% variance genes
    #    (SCOPE: gene - biological )
    gene_var <- timeSeries::colVars(t_mat_log)
    var_cut  <- quantile(
      gene_var, probs = var_quantile, na.rm = TRUE)

    t_mat_log <- t_mat_log[, gene_var > var_cut]

    # 4) Correlation calculation
    corr_mat <- Hmisc::rcorr(
      t_mat_log, type = method)

    # keep upper triangle only
    corr_mat$r[lower.tri(corr_mat$r, diag = TRUE)] <- NA
    corr_mat$P[lower.tri(corr_mat$P, diag = TRUE)] <- NA

    corr_df <- reshape2::melt(
      corr_mat$r, na.rm = TRUE)
    colnames(corr_df) <-
      c("symbol_A", "symbol_B", "r")


    # 5) Filter gene pairs based on
    #    non-zero overlap (SCOPE: gene pairs - biological)
    min_nz <- n_samples * min_nonzero_overlap

    # Precompute binary expression (samples x genes)
    expr_binary <- t_mat > 0

    nz_mat <- crossprod(expr_binary)
    corr_df$nonzero_overlap <- nz_mat[
          cbind(corr_df$symbol_A, corr_df$symbol_B)]

    # Filter pairs
    corr_df <- corr_df |>
      dplyr::filter(nonzero_overlap >= min_nz)

    pval_df <- reshape2::melt(
      corr_mat$P, na.rm = TRUE)
    colnames(pval_df) <-
      c("symbol_A", "symbol_B", "p_value")

    corr_df <- dplyr::inner_join(
      corr_df, pval_df,
      by = c("symbol_A", "symbol_B"))

    # 6) FDR correction and filtering (SCOPE: gene pairs - statistical)
    corr_df$fdr <- p.adjust(
      corr_df$p_value, method = "BH")
    corr_df <- corr_df |>
      dplyr::filter(
        fdr < fdr_cutoff &
          abs(r) >= correlation_threshold) |>
      dplyr::arrange(dplyr::desc(abs(r)))

    # === 7) Add Median [Q1 - Q3] per gene in log2(TPM + 1) ===
    gene_stats <- data.frame(
      symbol = colnames(t_mat_log),
      median = apply(t_mat_log, 2,
                     median, na.rm = TRUE),
      Q1 = apply(t_mat_log, 2,
                 quantile, probs = 0.25,
                 na.rm = TRUE),
      Q3 = apply(t_mat_log, 2,
                 quantile, probs = 0.75,
                 na.rm = TRUE)
    )

    # Create formatted string "Median [Q1 - Q3]"
    gene_stats$expr_summary <- paste0(
      round(gene_stats$median, 2),
      " [", round(gene_stats$Q1, 2),
      " - ", round(gene_stats$Q3, 2), "]"
    )

    # Join summaries for each gene in the pair
    corr_df <- corr_df |>
      dplyr::left_join(
        gene_stats[, c("symbol", "expr_summary")],
        by = c("symbol_A" = "symbol")
      ) |>
      dplyr::rename(expr_A = expr_summary) |>
      dplyr::left_join(
        gene_stats[, c("symbol", "expr_summary")],
        by = c("symbol_B" = "symbol")
      ) |>
      dplyr::rename(expr_B = expr_summary)

    saveRDS(
      corr_df,
      file=glue::glue("{output_dir}/coexpmat_{tumor}_{release}.rds"))

  }

}
#
# tcga_make_coexpression <- function(tumor = "BRCA", method = "pearson",
#                                    write_rds = T, release = "release27_20201029"){
#   rds_fname <- paste0('data/rnaseq/tcga_rnaseq_',tumor,'_',release,'.rds')
#   corr_df <- NULL
#   if(file.exists(rds_fname)){
#     expression_data <- readRDS(file=rds_fname)
#     ## get FPKM values and perform log transformation (add 1 to avoid log2 of 0)
#     mat <- log2(expression_data$fpkm$matrix + 1)
#
#     ## transpose matrix (rows are tumor samples, columns are genes)
#     t_mat <- t(mat)
#
#     rm(mat)
#     rm(expression_data)
#
#     ## remove genes with zero variance (co-exp correlation undefined)
#     t_mat <- t_mat[, timeSeries::colVars(t_mat) != 0]
#
#     corr_mat <- Hmisc::rcorr(t_mat, type = method)
#
#     ## keep upper triangle only
#     corr_mat$r[lower.tri(corr_mat$r)] <- NA
#     corr_mat$P[lower.tri(corr_mat$P)] <- NA
#
#     corr_df <- as.data.frame(setNames(reshape2::melt(corr_mat$r, na.rm = T), c('symbol_A', 'symbol_B', 'r')))
#     pval_df <- as.data.frame(setNames(reshape2::melt(corr_mat$P, na.rm = T), c('symbol_A', 'symbol_B', 'p_value')))
#
#     corr_df$symbol_A <- as.character(corr_df$symbol_A)
#     corr_df$symbol_B <- as.character(corr_df$symbol_B)
#     pval_df$symbol_A <- as.character(pval_df$symbol_A)
#     pval_df$symbol_B <- as.character(pval_df$symbol_B)
#
#     corr_df <- corr_df |> dplyr::inner_join(pval_df)
#
#     saveRDS(corr_df,file="data/rnaseq/co_expression/corr_mat_",tumor,"_20190326.rds")
#
#   }
#   return(corr_df)
# }

tile_aberration_plot <- function(qgenes, tcga_gene_data, cstrata = "site", vtype = "snv_indel", percentile = FALSE){

  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  title <- 'SNVs/InDels - TCGA'
  color <- 'steelblue'
  if(vtype == 'cna_ampl'){
    title <- 'Copy number amplifications (sCNA) - TCGA'
    color <- 'darkgreen'
  }
  if(vtype == 'cna_homdel'){
    title <- 'Homozygous deletions - (sCNA) - TCGA'
    color <- 'firebrick'
  }

  tcga_gene_stats <- tcga_gene_data |>
    dplyr::filter(clinical_strata == cstrata) |>
    dplyr::filter(primary_site != "Other/Unknown") |>
    dplyr::inner_join(dplyr::select(query_genes_df, symbol))


  gene_candidates_init <- data.frame()
  tcga_ttypes <- sort(unique(tcga_gene_stats$primary_site))
  for(i in 1:length(sort(unique(tcga_gene_stats$symbol)))){
    init <- data.frame('primary_site' <- tcga_ttypes, 'primary_diagnosis_very_simplified' = NA, 'symbol' = sort(unique(tcga_gene_stats$symbol))[i], 'genomic_strata' = 'gene', 'clinical_strata' = cstrata, 'percent_mutated' = 0, 'percentile' = 0, 'variant_type' = vtype, 'consensus_calls' = F, 'fp_driver_gene' = as.logical(NA), decile = 0, stringsAsFactors = F)
    colnames(init) <- c('primary_site','primary_diagnosis_very_simplified','symbol','genomic_strata','clinical_strata','percent_mutated','percentile','variant_type', 'consensus_calls','fp_driver_gene','decile')
    gene_candidates_init <- rbind(gene_candidates_init, init)
    i <- i + 1
  }
  gene_candidates_init <- gene_candidates_init |>
    dplyr::filter(primary_site != 'Other/Unknown' & primary_site != "Pancancer")


  gene_aberrations <- tcga_gene_stats |>
    dplyr::filter(variant_type == vtype) |>
    dplyr::filter(primary_site != "Pancancer" & primary_site != "Other/Unknown")


  site_stats_zero <- tcga_gene_stats |>
    dplyr::select(primary_site,tot_samples) |>
    dplyr::distinct() |>
    dplyr::mutate(samples_mutated = 0)

  pancan_order <- tcga_gene_stats |>
    dplyr::filter(variant_type == vtype & primary_site == "Pancancer") |>
    dplyr::mutate(pancancer_percent_mutated = percent_mutated) |>
    dplyr::mutate(pancancer_percentile = percentile) |>
    dplyr::select(symbol, pancancer_percent_mutated, pancancer_percentile)

  zero_frequency_genes <- dplyr::anti_join(gene_candidates_init, gene_aberrations,
                                           by=c("symbol","primary_site","variant_type")) |>
    #dplyr::select(-c(tot_samples,samples_mutated)) |>
    dplyr::left_join(site_stats_zero)
  gene_aberrations <- dplyr::left_join(dplyr::bind_rows(gene_aberrations, zero_frequency_genes), pancan_order)


  gene_aberrations <- gene_aberrations |>
    dplyr::arrange(pancancer_percent_mutated) |>   # rearrange the df in the order we want
    mutate(symbol = factor(symbol, unique(symbol)))

  p <- ggplot2::ggplot(gene_aberrations,ggplot2::aes(x=primary_site,y=symbol)) +
    ggplot2::geom_text(ggplot2::aes(label = round(percent_mutated,1)), color = "#E69F00") +
    ggplot2::geom_tile(ggplot2::aes(fill = percent_mutated), colour="black",size=0.40)

  if(percentile == T){
    p <- ggplot2::ggplot(gene_aberrations,ggplot2::aes(x=primary_site,y=symbol)) +
      ggplot2::geom_text(ggplot2::aes(label = round(percentile,1))) +
      ggplot2::geom_tile(ggplot2::aes(fill = percentile), colour="black",size=0.40)
  }
  p <- p +
    #add border white colour of line thickness 0.25
    #ggplot2::geom_tile(aes(fill = pct_mutated), colour="black",size=0.40)+
    ggplot2::scale_fill_gradient(low = "white", high = color) +
    #remove x and y axis labels
    ggplot2::labs(x="",y="")+
    ggplot2::ggtitle(title) +
    #remove extra space
    #scale_y_discrete(expand=c(0,0))+
    #ggplot2::coord_fixed(ratio = 1.3)+
    #set a base size for all fonts
    ggplot2::theme_grey(base_size=14)+
    #theme options
    ggplot2::theme(
      #bold font for both axis text
      #legend.title = ggplot2::element_text()
      axis.text = ggplot2::element_text(face="bold",family = "Helvetica", size = 14),
      #set thickness of axis ticks
      axis.ticks = ggplot2::element_line(size=0.2),
      #remove plot background
      plot.background = ggplot2::element_blank(),
      #remove plot border
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}

#
#
# tcga_significant_coexpression <- function(tumor = "BRCA", basedir = "/Volumes/sigven/tcga/data/rnaseq/co_expression"){
#   rds_fname <- paste0(basedir,"/corr_mat_",tumor,"_20190917.rds")
#   if(!(file.exists(rds_fname))){
#     return(0)
#   }
#   ## get moderate ([0.4,0.6]), strong ([0.6,0.8]), and very strong (>0.8) positively correlated genes
#   ## get moderate ([-0.4,-0.6]), strong ([-0.6,-0.8]), and very strong (< -0.8) negatively correlated genes
#   significant_co_expression <- readRDS(rds_fname) |>
#     dplyr::filter(p_value < 0.000001 & ((r > 0 & r >= 0.4) | (r < 0 & r <= -0.4))) |>
#     dplyr::mutate(tumor = tumor)
#
#   return(significant_co_expression)
# }
#
# tcga_geneset_expression <- function(qgenes, tumor = "BRCA",
#                                      rnaseq_datatype = "fpkm",
#                                      tcga_clinical = NULL,
#                                      tcga_release = "release27_20201029"){
#   query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
#   rnaseq_rds_file <- paste0("/Users/sigven/research/tcga/data/rnaseq/tcga_rnaseq_",tumor,"_",tcga_release,".rds")
#   if(!file.exists(rnaseq_rds_file)){
#     next
#   }
#   m <- readRDS(file=rnaseq_rds_file)
#   rnaseq_df <- dfrtopics::gather_matrix(m$fpkm$matrix, col_names = c('symbol','tumor_sample_barcode','expr'))
#   if(rnaseq_datatype == "counts"){
#     rnaseq_df <- dfrtopics::gather_matrix(m$counts$matrix, col_names = c('symbol','tumor_sample_barcode','expr'))
#   }
#   rnaseq_df <- rnaseq_df |>
#     dplyr::inner_join(query_genes_df,by=c("symbol")) |>
#     dplyr::mutate(measure = rnaseq_datatype) |>
#     dplyr::filter(stringr::str_detect(tumor_sample_barcode,"-0[0-9][A-Z]$")) |>
#     dplyr::mutate(sample_type = dplyr::if_else(!stringr::str_detect(tumor_sample_barcode,"-01[A-Z]$"),as.character(NA),"Primary Tumor")) |>
#     dplyr::mutate(sample_type = dplyr::if_else(stringr::str_detect(tumor_sample_barcode,"-06[A-Z]$"),"Metastic",as.character(sample_type))) |>
#     dplyr::mutate(bcr_patient_barcode = stringr::str_replace(tumor_sample_barcode,"-[0-9][0-9][A-Z]$","")) |>
#     dplyr::mutate(log2_expr = log2(expr + 1))
#
#   if(!is.null(tcga_clinical)){
#     rnaseq_df <- rnaseq_df |>
#       dplyr::left_join(dplyr::select(tcga_clinical, bcr_patient_barcode, primary_site,
#                                      primary_diagnosis_simplified, tumor), by=c("bcr_patient_barcode"))
#   }else{
#     rnaseq_df <- rnaseq_df |> dply::mutate(tumor = tumor)
#   }
#
#   rnaseq_df$median_expr <- median(rnaseq_df$expr)
#   rnaseq_df$median_log2_expr <- median(rnaseq_df$log2_expr)
#   rnaseq_df$lower_log2_expr <- quantile(rnaseq_df$log2_expr)[2]
#   rnaseq_df$upper_log2_expr <- quantile(rnaseq_df$log2_expr)[4]
#
#   return(rnaseq_df)
# }
#

# get_sample_vcf_data <- function(tcga_calls, sample_barcode = "TCGA-A2-A04W-01A", min_callers = 1){
#
#   sample_calls <- tcga_calls |>
#     dplyr::filter(Tumor_Sample_Barcode == sample_barcode) |>
#     dplyr::select(CHROM,POS,REF,ALT,t_depth,t_alt_count,n_depth,tumor,
#                   primary_site, algorithms, Tumor_Sample_Barcode) |>
#     dplyr::mutate(primary_site = stringr::str_replace_all(primary_site,"/| ","_")) |>
#     dplyr::mutate(num_callers = stringr::str_count(algorithms,pattern=",") + 1) |>
#     dplyr::filter(num_callers >= min_callers) |>
#     dplyr::mutate(FILTER = "PASS", QUAL = ".", ID = ".") |>
#     dplyr::mutate(TVAF = round(as.numeric(t_alt_count)/t_depth, digits = 4)) |>
#     dplyr::mutate(INFO = paste0("TDP=", t_depth,";TVAF=", TVAF,";CDP=", n_depth,";ALGS=",
#                                 algorithms,";COHORT=", tumor,";PRIMARY_SITE=", primary_site,
#                                 ";TUMOR_SAMPLE_BARCODE=",Tumor_Sample_Barcode)) |>
#     dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) |>
#     dplyr::distinct()
#
#   header_lines <- c("##fileformat=VCFv4.2",
#                     "##INFO=<ID=TUMOR_SAMPLE_BARCODE,Number=.,Type=String,Description=\"TCGA tumor smple barcode\">",
#                     "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Sequencing depth at variant site - tumor\">",
#                     "##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Sequencing depth at variant site - control\">",
#                     "##INFO=<ID=TVAF,Number=1,Type=Float,Description=\"Allelic fraction of alternate allele - tumor\">",
#                     "##INFO=<ID=ALGS,Number=.,Type=String,Description=\"List of variant callers/algorithms that detected the variant\">",
#                     "##INFO=<ID=COHORT,Number=1,Type=String,Description=\"TCGA tumor sequencing cohort\">",
#                     "##INFO=<ID=PRIMARY_SITE,Number=1,Type=String,Description=\"Primary tumor site\">",
#                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
#
#   return(list('vcf_content' = sample_calls, 'vcf_header_lines' = header_lines))
#
# }

make_pcgr_command <- function(pcgr_data_dir = NULL,
                              output_dir = NULL,
                              conf = NULL,
                              sample_id = NULL,
                              sample_id_tcga = NULL,
                              cna_tsv = NULL,
                              query_vcf = NULL,
                              tumor_only = FALSE,
                              tumor_dp_tag = "_NA_",
                              tumor_af_tag = "_NA",
                              call_conf_tag = "_NA_",
                              show_noncoding = FALSE,
                              assay = "WGS",
                              cpsr_json = FALSE,
                              preserved_info_tags = NULL,
                              include_artefact_signatures = FALSE,
                              prevalence_ref_signatures = 5,
                              vcf2maf = TRUE,
                              vep_regulatory = FALSE,
                              target_size_mb = 34,
                              exclude_pon = FALSE,
                              exclude_likely_hom_germline = FALSE,
                              exclude_likely_het_germline = FALSE,
                              exclude_dbsnp_nonsomatic = TRUE,
                              exclude_nonexonic = TRUE,
                              p_site = 9,
                              nonfloating_toc = FALSE,
                              tumor_purity = NULL,
                              tumor_ploidy = NULL,
                              assembly = 'grch38',
                              tumor_dp_min = 0,
                              tumor_af_min = 0){

  pcgr_command <- paste0("pcgr ",
                         " --tumor_site ", p_site,
                         " --input_vcf ",query_vcf,".gz",
                         " --assay ",assay,
                         " --genome_assembly ",assembly,
                         " --include_trials",
                         " --prevalence_reference_signatures ",
                         prevalence_ref_signatures,
                         " --tumor_dp_tag ",tumor_dp_tag,
                         " --tumor_af_tag ",tumor_af_tag,
                         " --call_conf_tag ",call_conf_tag,
                         " --sample_id ",sample_id,
                         " --output_dir ",output_dir,
                         " --pcgr_dir ",pcgr_data_dir,
                         " --tumor_dp_min ", tumor_dp_min,
                         " --estimate_tmb",
                         " --tumor_af_min ",tumor_af_min,
                         " --target_size_mb ", target_size_mb,
                         " --force_overwrite --no_docker")
  if(assay != "TARGETED"){
    pcgr_command <- paste0(pcgr_command, " --estimate_msi --estimate_signatures")
    if(include_artefact_signatures == T){
      pcgr_command <- paste0(pcgr_command," --include_artefact_signatures")
    }
  }

  if(cpsr_json == T){
    cpsr_gzipped_fname <- file.path(output_dir, paste0(sample_id_tcga,".cpsr.",assembly,".json.gz"))
    pcgr_command <- paste0(pcgr_command," --cpsr_report ",
                           cpsr_gzipped_fname)
  }

  if(!is.null(preserved_info_tags)){
    pcgr_command <- paste0(pcgr_command," --preserved_info_tags ",preserved_info_tags)

  }
  if(!is.null(tumor_purity)){
    pcgr_command <- paste0(pcgr_command," --tumor_purity ",tumor_purity)
  }
  if(!is.null(tumor_purity)){
    pcgr_command <- paste0(pcgr_command," --tumor_ploidy ",tumor_ploidy)
  }
  if(!is.null(cna_tsv)){
    pcgr_command <- paste0(pcgr_command, " --logr_gain 0.4 --logr_homdel -0.4 --input_cna ",cna_tsv)
  }
  if(vcf2maf == T){
    pcgr_command <- paste0(pcgr_command," --vcf2maf")
  }

  if(vep_regulatory == T){
    pcgr_command <- paste0(pcgr_command," --vep_regulatory")
  }

  show_coding_type <- "NOCODING"
  if(show_noncoding == T){
    pcgr_command <- paste0(pcgr_command," --show_noncoding")
    show_coding_type <- "ALL"
  }

  if(nonfloating_toc == T){
    pcgr_command <- paste0(pcgr_command," --report_nonfloating_toc")
  }

  seqtype = "TC"
  if(tumor_only){
    pcgr_command <- paste0(pcgr_command, " --tumor_only")
    seqtype = "TO"
    if(exclude_likely_het_germline == T){
      pcgr_command <- paste0(pcgr_command, " --exclude_likely_het_germline")

    }
    if(exclude_likely_hom_germline == T){
      pcgr_command <- paste0(pcgr_command, " --exclude_likely_hom_germline")

    }
    if(exclude_dbsnp_nonsomatic == T){
      pcgr_command <- paste0(pcgr_command, " --exclude_dbsnp_nonsomatic")

    }
    if(exclude_nonexonic == T){
      pcgr_command <- paste0(pcgr_command, " --exclude_nonexonic")

    }

  }
  log_file <- paste0(output_dir,"/",sample_id,".", assembly,".", assay,".",seqtype,".",show_coding_type,".pcgr.log")

  if(tumor_dp_min > 0 | tumor_af_min > 0){
    log_file <- paste0(output_dir,"/",sample_id,".", assembly,".", assay,".",seqtype,".",show_coding_type,".DP_AF.pcgr.log")
  }

  pcgr_command <- paste0(pcgr_command, " > ", log_file)

  return(pcgr_command)

}

make_cpsr_command <-
  function(pcgr_data_dir = NULL,
           output_dir = NULL,
           conf = NULL,
           sample_id = NULL,
           query_vcf = NULL,
           preserved_info_tags = NULL,
           classify_all = TRUE,
           nonfloating_toc = TRUE,
           diagnostic_grade_only = FALSE,
           clinvar_ignore_noncancer = FALSE,
           panel_id = "-1",
           pcgr_integration = F,
           assembly = 'grch38',
           pop_gnomad = 'nfe',
           table_display = "light",
           ignore_noncoding = F){

  sample_id2 <- paste0(sample_id, "_GL_PID_", stringr::str_replace_all(panel_id,",","_"))
  sample_id2 <- paste0(sample_id2,"_",toupper(pop_gnomad))
  cpsr_command <- paste0("cpsr ",
                         " --input_vcf ",query_vcf,".gz",
                         " --genome_assembly ",assembly,
                         " --secondary_findings --gwas_findings --panel_id ",panel_id,
                         " --report_table_display ",table_display,
                         " --pop_gnomad ",pop_gnomad,
                         " --output_dir ",output_dir,
                         " --pcgr_dir ", pcgr_data_dir)
  if(!is.null(preserved_info_tags)){
    cpsr_command <- paste0(cpsr_command," --preserved_info_tags ",preserved_info_tags)
  }

  #sample_id2 <- sample_id
  if(classify_all == T){
    cpsr_command <- paste0(cpsr_command," --classify_all")
    sample_id2 = paste0(sample_id2, "_CLASSIFY_ALL")
  }
  if(nonfloating_toc == T){
    cpsr_command <- paste0(cpsr_command," --report_nonfloating_toc")

  }
  if(clinvar_ignore_noncancer == T){
    cpsr_command <- paste0(cpsr_command," --clinvar_ignore_noncancer")
    sample_id22 = paste0(sample_id2, "_IGNORE_NONCANCER")

  }
  if(diagnostic_grade_only == T & panel_id != "0"){
    sample_id2 = paste0(sample_id2, "_DGRADE")
    cpsr_command <- paste0(cpsr_command," --diagnostic_grade_only")
  }

  if(pcgr_integration == T){
    sample_id2 = sample_id
  }

  cpsr_command <- paste0(cpsr_command, " --sample_id ", sample_id2)
  cpsr_command <- paste0(cpsr_command, " --force_overwrite --no_docker > ",output_dir,"/",sample_id,".", assembly,".cpsr.log")

  return(cpsr_command)

}

#
# write_vcf_df <- function(vcf_df, vcf_path, vcf_prefix, header_lines){
#
#   sample_vcf_fname <- paste0(vcf_path,'/',vcf_prefix,".vcf")
#   sample_vcf_content_fname <- paste0(vcf_path,"/",vcf_prefix,".vcf_content.tsv")
#   write(header_lines,file=sample_vcf_fname,sep="\n")
#
#   sample_vcf <- vcf_df[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]
#   write.table(sample_vcf, file=sample_vcf_content_fname,sep="\t",col.names = F,quote=F, row.names = F)
#
#   system(paste0("cat ",sample_vcf_content_fname," | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",sample_vcf_fname))
#   system(paste0("cat ",sample_vcf_content_fname," | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",sample_vcf_fname))
#   system(paste0("bgzip -c ",sample_vcf_fname," > ",sample_vcf_fname,".gz"))
#   system(paste0("tabix -p vcf ",sample_vcf_fname,".gz"))
#
#   system(paste0('rm -f ',sample_vcf_content_fname))
#   #system(paste0('rm -f ',sample_vcf_fname))
# }
#
# write_cna_bed <- function(cna_df, bed_path, bed_prefix){
#
#   sample_bed_fname <- paste0(bed_path,'/',bed_prefix,".bed")
#   sample_bed_content_fname <- paste0(bed_path,"/",bed_prefix,".vcf_content.tsv")
#   #write(header_lines,file=sample_vcf_fname,sep="\n")
#
#   assertable::assert_colnames(cna_df, c("Chromosome","Start","End","Segment_Mean"),
#                               only_colnames = F, quiet = T)
#   #write.table(sample_vcf, file=sample_vcf_content_fname,sep="\t",col.names = F,quote=F, row.names = F)

  # system(paste0("cat ",sample_vcf_content_fname," | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",sample_vcf_fname))
  # system(paste0("cat ",sample_vcf_content_fname," | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 | awk '{printf(\"%s\\n\",$0);}' >> ",sample_vcf_fname))
  # system(paste0("bgzip -c ",sample_vcf_fname," > ",sample_vcf_fname,".gz"))
  # system(paste0("tabix -p vcf ",sample_vcf_fname,".gz"))
  #
  # system(paste0('rm -f ',sample_vcf_content_fname))
  #system(paste0('rm -f ',sample_vcf_fname))
#}


# get_unambiguous_gene_aliases <- function(
#   basedir = "/Users/sigven/project_data/analysis__tcga/tcga/data-raw/"){
#
#   # gene_info_fname <- file.path(basedir,"ncbi_gene","gene_info.rds")
#
#   gene_info_fname <- file.path(basedir, "ncbi_gene",
#                                "Homo_sapiens.gene_info.gz")
#
#   primary_to_primary_all <-
#     readr::read_tsv(gene_info_fname,
#                     show_col_types = F) |>
#     janitor::clean_names() |>
#     dplyr::rename(entrezgene = gene_id,
#                   alias = symbol) |>
#     dplyr::select(entrezgene, alias) |>
#     dplyr::distinct()
#
#   # primary_to_primary_all <- readRDS(file=gene_info_fname) |>
#   #   dplyr::select(symbol, entrezgene) |>
#   #   dplyr::mutate(alias = symbol) |>
#   #   dplyr::distinct()
#
#   gene_synonyms <-
#     readr::read_tsv(gene_info_fname,
#                     show_col_types = F) |>
#     janitor::clean_names() |>
#     # readRDS(
#     # file = gene_info_fname) |>
#     dplyr::select(symbol, synonyms, gene_id) |>
#     dplyr::rename(alias = synonyms, entrezgene = gene_id) |>
#     tidyr::separate_rows(alias, sep="\\|") |>
#     dplyr::filter(!(symbol == "H3P10" & alias == "p16")) |>
#     dplyr::filter(nchar(alias) > 2) |>
#     dplyr::filter(stringr::str_detect(alias,"-|[0-9]")) |>
#     dplyr::anti_join(primary_to_primary_all, by = "alias")
#
#   unique_aliases <- as.data.frame(
#     gene_synonyms |>
#       dplyr::group_by(alias) |>
#       dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
#       dplyr::filter(n == 1) |>
#       dplyr::select(-n)
#   )
#
#   ambig_aliases <- as.data.frame(
#     gene_synonyms |>
#       dplyr::group_by(alias) |>
#       dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
#       dplyr::filter(n > 1) |>
#       dplyr::inner_join(gene_synonyms, by = "alias")
#   )
#
#   gene_dictionary <- as.data.frame(gene_synonyms |>
#     dplyr::inner_join(unique_aliases, by = "alias") |>
#     dplyr::bind_rows(primary_to_primary_all) |>
#     dplyr::filter(!(symbol == "HRAS" & (alias == "KRAS" | alias == "c-K-ras" | alias == "c-Ki-ras"))) |>
#     dplyr::arrange(symbol, alias) |>
#     dplyr::mutate(symbol = dplyr::if_else(
#       is.na(symbol) & !is.na(alias),
#       as.character(alias),
#       as.character(symbol)
#     ))
#   )
#
#   return(gene_dictionary)
# }


# parse_transvar_file <- function(transvar_output_fname, build = "hg19") {
#
#   transvar_output_raw <- as.data.frame(
#     readr::read_tsv(transvar_output_fname,
#                     col_names = T,
#                     show_col_types = F) |>
#       janitor::clean_names() |>
#       tidyr::separate(coordinates_g_dna_c_dna_protein,
#                       c("gdna", "cdna", "hgvsp"), sep = "/") |>
#       dplyr::rename(transcript_id = transcript, transvar_id = input) |>
#       dplyr::mutate(transcript_id = stringr::str_replace(
#         transcript_id, " \\(protein_coding\\)", "")) |>
#       dplyr::filter(stringr::str_detect(gdna, ">|del[A-Z]{1,}|ins[A-Z]{1,}")) |>
#       dplyr::distinct()
#   )
#
#   transvar_gdna_maps <- list()
#
#   ## Main gDNA mapping - SNVS
#   transvar_gdna_maps[['snv']] <- transvar_output_raw |>
#     dplyr::filter(stringr::str_detect(
#       gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}(A|G|C|T)>(A|G|C|T)"
#     )) |>
#     dplyr::select(gdna, transvar_id, transcript_id) |>
#     dplyr::distinct()
#
#   ## Alternative gDNA mappings - SNVs
#   transvar_gdna_maps[['snv_alt']] <- transvar_output_raw |>
#     dplyr::filter(stringr::str_detect(
#       gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}(A|G|C|T)>(A|G|C|T)"
#     )) |>
#     dplyr::select(info, transcript_id, transvar_id) |>
#     dplyr::mutate(
#       info =
#         stringr::str_replace_all(
#           info,
#           "CSQN=\\S+;reference_codon=\\S+;candidate_codons=((A|C|G|T){3},?){1,};",
#           "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, "source=GENCODE$", "NA")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, "candidate_mnv_variants=|candidate_snv_variants=", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";?aliases=ENSP[0-9]{1,}\\.[0-9]{1,}", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";source=GENCODE$", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";", ",")) |>
#     dplyr::filter(!is.na(info) & info != "NA" & nchar(info) > 0) |>
#     dplyr::rename(gdna = info) |>
#     tidyr::separate_rows(gdna, sep = ",") |>
#     dplyr::filter(gdna != "NA" & !is.na(gdna)) |>
#     dplyr::select(gdna, transvar_id, transcript_id) |>
#     dplyr::filter(nchar(gdna) > 0) |>
#     dplyr::distinct()
#
#
#
#   ## Main gDNA mappings - insertions
#   transvar_gdna_map_ins <-
#     dplyr::select(transvar_output_raw, transvar_id,
#                   transcript_id, gdna) |>
#     dplyr::filter(stringr::str_detect(
#       gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}_[0-9]{1,}ins(A|G|C|T){1,}$"
#     )) |>
#     dplyr::mutate(altbase = stringr::str_match(gdna,"[A-Z]{1,}$")[,1],
#                   chr_start = stringr::str_replace_all(
#                     stringr::str_match(gdna,"g\\.[0-9]{1,}_"),
#                     "g\\.|_","")) |>
#     dplyr::mutate(chromosome =
#                     stringr::str_split_fixed(gdna, ":", 2)[, 1]
#     ) |>
#     dplyr::mutate(chr_stop = as.integer(chr_start) + nchar(altbase) - 1)
#
#
#   seqinfo <- GenomeInfoDb::Seqinfo(
#     seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
#     seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
#     genome = 'hg19')
#
#   if(build == 'hg38'){
#     seqinfo <- GenomeInfoDb::Seqinfo(
#       seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)),
#       seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)),
#       genome = 'hg38')
#   }
#
#   if(NROW(transvar_gdna_map_ins) > 0){
#
#     vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(
#       transvar_gdna_map_ins,
#       keep.extra.columns = T,
#       seqinfo = seqinfo,
#       seqnames.field = 'chromosome',
#       start.field = 'chr_start',
#       end.field = 'chr_start',
#       ignore.strand = T,
#       starts.in.df.are.0based = F)
#
#     genome_seq <- BSgenome.Hsapiens.UCSC.hg19
#     if(build == "hg38"){
#       genome_seq <- BSgenome.Hsapiens.UCSC.hg38
#     }
#     refbase_seq <- Biostrings::getSeq(genome_seq, vcf_df_gr)
#     df_ref <- data.frame('refbase'=toupper(unlist(strsplit(toString(refbase_seq),", "))),
#                          stringsAsFactors=F)
#     transvar_gdna_map_ins <- transvar_gdna_map_ins |>
#       dplyr::bind_cols(df_ref) |>
#       dplyr::mutate(altbase = paste0(refbase,altbase)) |>
#       dplyr::mutate(chromosome = stringr::str_replace(chromosome,"chr",""))
#   }
#
#
#   ## Main gDNA mappings - block substitutions (delins)
#   transvar_gdna_mnv_all <-
#     dplyr::select(transvar_output_raw, transvar_id,
#                   transcript_id, gdna, info) |>
#     dplyr::filter(stringr::str_detect(
#       gdna, "chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}_[0-9]{1,}del(A|G|C|T){1,}ins(A|G|C|T){1,}$"
#     )) |>
#     dplyr::mutate(
#       info =
#         stringr::str_replace_all(
#           info,
#           "CSQN=\\S+;reference_codon=\\S+;candidate_codons=((A|C|G|T){3},?){1,};",
#           "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, "source=GENCODE$", "NA")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, "candidate_mnv_variants=|candidate_snv_variants=", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";?aliases=ENSP[0-9]{1,}\\.[0-9]{1,}", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";source=GENCODE$", "")) |>
#     dplyr::mutate(info = stringr::str_replace_all(
#       info, ";", ",")) |>
#     dplyr::filter(!is.na(info) & info != "NA" & nchar(info) > 0) |>
#     dplyr::rename(gdna_alts = info) |>
#     tidyr::separate_rows(gdna_alts, sep = ",") |>
#     dplyr::filter(gdna_alts != "NA" & !is.na(gdna_alts))
#
#   transvar_gdna_maps[['mnv']] <- transvar_gdna_mnv_all |>
#     dplyr::select(transvar_id, transcript_id, gdna) |>
#     dplyr::distinct()
#
#   ## Alternative gDNA mappings - block substitutions (delins)
#   transvar_gdna_maps[['mnv_alt']] <- transvar_gdna_mnv_all |>
#     dplyr::select(transvar_id, transcript_id, gdna_alts) |>
#     dplyr::filter(nchar(gdna_alts) > 0) |>
#     dplyr::rename(gdna = gdna_alts) |>
#     dplyr::distinct()
#
#
#   transvar_gdna_mapped_all <- do.call(rbind, transvar_gdna_maps) |>
#     magrittr::set_rownames(NULL) |>
#     dplyr::mutate(alleles = stringr::str_match(gdna,"[A-Z]>[A-Z]|(del[A-Z]{1,}ins[A-Z]{1,})")[,1],
#                   chr_start = stringr::str_replace(
#                     stringr::str_match(gdna,"g\\.[0-9]{1,}"),
#                     "g\\.","")) |>
#     dplyr::mutate(chromosome = stringr::str_replace(
#       stringr::str_split_fixed(gdna, ":", 2)[, 1], "chr", ""
#     )) |>
#     dplyr::mutate(refbase = stringr::str_replace(
#       stringr::str_split_fixed(alleles,">|ins",2)[,1],"del",""
#     )) |>
#     dplyr::mutate(altbase = stringr::str_split_fixed(
#       alleles,">|ins",2)[,2]) |>
#     dplyr::mutate(chr_stop =
#                     as.integer(chr_start) + nchar(altbase) - 1
#     ) |>
#     dplyr::select(-c(alleles)) |>
#     dplyr::bind_rows(transvar_gdna_map_ins) |>
#     dplyr::arrange(chromosome, chr_start) |>
#     dplyr::mutate(transcript_no_version = stringr::str_replace_all(
#       transcript_id,"\\.[0-9]{1,}$",""
#     ))
#
#   return(transvar_gdna_mapped_all)
#
# }

# get_tcga_driver_mutations <- function(
#     data_raw_dir = NA){
#
#   tcga_pcdm_raw_fname <-
#     file.path(data_raw_dir,
#               "driver_mutations",
#               "Mutation.CTAT.3D.Scores.txt.gz")
#
#   tcga_pcdm_raw <- as.data.frame(
#     readr::read_tsv(
#       file = tcga_pcdm_raw_fname,
#       show_col_types = F) |>
#       janitor::clean_names() |>
#       dplyr::select(gene,transcript,protein_change,
#                     recurrence, code, new_linear_functional_flag,
#                     new_linear_cancer_focused_flag,
#                     new_3d_mutational_hotspot_flag) |>
#       dplyr::filter(new_linear_functional_flag == 1 |
#                       new_linear_cancer_focused_flag == 1 |
#                       new_3d_mutational_hotspot_flag == 1) |>
#       dplyr::distinct() |>
#       dplyr::mutate(DISCOVERY_APPROACH_1 = dplyr::if_else(
#         new_linear_functional_flag == 1,'CTAT_POPULATION','','')) |>
#       dplyr::mutate(DISCOVERY_APPROACH_2 = dplyr::if_else(
#         new_linear_cancer_focused_flag == 1,'CTAT_CANCER','','')) |>
#       dplyr::mutate(DISCOVERY_APPROACH_3 = dplyr::if_else(
#         new_3d_mutational_hotspot_flag == 1,'3D_STRUCTURE','','')) |>
#       dplyr::mutate(driver_mutation = paste(
#         DISCOVERY_APPROACH_1, DISCOVERY_APPROACH_2,
#         DISCOVERY_APPROACH_3,sep="|")) |>
#       dplyr::mutate(driver_mutation = stringr::str_replace(
#         driver_mutation,"(^\\|)|(\\|$)",""))  |>
#       dplyr::rename(SYMBOL = gene,
#                     Feature = transcript,
#                     HGVSp_Short = protein_change,
#                     tumor = code) |>
#       #dplyr::mutate(transvar_id = paste(symbol,hgvsp_short,sep=":")) |>
#       dplyr::mutate(
#         num_approaches = new_linear_functional_flag +
#           new_linear_cancer_focused_flag +
#           new_3d_mutational_hotspot_flag) |>
#       dplyr::select(-c(DISCOVERY_APPROACH_1,DISCOVERY_APPROACH_2,DISCOVERY_APPROACH_3,
#                        new_linear_functional_flag,new_linear_cancer_focused_flag,
#                        new_3d_mutational_hotspot_flag)) |>
#       dplyr::filter(num_approaches >= 2)
#   )
#
#   driver_mutations <- list()
#   driver_mutations[['raw']] <- tcga_pcdm_raw
#   driver_mutations[['by_transcript_id']] <-
#     tcga_pcdm_raw |>
#     dplyr::select(-c(SYMBOL, recurrence,
#                      tumor, num_approaches)) |>
#     dplyr::distinct()
#
#   driver_mutations[['by_symbol']] <-
#     tcga_pcdm_raw |>
#     dplyr::select(-c(Feature, recurrence,
#                      tumor, num_approaches)) |>
#     dplyr::distinct()
#
#   return(driver_mutations)
# }

# get_curated_docm_mutations <- function(
#     data_raw_dir = NA){
#
#   docm_curated_mutations <-
#     as.data.frame(
#       readr::read_tsv(
#         file = file.path(
#           data_raw_dir,
#           "curated_mutations",
#           "docm.grch38.tsv.gz"),
#         show_col_types = F) |>
#         dplyr::select(
#           disease, pmid, docm_id,
#           chrom, pos, ref, alt) |>
#         dplyr::mutate(docm_id = stringr::str_replace_all(
#           docm_id, "\\*","X"
#         )) |>
#         dplyr::mutate(docm_id = stringr::str_replace_all(
#           docm_id, ">","_"
#         )) |>
#         dplyr::mutate(
#           ref = as.character(ref),
#           alt = as.character(alt)) |>
#         dplyr::mutate(
#           disease = stringr::str_replace_all(
#             disease," ","_")) |>
#         dplyr::rename(
#           CHROM = chrom,
#           POS = pos,
#           REF = ref,
#           ALT = alt) |>
#         dplyr::group_by(CHROM, POS,
#                         REF, ALT) |>
#         dplyr::summarise(
#           docm_disease = paste(unique(disease),
#                                collapse = ","),
#           docm_pmid = paste(unique(sort(pmid)),
#                             collapse = ","),
#           docm_id = paste(unique(docm_id),
#                           collapse = ","),
#           .groups = "drop")
#     )
#
#   return(docm_curated_mutations)
# }
#

# load_cancer_hotspots <- function(
#     data_raw_dir = NA,
#     gOncoX = NULL){
#
#   cancer_hotspots <- list()
#
#   cancer_hotspots[['indel']] <-
#     openxlsx::read.xlsx(
#       xlsxFile = file.path(
#         data_raw_dir,
#         "cancer_hotspots",
#         "hotspots_v2.xlsx"),
#       sheet = 2,
#       startRow = 1) |>
#     janitor::clean_names() |>
#     dplyr::select(
#       hugo_symbol, qvalue,
#       variant_amino_acid, samples, tm) |>
#     dplyr::mutate(qvalue = stringr::str_trim(
#       format(as.numeric(
#         as.character(qvalue)),
#         scientific = TRUE, digits = 2))) |>
#     dplyr::mutate(var_aa = stringr::str_replace(
#       variant_amino_acid,"\\*","X")) |>
#     dplyr::mutate(hgvsp = paste0(
#       "p.",stringr::str_replace(
#         var_aa, ":[0-9]{1,}$",""))) |>
#     dplyr::mutate(amino_acid_position = stringr::str_split_fixed(
#       tm, " ",2)[,2]) |>
#     dplyr::select(-c(var_aa,
#                      variant_amino_acid, tm)) |>
#     tidyr::separate_rows(samples, sep="\\|") |>
#     dplyr::rename(tumor_type_freq = samples) |>
#     rename_hotspot_tumor_types() |>
#     tidyr::separate(tumor_type_freq,
#                     into = c("ttype","freq"),
#                     sep = ":",
#                     remove = T) |>
#     dplyr::mutate(freq = as.integer(freq)) |>
#     dplyr::mutate(reference_amino_acid = as.character(NA),
#                   codon = as.character(NA)) |>
#     dplyr::select(
#       hugo_symbol, qvalue, ttype, freq, hgvsp, amino_acid_position,
#       reference_amino_acid, ttype) |>
#     dplyr::distinct()
#
#   cancer_hotspots[['indel']] <- resolve_gene_symbol(
#     df = cancer_hotspots[['indel']], gOncoX = gOncoX
#   ) |>
#     dplyr::mutate(MUTATION_HOTSPOT = paste0(
#       symbol, "|", entrezgene, "|",
#       amino_acid_position, "|",
#       qvalue)) |>
#     dplyr::mutate(MUTATION_HOTSPOT2 = stringr::str_replace_all(
#       MUTATION_HOTSPOT, "\\*","X"
#     )) |>
#     dplyr::mutate(hgvsp2 = stringr::str_replace_all(
#       hgvsp, "\\*","X"
#     ))
#
#   cancer_hotspots[['snv']] <-  openxlsx::read.xlsx(
#     xlsxFile = file.path(
#       data_raw_dir,
#       "cancer_hotspots",
#       "hotspots_v2.xlsx"),
#     sheet = 1,
#     startRow = 1) |>
#     janitor::clean_names() |>
#     dplyr::select(
#       hugo_symbol, mutation_count, amino_acid_position,
#       reference_amino_acid, qvalue, qvalue_pancan,
#       qvaluect, detailed_cancer_types, variant_amino_acid,
#       samples, total_samples) |>
#     dplyr::mutate(reference_amino_acid = stringr::str_replace(
#       reference_amino_acid,":[0-9]{1,}","")) |>
#     dplyr::mutate(variant_amino_acid = stringr::str_replace(
#       variant_amino_acid,":[0-9]{1,}","")) |>
#     dplyr::mutate(qvalue = stringr::str_trim(
#       format(as.numeric(
#         as.character(qvalue)),
#         scientific = TRUE, digits = 2))) |>
#     dplyr::mutate(hgvsp = paste0(
#       "p.", reference_amino_acid,
#       amino_acid_position, variant_amino_acid)) |>
#     dplyr::mutate(codon = paste0(
#       "p.", reference_amino_acid,
#       amino_acid_position)) |>
#     tidyr::separate_rows(samples, sep="\\|") |>
#     dplyr::rename(tumor_type_freq = samples) |>
#     rename_hotspot_tumor_types() |>
#     tidyr::separate(tumor_type_freq,
#                     into = c("ttype","freq"),
#                     sep = ":",
#                     remove = T) |>
#     dplyr::mutate(freq = as.integer(freq)) |>
#     dplyr::select(
#       hugo_symbol, qvalue, ttype, freq, hgvsp, codon,
#       amino_acid_position,
#       reference_amino_acid, variant_amino_acid, ttype) |>
#     dplyr::distinct()
#
#
#   cancer_hotspots[['snv']] <- resolve_gene_symbol(
#     df = cancer_hotspots[['snv']], gOncoX = gOncoX) |>
#     dplyr::mutate(MUTATION_HOTSPOT = paste0(
#       symbol, "|", entrezgene, "|",
#       reference_amino_acid,
#       amino_acid_position, "|",
#       variant_amino_acid,"|",
#       as.character(qvalue))) |>
#     dplyr::mutate(MUTATION_HOTSPOT2 = stringr::str_replace_all(
#       MUTATION_HOTSPOT, "\\*","X"
#     )) |>
#     dplyr::mutate(hgvsp2 = stringr::str_replace_all(
#       hgvsp, "\\*","X"
#     ))
#
#
#   site_freqs <- list()
#   site_freqs[['snv']] <- as.data.frame(
#     cancer_hotspots[['snv']] |>
#       dplyr::select(
#         symbol,
#         entrezgene,
#         amino_acid_position,
#         reference_amino_acid,
#         ttype,
#         freq) |>
#       dplyr::group_by(
#         symbol,
#         entrezgene,
#         amino_acid_position,
#         reference_amino_acid,
#         ttype) |>
#       dplyr::summarise(
#         ttype_site_freq = sum(as.integer(freq)),
#         .groups = "drop")
#   )
#
#   site_freqs[['indel']] <- as.data.frame(
#     cancer_hotspots[['indel']] |>
#       dplyr::select(
#         symbol,
#         entrezgene,
#         amino_acid_position,
#         ttype,
#         freq) |>
#       dplyr::group_by(
#         symbol,
#         entrezgene,
#         amino_acid_position,
#         ttype) |>
#       dplyr::summarise(
#         ttype_site_freq = sum(as.integer(freq)),
#         .groups = "drop")
#   )
#
#   for(t in c('snv','indel')){
#     cancer_hotspots[[t]] <- cancer_hotspots[[t]] |>
#       dplyr::left_join(site_freqs[[t]]) |>
#       dplyr::mutate(vartype = t)
#   }
#
#   hotspots_long <- dplyr::bind_rows(
#     cancer_hotspots$snv,
#     cancer_hotspots$indel
#   )
#
#   hotspots_wide <-  as.data.frame(
#     dplyr::bind_rows(
#       cancer_hotspots$snv,
#       cancer_hotspots$indel) |>
#       dplyr::mutate(
#         MUTATION_HOTSPOT_CANCERTYPE = paste(
#           ttype, ttype_site_freq, freq, sep="|"
#         )) |>
#       dplyr::group_by(
#         symbol,
#         entrezgene,
#         amino_acid_position,
#         reference_amino_acid,
#         vartype,
#         codon,
#         hgvsp,
#         hgvsp2,
#         MUTATION_HOTSPOT,
#         MUTATION_HOTSPOT2) |>
#       dplyr::summarise(
#         MUTATION_HOTSPOT_CANCERTYPE = paste(
#           sort(MUTATION_HOTSPOT_CANCERTYPE), collapse=","
#         ),
#         .groups = "drop")
#   )
#
#   return(list('wide' = hotspots_wide, 'long' = hotspots_long))
#
# }


# rename_hotspot_tumor_types <- function(df){
#
#   df <- df |>
#     dplyr::mutate(tumor_type_freq = stringr::str_to_title(tumor_type_freq)) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Adrenalgland","Adrenal_Gland"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Biliarytract","Biliary_Tract"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Bladder","Bladder@Urinary_Tract"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Blood","Myeloid"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Cnsbrain","CNS@Brain"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Bowel","Colon@Rectum"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Esophagusstomach","Esophagus@Stomach"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Ovaryfallopiantube","Ovary@Fallopian_Tube"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Headandneck","Head_and_Neck"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Ampullaofvater","Ampulla_of_Vater"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Vulvavagina","Vulva@Vagina"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Softtissue","Soft_Tissue"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Unk","Unknown"
#     )) |>
#     dplyr::mutate(tumor_type_freq = stringr::str_replace(
#       tumor_type_freq, "Lymph","Lymphoid"
#     ))
#
#   return(df)
# }
#
# resolve_gene_symbol <- function(df, gOncoX = NULL){
#
#   df <- df |>
#     dplyr::left_join(
#       dplyr::select(
#         gOncoX$basic$records, symbol, entrezgene),
#       by = c("hugo_symbol" = "symbol")
#     )
#
#   df1 <- df |>
#     dplyr::filter(!is.na(entrezgene)) |>
#     dplyr::rename(symbol = hugo_symbol)
#
#   df2 <- df |>
#     dplyr::filter(is.na(entrezgene)) |>
#     dplyr::select(-entrezgene) |>
#     dplyr::left_join(
#       dplyr::select(
#         gOncoX$alias, value, entrezgene),
#       by = c("hugo_symbol" = "value")
#     ) |>
#     dplyr::select(-hugo_symbol) |>
#     dplyr::filter(!is.na(entrezgene)) |>
#     dplyr::left_join(
#       dplyr::select(
#         gOncoX$basic$records, symbol, entrezgene),
#       by = "entrezgene")
#
#   df <- dplyr::bind_rows(
#     df1, df2
#   )
#
#   return(df)
#
#
# }

map_cancer_hotspots <- function(
    raw_maf = NULL){


  mapped_hotspot_mutations <- list()
  mapped_hotspot_mutations[['by_mutation']] <-
    data.frame()
  mapped_hotspot_mutations[['by_codon']] <-
    data.frame()

  raw_maf2 <- raw_maf |>
    dplyr::left_join(
      dplyr::select(
        cancerHotspots::cancer_hotspots$wide,
        entrezgene,
        hgvsp,
        MUTATION_HOTSPOT,
        MUTATION_HOTSPOT2,
        MUTATION_HOTSPOT_CANCERTYPE),
      by = c("Entrez_Gene_Id" = "entrezgene",
             "tmp_hgvsp" = "hgvsp"),
      relationship = "many-to-many"
    )

  mapped_hotspot_mutations[['by_mutation']] <-
    raw_maf2 |>
    dplyr::filter(!is.na(MUTATION_HOTSPOT))

  if(NROW(mapped_hotspot_mutations[['by_mutation']]) > 0){
    mapped_hotspot_mutations[['by_mutation']] <-
      mapped_hotspot_mutations[['by_mutation']] |>
      dplyr::select(
        Chromosome,
        Start_Position,
        End_Position,
        Variant_Type,
        Reference_Allele,
        Tumor_Seq_Allele2,
        tumor_sample_barcode,
        Entrez_Gene_Id,
        principal_hgvsp,
        MUTATION_HOTSPOT,
        MUTATION_HOTSPOT2,
        MUTATION_HOTSPOT_CANCERTYPE
      ) |>
      dplyr::distinct() |>
      dplyr::mutate(
        MUTATION_HOTSPOT_MATCH = "by_mutation"
      )

  }

  ## vars at mutation hotspots
  #maf_mutation_hotspots <- raw_maf |>
  #  dplyr::filter(!is.na(MUTATION_HOTSPOT))

  ## vars not at hotspots - map by codon as well
  mapped_hotspot_mutations[['by_codon']] <-
    raw_maf2 |>
    dplyr::filter(
      is.na(MUTATION_HOTSPOT) &
        (Variant_Classification == "Missense_Mutation" |
           Variant_Classification == "Nonsense_Mutation" |
           Variant_Classification == "Silent"))


  if(NROW(mapped_hotspot_mutations[['by_codon']]) > 0){
    mapped_hotspot_mutations[['by_codon']] <-
      mapped_hotspot_mutations[['by_codon']] |>
      dplyr::select(
        -c(MUTATION_HOTSPOT,MUTATION_HOTSPOT2,
           MUTATION_HOTSPOT_CANCERTYPE)
      ) |>
      dplyr::select(
        Chromosome,
        Start_Position,
        End_Position,
        Variant_Type,
        Reference_Allele,
        Tumor_Seq_Allele2,
        tumor_sample_barcode,
        Entrez_Gene_Id,
        principal_hgvsp,
        tmp_hgvsp) |>
      dplyr::mutate(
        codon = stringr::str_match(
          tmp_hgvsp, "p.[A-Z]{1}[0-9]{1,}")[,1]) |>
      dplyr::filter(!is.na(codon)) |>
      dplyr::rename(entrezgene = Entrez_Gene_Id) |>
      dplyr::left_join(
        dplyr::select(
          cancerHotspots::cancer_hotspots$wide,
          codon, entrezgene,
          MUTATION_HOTSPOT,
          MUTATION_HOTSPOT2,
          MUTATION_HOTSPOT_CANCERTYPE),
        by = c("codon","entrezgene"),
        relationship = "many-to-many") |>
      dplyr::filter(!is.na(MUTATION_HOTSPOT)) |>
      dplyr::select(-c(tmp_hgvsp, codon)) |>
      dplyr::distinct() |>
      dplyr::rename(Entrez_Gene_Id = entrezgene) |>
      dplyr::mutate(MUTATION_HOTSPOT_MATCH = "by_codon_only") |>
      tidyr::separate(
        MUTATION_HOTSPOT,
        into = c('tmp_sym','tmp_geneid','aapos',
                 'alt_aa','qvalue'),
        sep = '\\|',
        remove = T
      ) |>
      dplyr::mutate(
        MUTATION_HOTSPOT = paste(
          tmp_sym, tmp_geneid,
          aapos, '', qvalue, sep="|"
        )
      ) |>
      dplyr::select(
        -c(tmp_sym, tmp_geneid, aapos,
           alt_aa, qvalue)
      ) |>
      tidyr::separate(
        MUTATION_HOTSPOT2,
        into = c('tmp_sym','tmp_geneid','aapos',
                 'alt_aa','qvalue'),
        sep = "\\|",
        remove = T
      ) |>
      dplyr::mutate(
        MUTATION_HOTSPOT2 = paste(
          tmp_sym, tmp_geneid,
          aapos, '', qvalue, sep="|"
        )
      ) |>
      dplyr::select(
        -c(tmp_sym, tmp_geneid, aapos,
           alt_aa, qvalue)
      ) |>
      tidyr::separate_rows(
        MUTATION_HOTSPOT_CANCERTYPE, sep = ","
      ) |>
      dplyr::mutate(
        MUTATION_HOTSPOT_CANCERTYPE = stringr::str_replace(
          MUTATION_HOTSPOT_CANCERTYPE, "\\|[0-9]{1,}$",""
        )
      ) |>
      dplyr::distinct() |>
      dplyr::group_by(
        dplyr::across(c(-MUTATION_HOTSPOT_CANCERTYPE))) |>
      dplyr::summarise(
        MUTATION_HOTSPOT_CANCERTYPE = paste(
          sort(MUTATION_HOTSPOT_CANCERTYPE), collapse=","
        ), .groups = "drop"
      )

  }


  hotspot_map <-
    dplyr::bind_rows(
      mapped_hotspot_mutations[['by_mutation']],
      mapped_hotspot_mutations[['by_codon']]
    )


  ## some variants are mapped to multiple hotspots (transcript-specific)

  num_hotspots_per_variant <- hotspot_map |>
    dplyr::group_by(tumor_sample_barcode,
                    Chromosome,
                    Start_Position,
                    End_Position,
                    Reference_Allele,
                    Tumor_Seq_Allele2) |>
    dplyr::summarise(n = dplyr::n(),
                     .groups = "drop")

  principal_dup_hotspots <- num_hotspots_per_variant |>
    dplyr::filter(n > 1) |>
    dplyr::inner_join(
      hotspot_map,
      by = c("tumor_sample_barcode",
             "Chromosome",
             "Start_Position",
             "End_Position",
             "Reference_Allele",
             "Tumor_Seq_Allele2"),
      relationship = "many-to-many") |>
    dplyr::filter(principal_hgvsp == T)

  other_hotspots <- num_hotspots_per_variant |>
    dplyr::filter(n == 1) |>
    dplyr::inner_join(
      hotspot_map,
      by = c("tumor_sample_barcode",
             "Chromosome",
             "Start_Position",
             "End_Position",
             "Reference_Allele",
             "Tumor_Seq_Allele2"),
      relationship = "many-to-many")

  hotspot_map_nonredundant <-
    dplyr::bind_rows(
      principal_dup_hotspots,
      other_hotspots
    ) |>
    dplyr::select(-principal_hgvsp)

  if(NROW(hotspot_map_nonredundant) > 0){
    raw_maf <- raw_maf |>
      dplyr::left_join(
        hotspot_map_nonredundant,
        by = c(
          "Chromosome",
          "Start_Position",
          "End_Position",
          "Variant_Type",
          "Reference_Allele",
          "Tumor_Seq_Allele2",
          "tumor_sample_barcode",
          "Entrez_Gene_Id"
        ), relationship = "many-to-many")
  }else{
    raw_maf <- raw_maf |>
      dplyr::mutate(
        MUTATION_HOTSPOT = NA,
        MUTATION_HOTSPOT2 = NA,
        MUTATION_HOTSPOT_CANCERTYPE = NA,
        MUTATION_HOTSPOT_MATCH = NA,
      )
  }

  return(raw_maf)

}
