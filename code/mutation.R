suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

get_tcga_snv <- function(
    gdc_projects = NA,
    tcga_clinical_info = NA,
    tcga_release = NA,
    data_raw_dir = NA,
    overwrite = FALSE,
    gOncoX = NA,
    gdc_files_per_chunk = 50,
    output_dir = NA){

  tcga_full_rds <-
    file.path(
      output_dir,
      tcga_release,
      "snv_indel",
      "tcga_mutation_grch38.rds"
    )

  if(file.exists(tcga_full_rds) & overwrite == FALSE){
    tcga_snv_calls <- readRDS(file = tcga_full_rds)
    return(tcga_snv_calls)
  }

  tcga_driver_mutations <- get_tcga_driver_mutations(
    data_raw_dir = data_raw_dir
  )

  docm_mutations <- get_curated_docm_mutations(
    data_raw_dir = data_raw_dir
  )

  tcga_clinical <- tcga_clinical_info[['slim']] |>
    dplyr::filter(!is.na(primary_site) &
                    primary_site != 'Other/Unknown')

  if(!dir.exists(
    file.path(output_dir, tcga_release))){
    dir.create(
      file.path(output_dir, tcga_release))
  }

  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "snv_indel"))){
    dir.create(
      file.path(output_dir, tcga_release, "snv_indel"))
  }

  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = GenomeInfoDb::seqlevels(
      GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)),
    seqlengths = GenomeInfoDb::seqlengths(
      GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)),
    genome = 'hg38')


  j <- 0
  tcga_calls <- data.frame()

  ucsc_url <-
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"

  repeat_tracks <- list()
  for(e in c('simpleRepeat','windowmaskerSdust')){
    repeat_tracks[[e]] <-
      GenomicRanges::makeGRangesFromDataFrame(
        as.data.frame(readr::read_tsv(
          file = file.path(
            ucsc_url,
            paste0(e,".txt.gz")),
          col_names = c('bin','chrom','start','end','name'),
          show_col_types = F) |>
          dplyr::select(chrom,start,end) |>
          dplyr::filter(!stringr::str_detect(chrom,"_")) |>
          dplyr::mutate(chrom = stringr::str_replace(chrom,"chr","")) |>
          dplyr::distinct()
        ),
        keep.extra.columns = TRUE
      )
  }


  for(i in 1:nrow(gdc_projects)){

    project_code <- gdc_projects[i,]$project_id
    tumor <- gdc_projects[i,]$tumor

    cat(i," - ",gdc_projects[i,]$project_id,'\n')

    rds_file <- file.path(
      output_dir,
      tcga_release,
      "snv_indel",
      paste0(
        "tcga_mutation_",
        tumor,
        "_grch38",
        ".rds"
      )
    )

    project_calls <- NULL
    if(file.exists(rds_file) & overwrite == F){
      project_calls <- readRDS(rds_file)
    }else{

      query <- suppressMessages(TCGAbiolinks::GDCquery(
        project = project_code,
        data.category = "Simple Nucleotide Variation",
        access = "open",
        data.type = "Masked Somatic Mutation",
        workflow.type =
          "Aliquot Ensemble Somatic Variant Merging and Masking"
      ))

      TCGAbiolinks::GDCdownload(
        query,
        directory = "GDCdata",
        method = "api",
        files.per.chunk = gdc_files_per_chunk)

      raw_maf <- as.data.frame(
        TCGAbiolinks::GDCprepare(query, directory = "GDCdata")) |>
        dplyr::select(Chromosome,
                      Start_Position,
                      End_Position,
                      Variant_Type,
                      Reference_Allele,
                      Tumor_Seq_Allele2,
                      Tumor_Sample_Barcode,
                      callers,
                      t_depth,
                      t_ref_count,
                      t_alt_count,
                      n_depth,
                      NCBI_Build,
                      Strand,
                      Variant_Classification,
                      Gene,
                      HGNC_ID,
                      Entrez_Gene_Id,
                      SYMBOL,
                      Feature,
                      Feature_type,
                      One_Consequence,
                      Consequence,
                      all_effects,
                      HGVSc,
                      HGVSp_Short) |>
        dplyr::rename(Algorithms = callers,
                      All_variant_effects = all_effects) |>
        dplyr::mutate(Reference_Allele = dplyr::if_else(
          Reference_Allele == "TRUE",
          as.character("T"),
          as.character(Reference_Allele)
        )) |>
        dplyr::mutate(Tumor_Seq_Allele2 = dplyr::if_else(
          Tumor_Seq_Allele2 == "TRUE",
          as.character("T"),
          as.character(Tumor_Seq_Allele2)
        )) |>
        dplyr::mutate(tumor_sample_barcode = stringr::str_extract(
          Tumor_Sample_Barcode,
          "TCGA-([:alnum:]){2}-([:alnum:]){4}-0[0-9][A-Z]")) |>
        dplyr::select(-c(Tumor_Sample_Barcode)) |>
        dplyr::group_by(
          dplyr::across(-c(t_depth, t_ref_count,
                           t_alt_count, n_depth))
        ) |>
        dplyr::summarise(
          t_depth = floor(median(t_depth)),
          t_ref_count = floor(median(t_ref_count)),
          t_alt_count = floor(median(t_alt_count)),
          n_depth = floor(median(n_depth)),
          #num = dplyr::n(),
          .groups = "drop"
        ) |>
        tidyr::separate_rows(All_variant_effects, sep=";") |>
        tidyr::separate(All_variant_effects,
                        into = c("tmp_genesymbol",
                                 "tmp_consequence",
                                 "tmp_hgvsp",
                                 "tmp_ensembl_trans",
                                 "tmp_refseq_trans",
                                 "tmp_hgvsc",
                                 "tmp_impact",
                                 "tmp_canonical",
                                 "tmp_sift",
                                 "tmp_polyphen",
                                 "tmp_unknown"),
                        sep = ",") |>
        dplyr::mutate(principal_hgvsp = dplyr::if_else(
          tmp_hgvsp == HGVSp_Short,
          TRUE,
          FALSE
        ))


      raw_maf_drivers1 <- raw_maf |>
        dplyr::left_join(tcga_driver_mutations[['by_symbol']],
                         by = c("SYMBOL","HGVSp_Short")) |>
        dplyr::filter(!is.na(driver_mutation)) |>
        dplyr::distinct()

      raw_maf_drivers2 <- raw_maf |>
        dplyr::left_join(tcga_driver_mutations[['by_symbol']],
                         by = c("SYMBOL","HGVSp_Short")) |>
        dplyr::filter(is.na(driver_mutation)) |>
        dplyr::select(-driver_mutation) |>
        dplyr::left_join(tcga_driver_mutations[['by_transcript_id']],
                         by = c("HGVSp_Short","Feature")) |>
        dplyr::distinct()


      raw_maf <- as.data.frame(
        dplyr::bind_rows(
          raw_maf_drivers1,
          raw_maf_drivers2
        )
      )

      raw_maf_with_hotspots <- as.data.frame(map_cancer_hotspots(
        raw_maf = raw_maf)
      )

      if(NROW(raw_maf) != NROW(raw_maf_with_hotspots)){
        cat("WARNING: non-unique variant set")
      }

      raw_maf <- as.data.frame(raw_maf_with_hotspots |>
        dplyr::mutate(Alternate_effects = paste(
          tmp_genesymbol, tmp_consequence,
          tmp_hgvsp, tmp_ensembl_trans, sep=":"
        )) |>
        dplyr::group_by(Chromosome,
                      Start_Position,
                      End_Position,
                      Variant_Type,
                      Reference_Allele,
                      Tumor_Seq_Allele2,
                      tumor_sample_barcode,
                      Algorithms,
                      t_depth,
                      t_ref_count,
                      t_alt_count,
                      n_depth,
                      NCBI_Build,
                      Strand,
                      Variant_Classification,
                      Gene,
                      HGNC_ID,
                      Entrez_Gene_Id,
                      SYMBOL,
                      Feature,
                      Feature_type,
                      One_Consequence,
                      Consequence,
                      HGVSc,
                      HGVSp_Short,
                      driver_mutation,
                      MUTATION_HOTSPOT,
                      MUTATION_HOTSPOT2,
                      MUTATION_HOTSPOT_CANCERTYPE,
                      MUTATION_HOTSPOT_MATCH) |>
        dplyr::summarise(
          Alternate_effects = paste(
            sort(unique(Alternate_effects)), collapse=";"),
          .groups = "drop")
      )

      #raw_maf <- as.data.frame(raw_maf2)
      coding_csq_pattern <-
        paste0(
          "^(stop_|start_lost|frameshift_|missense_variant|",
          "splice_donor|splice_acceptor|protein_altering|inframe_)")
      exonic_csq_pattern <-
        paste0(
          "^(stop_|start_lost|frameshift_|missense_variant|splice_donor",
          "|splice_acceptor|synonymous|protein_altering|inframe_)")

      project_calls <- get_proper_maf_alleles(
        maf_df = raw_maf,
        genome_seq = genome_seq,
        seqinfo = seqinfo) |>
        dplyr::left_join(
          docm_mutations, by = c("CHROM","POS","REF","ALT")) |>
        dplyr::mutate(CODING_STATUS = dplyr::if_else(
          stringr::str_detect(
          One_Consequence,coding_csq_pattern),
          "coding","noncoding","noncoding")) |>
        dplyr::mutate(EXONIC_STATUS = dplyr::if_else(stringr::str_detect(
          One_Consequence,exonic_csq_pattern),
          "exonic","nonexonic","nonexonic")) |>
        dplyr::mutate(bcr_patient_barcode = stringr::str_replace(
          tumor_sample_barcode,
          "-0[1-9]{1}[A-Z]{1}$","")) |>
        dplyr::mutate(sample_type = dplyr::if_else(
          stringr::str_detect(
            tumor_sample_barcode,
            "-0(1|5)[A-Z]$"),"Solid Tumor - Primary",as.character(NA))) |>
        dplyr::mutate(sample_type = dplyr::if_else(
          stringr::str_detect(
            tumor_sample_barcode,"-02[A-Z]$"),
          "Solid Tumor - Recurrent",as.character(sample_type))) |>
        dplyr::mutate(sample_type = dplyr::if_else(
          stringr::str_detect(
            tumor_sample_barcode,"-0(6|7)[A-Z]$"),
          "Metastatic",as.character(sample_type))) |>
        dplyr::mutate(sample_type = dplyr::if_else(
          stringr::str_detect(
            tumor_sample_barcode,"-0(3|4|9)[A-Z]$"),
          "Blood-Derived Cancer",as.character(sample_type))) |>
        dplyr::rename(Hugo_Symbol = SYMBOL) |>
        dplyr::mutate(tumor = tumor) |>
        dplyr::left_join(
          dplyr::select(tcga_clinical,
                        bcr_patient_barcode,
                        site_diagnosis_code,
                        primary_site,
                        primary_diagnosis,
                        primary_diagnosis_very_simplified),
          by = "bcr_patient_barcode",
          relationship = "many-to-many") |>
        dplyr::distinct() |>
        dplyr::mutate(
          REF = as.character(REF),
          ALT = as.character(ALT)
        ) |>
        dplyr::select(bcr_patient_barcode,
                      tumor_sample_barcode,
                      sample_type,
                      primary_site,
                      primary_diagnosis,
                      primary_diagnosis_very_simplified,
                      site_diagnosis_code,
                      CHROM,
                      POS,
                      REF,
                      ALT,
                      GENOMIC_CHANGE,
                      dplyr::everything()
                )

      calls_simple <-
        dplyr::select(project_calls, CHROM, POS, ALT) |>
        dplyr::mutate(end = POS) |>
        dplyr::rename(chrom = CHROM,
                      start = POS) |>
        dplyr::select(chrom, start, end)

      gr_calls <- GenomicRanges::makeGRangesFromDataFrame(
        calls_simple,
        keep.extra.columns = TRUE
      )

      repeat_hits <- list()
      for(e in c('simpleRepeat',
                 'windowmaskerSdust')){

        hits <- GenomicRanges::findOverlaps(
          gr_calls, repeat_tracks[[e]],
          type = "any", select = "all")

        var_hits <-
          unique(S4Vectors::queryHits(hits))

        if(e == 'simpleRepeat'){
          project_calls$SIMPLEREPEATS_HIT <- FALSE
          if(length(var_hits) > 0){
            project_calls[var_hits,]$SIMPLEREPEATS_HIT = TRUE
          }
        }else{
          project_calls$WINMASKER_HIT <- FALSE
          if(length(var_hits) > 0){
            project_calls[var_hits,]$WINMASKER_HIT = TRUE
          }
        }

      }

      project_calls <- project_calls |>
        dplyr::mutate(WINMASKER_HIT = dplyr::if_else(
          is.na(WINMASKER_HIT),
          as.logical(FALSE),
          as.logical(WINMASKER_HIT)
        )) |>
        dplyr::mutate(SIMPLEREPEATS_HIT = dplyr::if_else(
          is.na(SIMPLEREPEATS_HIT),
          as.logical(FALSE),
          as.logical(SIMPLEREPEATS_HIT)
        ))

      saveRDS(project_calls, file=rds_file)
    }
    tcga_calls <- dplyr::bind_rows(
      tcga_calls, project_calls)

  }

  tcga_calls_simple <- tcga_calls |>
    dplyr::select(CHROM, POS, REF, ALT, GENOMIC_CHANGE,
                  tumor_sample_barcode, tumor) |>
    dplyr::distinct()

  tcga_samples_pr_cohort <- as.data.frame(
    dplyr::select(tcga_calls_simple,
                  tumor_sample_barcode, tumor) |>
      dplyr::distinct()  |>
      dplyr::group_by(tumor) |>
      dplyr::summarise(COHORT_SIZE = dplyr::n(),
                       .groups = "drop")
  )

  ## get variant recurrence pr tumor type, and pancancer
  tcga_calls_simple_recurrence <- as.data.frame(
    dplyr::group_by(tcga_calls_simple,
                    GENOMIC_CHANGE, tumor) |>
      dplyr::summarise(AFFECTED = dplyr::n(),
                       .groups = "drop")
  )

  tcga_calls_simple_recurrence_pancancer <- as.data.frame(
    dplyr::group_by(tcga_calls_simple, GENOMIC_CHANGE) |>
      dplyr::summarise(pancancer_frequency = dplyr::n(),
                       .groups = "drop")
  )

  tcga_calls_recurrence_stats <- as.data.frame(
    tcga_calls_simple_recurrence |>
      dplyr::left_join(
        tcga_samples_pr_cohort, by = "tumor") |>
      dplyr::mutate(
        PERCENT =
          round(as.numeric(AFFECTED /COHORT_SIZE) * 100, digits = 1)) |>
      dplyr::mutate(
        TCGA_FREQUENCY =
          paste(tumor, PERCENT, AFFECTED,
                COHORT_SIZE, sep = "|")) |>
      dplyr::group_by(GENOMIC_CHANGE) |>
      dplyr::summarise(
        cohort_frequency =
          paste(TCGA_FREQUENCY, collapse=","), .groups = "drop") |>
      dplyr::left_join(
        tcga_calls_simple_recurrence_pancancer,
        by = "GENOMIC_CHANGE")
  )

  tcga_calls <- tcga_calls |>
    dplyr::left_join(tcga_calls_recurrence_stats,
                     by = "GENOMIC_CHANGE")

  saveRDS(
    tcga_calls,
    file = tcga_full_rds
  )

  readr::write_tsv(
    tcga_calls,
    file =
      file.path(
        output_dir,
        tcga_release,
        "snv_indel",
        "tcga_mutation_grch38.tsv.gz"
      )
  )

  # system(
  #   paste0(
  #     'gzip -f ',
  #     file.path(
  #       output_dir,
  #       tcga_release,
  #       "snv_indel",
  #       "tcga_mutation_grch38.tsv"
  #     )
  #   )
  # )


  ## Create MAF files per TCGA cohort
  for(site in unique(tcga_calls$primary_site)){
    if(is.na(site) | site == "Other/Unknown"){
      next
    }
    tcga_maf <- tcga_calls |>
      dplyr::filter(primary_site == site) |>
      dplyr::select(Hugo_Symbol,
                    Chromosome,
                    Start_Position,
                    End_Position,
                    Variant_Type,
                    Reference_Allele,
                    Tumor_Seq_Allele2,
                    Variant_Classification,
                    NCBI_Build,
                    tumor_sample_barcode,
                    site_diagnosis_code,
                    primary_site) |>
      dplyr::rename(Tumor_Sample_Barcode = tumor_sample_barcode) |>
      dplyr::distinct() |>
      dplyr::filter(!is.na(
        site_diagnosis_code
      ))

    site_code <- stringr::str_replace(
      unique(tcga_maf$site_diagnosis_code),"_[0-9]","")[1]
    header_line_maf <- c("#version 2.4")
    maf_fname_all <-
      file.path(
        output_dir,
        tcga_release,
        "snv_indel",
        paste0(
          "tcga_mutation_grch38_",
          site_code,
          "_0.maf"
        )
      )
    write(header_line_maf,
          file = maf_fname_all,sep="\n")
    write.table(
      tcga_maf, file = maf_fname_all,
      na="", quote=F, row.names = F,
      col.names = T, sep="\t")
    system(paste0("gzip -f ", maf_fname_all))
    cat(site,'\n')

  }
}


write_tcga_vcf <- function(
    tcga_calls = NA,
    output_dir = NA,
    tcga_release = NA){


  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "vcf"))){
    dir.create(
      file.path(output_dir, tcga_release, "vcf"))
  }

  header_lines <- c(
    "##fileformat=VCFv4.2",
    paste0("##SOURCE_TCGA=",tcga_release),
    "##INFO=<ID=TCGA_SAMPLES,Number=.,Type=String,Description=\"TCGA barcodes\">",
    "##INFO=<ID=TCGA_PANCANCER_COUNT,Number=1,Type=Integer,Description=\"Raw variant count across all tumor types\">",
    "##INFO=<ID=TCGA_FREQUENCY,Number=.,Type=String,Description=\"Frequency of variant across TCGA cancer subtypes. Format: subtype|percent affected|affected cases|total cases\">",
    "##INFO=<ID=TCGA_DRIVER_MUTATION,Number=.,Type=String,Description=\"Putative cancer driver mutation discovered by multiple approaches (Bailey et al., Cell, 2018). Format: symbol:hgvsp:ensembl_transcript_ids:discovery_approaches\">",
    "##INFO=<ID=WINMASKER_HIT,Number=0,Type=Flag,Description=\"Overlap with UCSC repeat track windowmaskerSdust\">",
    "##INFO=<ID=SIMPLEREPEATS_HIT,Number=0,Type=Flag,Description=\"Overlap with UCSC simpleRepeats track\">",
    "##INFO=<ID=DOCM_DISEASE,Number=.,Type=String,Description=\"Disease associated with variant as found in the DoCM database\">",
    "##INFO=<ID=DOCM_ID,Number=.,Type=String,Description=\"DoCM variant identifier (HGVS)\">",
    "##INFO=<ID=DOCM_PMID,Number=.,Type=String,Description=\"PubMed IDs supporting variant-disease association in DoCM database\">")

  write(
    header_lines,
    file = file.path(
      output_dir,
      tcga_release,
      "vcf",
      "tcga.grch38.vcf"
    ), sep = "\n")


  tcga_samples_per_var <-
    as.data.frame(
      tcga_calls |>
        dplyr::group_by(
          CHROM, POS, REF, ALT) |>
        dplyr::summarise(TCGA_SAMPLES = paste(
          tumor_sample_barcode,collapse=","),
          .groups = "drop") |>
        dplyr::ungroup()
    )

  tcga_vcf_records <- as.data.frame(
    tcga_calls |>
      dplyr::left_join(
        tcga_samples_per_var,
        by = c("CHROM","POS","REF","ALT")) |>
      dplyr::mutate(driver_mutation = dplyr::if_else(
        is.na(driver_mutation),
        ".",
        as.character(driver_mutation)
      )) |>
      dplyr::mutate(
        QUAL = ".",
        FILTER = 'PASS',
        ID = '.') |>
      dplyr::mutate(winmasker_hit = dplyr::if_else(
        WINMASKER_HIT == F,
        "",
        as.character("WINMASKER_HIT")
      )) |>
      dplyr::mutate(simplerepeats_hit = dplyr::if_else(
        SIMPLEREPEATS_HIT == F,
        "",
        as.character("SIMPLEREPEATS_HIT")
      )) |>
      dplyr::mutate(INFO = paste0(
        "TCGA_SAMPLES=",TCGA_SAMPLES,
        ";TCGA_FREQUENCY=",cohort_frequency,
        ";TCGA_PANCANCER_COUNT=",pancancer_frequency,
        ";TCGA_DRIVER_MUTATION=", driver_mutation,
        ";DOCM_DISEASE=", docm_disease,
        ";DOCM_ID=", docm_id,
        ";DOCM_PMID=", docm_pmid,
        ";", winmasker_hit,
        ";", simplerepeats_hit)) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";TCGA_DRIVER_MUTATION=\\.","")) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";DOCM_DISEASE=NA","")) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";DOCM_ID=NA","")) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";DOCM_PMID=NA","")) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";{1,}$", ""
      )) |>
      dplyr::select(
        CHROM, POS, ID, REF,
        ALT, QUAL, FILTER, INFO) |>
      dplyr::distinct()
  )

  ## Write to VCF (temporary)
  vcfhelpR::write_vcf_records(
    vcf_records = tcga_vcf_records,
    output_dir =
      file.path(
        output_dir,
        tcga_release,
        "vcf"),
    vcf_fname_prefix = "tcga",
    header_lines = header_lines,
    genome_build = "grch38"
  )

  ## CrossMap to grch37
  vcfhelpR::crossmap_vcf(
    source_vcf =
      file.path(
        output_dir,
        tcga_release,
        "vcf",
        "tcga_grch38.vcf.gz"
      ),
    target_vcf =
      file.path(
        output_dir,
        tcga_release,
        "vcf",
        "tcga_grch37.vcf"
      ),
  )

  write(c(header_lines[2],
          header_lines[4:6],
          header_lines[9:11]),
        file = file.path(
          output_dir,
          tcga_release,
          "vcf",
          "tcga.vcfanno.vcf_info_tags.txt"
        ), sep = "\n")

}


calculate_sample_tmb <- function(
    tcga_calls = NA,
    tcga_release = NA,
    tcga_clinical_info = NA,
    output_dir = NA,
    overwrite = FALSE,
    exome_target_size_mb = 34.0){

  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "tmb"))){
    dir.create(
      file.path(
        output_dir, tcga_release, "tmb"))
  }

  tcga_tmb_rds <-
    file.path(
      output_dir,
      tcga_release,
      "tmb",
      "tcga_tmb.rds"
    )

  if(file.exists(tcga_tmb_rds) & overwrite == FALSE){
    tcga_tmb <- readRDS(file = tcga_tmb_rds)
    return(tcga_tmb)
  }

  ## default target size region of protein-coding portion for TCGA WES samples
  ## Approximate solution (34Mb) based on coding content of GENCODE protein-coding exons (grch38/v33)
  ## Notably, Most TCGA samples are apparently sequenced with Agilent Custom V2 Exome Bait (Broad Institute)
  ## (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6169918/), not been able to locate BED file for this target
  #exome_target_size_mb <- as.numeric(34.0)

  ## Approaches for TMB calculation (proxy for neoantigen load in tumor)
  # 1. Chalmers et al. Genome Medicine, 2017 - all coding mutations (SNVs + InDels), also including synonymous mutations
  # 2. Fernandez et al. JCO Precision Oncology, 2019 - all non-synonymous mutations only
  # 3. VEP - HIGH/MODERATE - Chalmers except for silent mutations


  vep_high_moderate_csq_pattern <- paste0(
    "^(stop_lost|start_lost|stop_gained|frameshift_|",
    "missense_|splice_donor|splice_acceptor|protein_altering|inframe_)"
  )


  tcga_tmb_measures <- list()

  tcga_tmb_measures[['chalmers']] <- as.data.frame(
    tcga_calls |>
      dplyr::filter(EXONIC_STATUS == "exonic") |>
      dplyr::select(GENOMIC_CHANGE, bcr_patient_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode) |>
      dplyr::summarise(n_coding_mutations = dplyr::n(),
                       .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(tmb = round(
        as.numeric(n_coding_mutations) /
          exome_target_size_mb,digits = 5)) |>
      dplyr::mutate(tmb_log10 = log10(tmb))
  )

  tcga_tmb_measures[['vep_high_moderate']] <- as.data.frame(
    tcga_calls |>
      dplyr::filter(stringr::str_detect(
        Consequence, vep_high_moderate_csq_pattern)) |>
      dplyr::select(GENOMIC_CHANGE, bcr_patient_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode) |>
      dplyr::summarise(n_coding_mutations_vep_high_moderate = dplyr::n(),
                       .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(tmb_vep_high_moderate = round(
        as.numeric(n_coding_mutations_vep_high_moderate) /
          exome_target_size_mb,digits = 5)) |>
      dplyr::mutate(tmb_vep_high_moderate_log10 = log10(tmb_vep_high_moderate))
  )

  tcga_tmb_measures[['fernandez']] <- as.data.frame(
    tcga_calls |>
      dplyr::filter(startsWith(Consequence,"missense_variant")) |>
      dplyr::select(GENOMIC_CHANGE, bcr_patient_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode) |>
      dplyr::summarise(n_nsyn_mutations = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(tmb_ns = round(
        as.numeric(n_nsyn_mutations) /
          exome_target_size_mb,digits = 5)) |>
      dplyr::mutate(tmb_ns_log10 = log10(tmb_ns))
  )

  tcga_tmb <- tcga_tmb_measures[['chalmers']] |>
    dplyr::left_join(tcga_tmb_measures[['fernandez']],
                     by = "bcr_patient_barcode") |>
    dplyr::left_join(tcga_tmb_measures[['vep_high_moderate']],
                     by = "bcr_patient_barcode") |>
    dplyr::filter(!is.na(tmb_ns)) |>
    dplyr::left_join(dplyr::select(
      tcga_clinical_info$slim, tumor,
      primary_site,
      primary_diagnosis_simplified,
      primary_diagnosis_very_simplified,
      MSI_status, year_of_birth,
      age_at_index, gender,
      ER_status, PR_status,
      HER2_status, pancan_subtype_selected,
      Gleason_score,
      bcr_patient_barcode),
      by = "bcr_patient_barcode") |>
    dplyr::mutate(cds_target_size_mb = exome_target_size_mb)


  saveRDS(
    tcga_tmb,
    file = file.path(
      output_dir,
      tcga_release,
      "tmb",
      "tcga_tmb.rds"
    )
  )

  readr::write_tsv(
    tcga_tmb,
    file = file.path(
      output_dir,
      tcga_release,
      "tmb",
      "tcga_tmb.tsv.gz"
    )
  )

  return(tcga_tmb)
}

calculate_gene_mutation_rate <- function(
    output_dir = NA,
    tcga_calls = NA,
    tcga_clinical_info = NA,
    tcga_release = NA,
    overwrite = FALSE){

  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "gene"))){
    dir.create(
      file.path(
        output_dir, tcga_release, "gene"))
  }

  tcga_gene_mut_rate_rds <-
    file.path(
      output_dir,
      tcga_release,
      "gene",
      "tcga_gene_aberration_rate.rds"
    )

  if(file.exists(tcga_gene_mut_rate_rds) &
     overwrite == FALSE){
    tcga_gene_mut_rate <-
      readRDS(file = tcga_gene_mut_rate_rds)
    return(tcga_gene_mut_rate)
  }

  cna_calls <- readRDS(
    file = file.path(
      output_dir,
      tcga_release,
      "cna",
      "tcga_cna_gistic2.ampl_homdel.rds")
  )

  tcga_calls <- tcga_calls |>
    dplyr::rename(symbol = Hugo_Symbol)

  cohort_stats1 <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "cna_homdel",
    clinical_strata = "site") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats2 <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "cna_homdel",
    clinical_strata = "site_diagnosis") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats3 <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "cna_ampl",
    clinical_strata = "site") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats4 <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "cna_ampl",
    clinical_strata = "site_diagnosis") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats5 <- cohort_mutation_stats(
    calls_df = tcga_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "snv_indel",
    clinical_strata = "site",
    genomic_strata = "gene") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats7 <- cohort_mutation_stats(
    calls_df = tcga_calls,
    clinical_df = tcga_clinical_info$slim,
    vartype = "snv_indel",
    clinical_strata = "site_diagnosis",
    genomic_strata = "gene") |>
    dplyr::filter(!is.na(symbol))


  tcga_gene_stats <- dplyr::bind_rows(
    cohort_stats7, cohort_stats5,
    cohort_stats4, cohort_stats3,
    cohort_stats2, cohort_stats1)

  saveRDS(
    tcga_gene_stats,
    file = tcga_gene_mut_rate_rds
  )

  return(tcga_gene_stats)
}
