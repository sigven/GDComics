suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

#' Get or make TCGA SNV/InDel data from GDC
#'
#' @param output_dir Directory to save/load processed SNV/InDel data
#' @param gencode_xref Data frame with GENCODE gene identifiers
#' @param data_raw_dir Directory containing raw GDC data
#' @param gdc_projects Character - GDC/TCGA project IDs
#' @param gdc_release TCGA data release version
#' @param overwrite Whether to overwrite existing processed data
#'
#' @export
#'
gdc_tcga_snv <- function(
    output_dir = NULL,
    gencode_xref = NULL,
    data_raw_dir = NULL,
    gdc_projects = NULL,
    gdc_release = "release45_20251204",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    !is.null(gdc_release)
  )

  tcga_full_rds <-
    file.path(
      output_dir,
      gdc_release,
      "snv_indel",
      "tcga_mutation_grch38.rds"
    )

  if(file.exists(tcga_full_rds) & overwrite == FALSE){
    tcga_snv_calls <- readRDS(file = tcga_full_rds)
    return(tcga_snv_calls)
  }

  tcga_clinical <- gdc_tcga_clinical(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir) |>
    dplyr::select(
      bcr_patient_barcode,
      tumor_sample_barcode,
      tumor,
      site_diagnosis_code,
      sample_type,
      primary_site,
      primary_diagnosis,
      primary_diagnosis_very_simplified
    )

  if(!dir.exists(
    file.path(
      output_dir, gdc_release, "snv_indel"))){
    dir.create(
      file.path(output_dir, gdc_release, "snv_indel"))
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

  assertthat::assert_that(
    !is.null(gdc_projects),
    !is.null(data_raw_dir)
  )
  assertthat::assert_that(
    dir.exists(data_raw_dir)
  )

  i <- 1
  for(project_code in gdc_projects){
    tumor <- stringr::str_replace(project_code, "TCGA-","")

    cat(i," - ",project_code,'\n')

    rds_file <- file.path(
      output_dir,
      gdc_release,
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

      raw_project_maf <- as.data.frame(load_ssm_maf(
        gdc_maf_cache_path =
          file.path(
            data_raw_dir,
            "GDCdata",
            "maf"
          ),
        gdc_projects = project_code))

      if(NROW(raw_project_maf) == 0){
        next
      }

      raw_maf <- raw_project_maf |>
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
                      HGVSp_Short,
                      RNA_Support,
                      RNA_depth,
                      RNA_ref_count,
                      RNA_alt_count,
                      gdc_maf_fname,
                      project_id) |>
        dplyr::rename(Algorithms = callers,
                      All_variant_effects = all_effects) |>
        dplyr::mutate(Start_Position = as.numeric(Start_Position),
                      End_Position = as.numeric(End_Position),
                      Entrez_Gene_Id = as.integer(Entrez_Gene_Id),
                      t_depth = as.integer(t_depth),
                      n_depth = as.integer(n_depth),
                      t_ref_count = as.integer(t_ref_count),
                      t_alt_count = as.integer(t_alt_count),
                      RNA_depth = as.integer(RNA_depth),
                      RNA_ref_count = as.integer(RNA_ref_count),
                      RNA_alt_count = as.integer(RNA_alt_count)) |>
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
        dplyr::mutate(bcr_patient_barcode = stringr::str_replace(
          tumor_sample_barcode,
          "-0[1-9]{1}[A-Z]{1}$","")) |>
        dplyr::arrange(
          .data$Chromosome,
          .data$Start_Position,
          .data$End_Position,
          .data$Reference_Allele,
          .data$Tumor_Seq_Allele2,
          .data$tumor_sample_barcode,
          dplyr::desc(t_depth)
        ) |>
        dplyr::group_by(
          dplyr::across(-c("t_depth", "t_ref_count","gdc_maf_fname",
                           "t_alt_count","n_depth","Algorithms"))
        ) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
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

      raw_maf_with_hotspots <- as.data.frame(map_cancer_hotspots(
        raw_maf = raw_maf)
      )

      if(NROW(raw_maf) != NROW(raw_maf_with_hotspots)){
        cat("WARNING: non-unique variant set")
      }

      raw_maf <- as.data.frame(
        raw_maf_with_hotspots |>
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
                          bcr_patient_barcode,
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
                          RNA_Support,
                          RNA_depth,
                          RNA_ref_count,
                          RNA_alt_count,
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
          "^(stop_(gained|lost)|start_lost|frameshift_|missense_variant|",
          "splice_donor|splice_acceptor|protein_altering|inframe_)")
      exonic_csq_pattern <-
        paste0(
          "^(stop_|start_lost|frameshift_|missense_variant|splice_donor|",
          "splice_acceptor|synonymous|splice_region_variant;synonymous|",
          "protein_altering|inframe_)")

      project_calls <- get_proper_maf_alleles(
        maf_df = raw_maf,
        genome_seq = genome_seq,
        seqinfo = seqinfo) |>
        dplyr::mutate(CODING_STATUS = dplyr::if_else(
          stringr::str_detect(
            One_Consequence,coding_csq_pattern),
          "coding","noncoding","noncoding")) |>
        dplyr::mutate(EXONIC_STATUS = dplyr::if_else(stringr::str_detect(
          One_Consequence,exonic_csq_pattern),
          "exonic","nonexonic","nonexonic")) |>
        dplyr::rename(Hugo_Symbol = SYMBOL) |>
        dplyr::inner_join(
          tcga_clinical,
          by = c("bcr_patient_barcode",
                 "tumor_sample_barcode")
        ) |>
        dplyr::distinct() |>
        dplyr::mutate(
          REF = as.character(REF),
          ALT = as.character(ALT)
        ) |>
        dplyr::select(bcr_patient_barcode,
                      tumor_sample_barcode,
                      sample_type,
                      tumor,
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

    i <- i + 1
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
        gdc_release,
        "snv_indel",
        "tcga_mutation_grch38.tsv.gz"
      )
  )

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
        gdc_release,
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
  return(tcga_calls)

}


write_tcga_sample_vcf <- function(
    tcga_calls = NA,
    output_dir = NA,
    tcga_release = NA){

    if(!dir.exists(
      file.path(
        output_dir, tcga_release, "vcf"))){
      dir.create(
        file.path(output_dir, tcga_release, "vcf"))
    }
    if(!dir.exists(
      file.path(
        output_dir, tcga_release, "vcf","per_sample"))){
      dir.create(
        file.path(output_dir, tcga_release, "vcf","per_sample"))
    }

  all_sample_calls <- tcga_calls |>
    dplyr::mutate(
      ALGORITHMS = stringr::str_replace_all(Algorithms,";",","),
      SAMPLE_TYPE = tolower(stringr::str_replace_all(
        stringr::str_replace_all(sample_type," - ","-"), " |-", "_")),
      TVAF = round(t_alt_count/(t_alt_count + t_ref_count), digits = 3),
      TDP = t_alt_count + t_ref_count,
      CDP = n_depth) |>
    dplyr::select(CHROM, POS, REF, ALT, tumor, ALGORITHMS, SAMPLE_TYPE,
                  TVAF, TDP, CDP, tumor_sample_barcode) |>
    dplyr::mutate(
      QUAL = ".",
      FILTER = 'PASS',
      ID = '.') |>
    dplyr::distinct()

  unique_samples <- unique(all_sample_calls$tumor_sample_barcode)
  i <- 1
  for(sample in unique_samples){
    sample_calls <- all_sample_calls |>
      dplyr::filter(tumor_sample_barcode == sample)

    if(i %% 50 == 0){
      cat(i, "samples processed\n")
    }
    if(NROW(sample_calls) > 0){
      tumor <- unique(sample_calls$tumor)

      if(!dir.exists(
        file.path(
          output_dir, tcga_release,
          "vcf","per_sample", paste0("TCGA-",tumor)))){
        dir.create(
          file.path(output_dir, tcga_release,
                    "vcf","per_sample", paste0("TCGA-",tumor)))
      }

      output_vcf_fname <- file.path(
        output_dir, tcga_release,
        "vcf","per_sample", paste0("TCGA-",tumor),
        paste0(sample,".grch38.vcf.gz")
      )

      if(!file.exists(output_vcf_fname)){

        #cat("balle\t",NROW(sample_calls),"\n")
        sample_calls <- sample_calls |>
          dplyr::mutate(INFO = paste0(
            "TDP=",TDP,
            ";TVAF=",TVAF,
            ";CDP=",CDP,
            ";ALGORITHMS=",ALGORITHMS,
            ";SAMPLE_TYPE=",SAMPLE_TYPE,
            ";TUMOR_SAMPLE_BARCODE=", tumor_sample_barcode)) |>
          dplyr::mutate(
            QUAL = ".",
            FILTER = 'PASS',
            ID = '.') |>
          dplyr::select(CHROM, POS, REF, ALT, QUAL, FILTER, ID, INFO)

        header_lines <- c(
          "##fileformat=VCFv4.2",
          "##INFO=<ID=TVAF,Number=1,Type=Float,Description=\"Fraction of reads supporting alternate allele - tumor\">",
          "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Sequencing read depth at variant site - tumor\">",
          "##INFO=<ID=TUMOR_SAMPLE_BARCODE,Number=1,Type=String,Description=\"TCGA tumor sample barcode\">",
          "##INFO=<ID=SAMPLE_TYPE,Number=1,Type=String,Description=\"Sample type\">",
          "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms supporting the variant call\">",
          "##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Sequencing read depth at variant site - control\">")

        ## Write to VCF
        vcfhelpR::write_vcf_records(
          vcf_records = sample_calls,
          output_dir =
            file.path(
              output_dir,
              tcga_release,
              "vcf", "per_sample", paste0("TCGA-",tumor)),
          vcf_fname_prefix = sample,
          header_lines = header_lines,
          genome_build = "grch38",
          fsep = "."
        )

        source_vcf_fname <-
          file.path(
            output_dir,
            tcga_release,
            "vcf", "per_sample", paste0("TCGA-",tumor),
            paste0(sample,".grch38.vcf.gz"))

        if(file.exists(source_vcf_fname)){
          ## CrossMap to grch37
          vcfhelpR::crossmap_vcf(
            source_vcf =
              source_vcf_fname,
            fsep = ".",
            target_vcf =
              file.path(
                output_dir,
                tcga_release,
                "vcf", "per_sample", paste0("TCGA-",tumor),
                paste0(sample,".grch37.vcf")
              ),
            direction = "hg38Tohg19"
          )
        }
      }
    }
    i <- i + 1
  }


}


#' Write TCGA VCF
#'
#' Function to write TCGA mutation calls to VCF format, including INFO fields for
#' sample barcodes, variant frequencies across TCGA cohorts, and flags for
#' overlaps with repeat regions.
#'
#' @param tcga_calls Data frame of TCGA mutation calls, as returned by
#'  `gdc_tcga_snv()`.
#' @param output_dir Directory to write the VCF files to.
#' @param gdc_release GDC release version
#'
#' @return None. VCF files are written to the specified output directory.
#'
write_tcga_vcf <- function(
    tcga_calls = NA,
    output_dir = NA,
    gdc_release = NA){


  if(!dir.exists(
    file.path(
      output_dir, gdc_release, "vcf"))){
    dir.create(
      file.path(output_dir, gdc_release, "vcf"))
  }

  header_lines <- c(
    "##fileformat=VCFv4.2",
    paste0("##SOURCE_TCGA=",gdc_release),
    "##INFO=<ID=TCGA_SAMPLES,Number=.,Type=String,Description=\"TCGA samples with mutation. Format: <sample_barcode>:<t_depth><t_ref_count><t_alt_count><n_depth>\">",
    "##INFO=<ID=TCGA_PANCANCER_COUNT,Number=1,Type=Integer,Description=\"Raw variant count across all tumor types\">",
    "##INFO=<ID=TCGA_FREQUENCY,Number=.,Type=String,Description=\"Frequency of variant across TCGA cancer subtypes. Format: subtype|percent affected|affected cases|total cases\">",
    "##INFO=<ID=WINMASKER_HIT,Number=0,Type=Flag,Description=\"Overlap with UCSC repeat track windowmaskerSdust\">",
    "##INFO=<ID=SIMPLEREPEATS_HIT,Number=0,Type=Flag,Description=\"Overlap with UCSC simpleRepeats track\">")

  tcga_samples_per_var <-
    as.data.frame(
      tcga_calls |>
        dplyr::mutate(
          tumor_sample_barcode_af_count =
            paste(
              tumor_sample_barcode, t_depth, t_ref_count,
              t_alt_count, n_depth, sep = ":"
            )
        ) |>
        dplyr::group_by(
          CHROM, POS, REF, ALT) |>
        dplyr::summarise(TCGA_SAMPLES = paste(
          tumor_sample_barcode_af_count,collapse=","),
          .groups = "drop") |>
        dplyr::ungroup()
    )

  tcga_vcf_records <- as.data.frame(
    tcga_calls |>
      dplyr::left_join(
        tcga_samples_per_var,
        by = c("CHROM","POS","REF","ALT")) |>
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
        ";", winmasker_hit,
        ";", simplerepeats_hit)) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";;", ";"
      )) |>
      dplyr::mutate(INFO = stringr::str_replace(
        INFO, ";{1,}$", ""
      )) |>
      dplyr::select(
        CHROM, POS, ID, REF,
        ALT, QUAL, FILTER, INFO) |>
      dplyr::distinct()
  )

  ## Write to VCF
  vcfhelpR::write_vcf_records(
    vcf_records = tcga_vcf_records,
    output_dir =
      file.path(
        output_dir,
        gdc_release,
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
        gdc_release,
        "vcf",
        "tcga_grch38.vcf.gz"
      ),
    target_vcf =
      file.path(
        output_dir,
        gdc_release,
        "vcf",
        "tcga_grch37.vcf"
      ),
  )

  write(c(header_lines[2],
          header_lines[4:5]),
        file = file.path(
          output_dir,
          gdc_release,
          "vcf",
          "tcga.vcfanno.vcf_info_tags.txt"
        ), sep = "\n")


}

#' Calculate TCGA Tumor Mutational Burden (TMB)
#'
#' Function to calculate tumor mutational burden (TMB) for TCGA samples
#' based on mutation calls from `gdc_tcga_snv()`.
#'
#' @param tcga_calls Data frame of TCGA mutation calls, as returned by
#' `gdc_tcga_snv()`.
#' @param t_depth_min Minimum tumor sequencing depth to consider a variant call.
#' Default is 30.
#' @param t_vaf_min Minimum tumor variant allele frequency (VAF) to consider
#' a variant call. Default is 0.05.
#' @param gdc_release GDC release version.
#' @param output_dir Directory to write output files to.
#' @param overwrite Logical flag indicating whether to overwrite existing TMB
#' results. Default is FALSE.
#' @param exome_target_size_mb Size of the exome target region in megabases (Mb).
#' Default is 34.0 Mb, approximating the coding content of GENCODE
#' protein-coding exons (grch38/v33).
#' @return Data frame with TMB measures for TCGA samples.
#'
#' @export
#'
calculate_sample_tmb <- function(
    tcga_calls = NA,
    t_depth_min = 30,
    t_vaf_min = 0.05,
    gdc_release = NA,
    output_dir = NA,
    overwrite = FALSE,
    exome_target_size_mb = 34.0){

  assertthat::assert_that(
    !is.na(gdc_release),
    msg = "gdc_release is missing"
  )

  if(!dir.exists(
    file.path(
      output_dir, gdc_release, "tmb"))){
    dir.create(
      file.path(
        output_dir, gdc_release, "tmb"))
  }

  tcga_tmb_rds <-
    file.path(
      output_dir,
      gdc_release,
      "tmb",
      "tcga_tmb.rds"
    )

  if(file.exists(tcga_tmb_rds) & overwrite == FALSE){
    tcga_tmb <- readRDS(file = tcga_tmb_rds)
    return(tcga_tmb)
  }

  assertable::assert_colnames(
    tcga_calls,
    c("bcr_patient_barcode",
      "tumor_sample_barcode",
      "t_depth",
      "t_alt_count",
      "CODING_STATUS",
      "GENOMIC_CHANGE",
      "EXONIC_STATUS"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  tcga_clinical <- gdc_tcga_clinical(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir) |>
    dplyr::select(
      bcr_patient_barcode,
      tumor_sample_barcode,
      tumor,
      sample_type,
      primary_site,
      primary_diagnosis,
      primary_diagnosis_very_simplified
    )

  ## Filter calls based on depth and VAF
  tcga_calls_filtered <- tcga_calls |>
    dplyr::mutate(
      t_vaf = as.numeric(t_alt_count) / as.numeric(t_depth)
    ) |>
    dplyr::filter(
      as.numeric(t_depth) >= t_depth_min &
        as.numeric(t_vaf) >= t_vaf_min
    )


  ## default target size region of protein-coding portion for TCGA WES samples
  ## Approximate solution (34Mb) based on coding content of GENCODE protein-coding exons (grch38/v33)
  ## Notably, Most TCGA samples are apparently sequenced with Agilent Custom V2 Exome Bait (Broad Institute)
  ## (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6169918/), not been able to locate BED file for this target
  #exome_target_size_mb <- as.numeric(34.0)

  ## Approaches for TMB calculation (proxy for neoantigen load in tumor)
  # 1. Chalmers et al. Genome Medicine, 2017 - all coding mutations (SNVs + InDels),
  #.   also including synonymous/silent mutations
  # 2. Fernandez et al. JCO Precision Oncology, 2019 - missense mutations only
  # 3. Chalmers except for silent mutations


  tcga_tmb_measures <- list()

  tcga_tmb_measures[['TMB_coding_and_silent']] <- as.data.frame(
    tcga_calls_filtered |>
      dplyr::filter(
        EXONIC_STATUS == "exonic") |>
      dplyr::select(GENOMIC_CHANGE,
                    bcr_patient_barcode,
                    tumor_sample_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode, tumor_sample_barcode) |>
      dplyr::summarise(n_coding_silent_muts = dplyr::n(),
                       .groups = "drop")
  )

  tcga_tmb_measures[['TMB_coding_non_silent']] <- as.data.frame(
    tcga_calls_filtered |>
      dplyr::filter(CODING_STATUS == "coding") |>
      dplyr::select(GENOMIC_CHANGE,
                    bcr_patient_barcode,
                    tumor_sample_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode,
                      tumor_sample_barcode) |>
      dplyr::summarise(n_coding_non_silent_muts = dplyr::n(),
                       .groups = "drop")
  )

  tcga_tmb_measures[['TMB_missense_only']] <- as.data.frame(
    tcga_calls_filtered |>
      dplyr::filter(
        startsWith(Consequence,"missense_variant")) |>
      dplyr::select(GENOMIC_CHANGE,
                    bcr_patient_barcode,
                    tumor_sample_barcode) |>
      dplyr::distinct() |>
      dplyr::group_by(bcr_patient_barcode,
                      tumor_sample_barcode) |>
      dplyr::summarise(n_missense_only_muts = dplyr::n(),
                       .groups = "drop")
  )

  tcga_tmb <- tcga_tmb_measures[['TMB_coding_and_silent']] |>
    dplyr::mutate(n_coding_silent_muts = dplyr::if_else(
      is.na(n_coding_silent_muts),
      as.integer(0),
      as.integer(n_coding_silent_muts))) |>
    dplyr::mutate(TMB_coding_and_silent = round(
      as.numeric(n_coding_silent_muts) /
        exome_target_size_mb,digits = 5)) |>
    dplyr::left_join(
      tcga_tmb_measures[['TMB_missense_only']],
      by = c("bcr_patient_barcode",
             "tumor_sample_barcode")) |>
    dplyr::mutate(n_missense_only_muts = dplyr::if_else(
      is.na(n_missense_only_muts),
      as.integer(0),
      as.integer(n_missense_only_muts))) |>
    dplyr::mutate(TMB_missense_only = round(
      as.numeric(n_missense_only_muts) /
        exome_target_size_mb, digits = 5)) |>
    dplyr::left_join(
      tcga_tmb_measures[['TMB_coding_non_silent']],
      by = c("bcr_patient_barcode",
             "tumor_sample_barcode")) |>
    dplyr::mutate(n_coding_non_silent_muts = dplyr::if_else(
      is.na(n_coding_non_silent_muts),
      as.integer(0),
      as.integer(n_coding_non_silent_muts))) |>
    dplyr::mutate(TMB_coding_non_silent = round(
      as.numeric(n_coding_non_silent_muts) /
        exome_target_size_mb, digits = 5)) |>
    dplyr::filter(!is.na(TMB_coding_and_silent)) |>
    dplyr::inner_join(
      tcga_clinical,
      by = c("bcr_patient_barcode",
             "tumor_sample_barcode")) |>
    dplyr::mutate(cds_target_size_mb = exome_target_size_mb)

  sample_stats <- tcga_tmb |>
    dplyr::group_by(primary_site, primary_diagnosis_very_simplified) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(for_subtype_tmb_comparison = dplyr::if_else(
      n > 15,
      TRUE,
      FALSE
    )) |>
    dplyr::select(-c("n"))

  tcga_tmb <- tcga_tmb |>
    dplyr::left_join(
      sample_stats,
      by = c("primary_site",
             "primary_diagnosis_very_simplified")
    ) |>
    dplyr::distinct()

  saveRDS(
    tcga_tmb,
    file = tcga_tmb_rds
  )

  readr::write_tsv(
    tcga_tmb,
    file = file.path(
      output_dir,
      gdc_release,
      "tmb",
      "tcga_tmb.tsv.gz"
    )
  )

  return(tcga_tmb)
}



# calculate_gene_mutation_rate <- function(
#     output_dir = NA,
#     tcga_calls = NA,
#     tcga_clinical_info = NA,
#     tcga_release = NA,
#     overwrite = FALSE){
#
#   if(!dir.exists(
#     file.path(
#       output_dir, gdc_release, "gene"))){
#     dir.create(
#       file.path(
#         output_dir, gdc_release, "gene"))
#   }
#
#   tcga_gene_mut_rate_rds <-
#     file.path(
#       output_dir,
#       gdc_release,
#       "gene",
#       "tcga_gene_aberration_rate.rds"
#     )
#
#   if(file.exists(tcga_gene_mut_rate_rds) &
#      overwrite == FALSE){
#     tcga_gene_mut_rate <-
#       readRDS(file = tcga_gene_mut_rate_rds)
#     return(tcga_gene_mut_rate)
#   }
#
#   cna_calls <- readRDS(
#     file = file.path(
#       output_dir,
#       tcga_release,
#       "cna",
#       "tcga_cna_gistic2.ampl_homdel.rds")
#   )
#
#   tcga_calls <- tcga_calls |>
#     dplyr::rename(symbol = Hugo_Symbol)
#
#   cohort_stats1 <- cohort_mutation_stats(
#     calls_df = cna_calls,
#     clinical_df = tcga_clinical,
#     vartype = "cna_homdel",
#     clinical_strata = "site") |>
#     dplyr::filter(!is.na(symbol))
#
#   cohort_stats2 <- cohort_mutation_stats(
#     calls_df = cna_calls,
#     clinical_df = tcga_clinical,
#     vartype = "cna_homdel",
#     clinical_strata = "site_diagnosis") |>
#     dplyr::filter(!is.na(symbol))
#
#   cohort_stats3 <- cohort_mutation_stats(
#     calls_df = cna_calls,
#     clinical_df = tcga_clinical,
#     vartype = "cna_ampl",
#     clinical_strata = "site") |>
#     dplyr::filter(!is.na(symbol))
#
#   cohort_stats4 <- cohort_mutation_stats(
#     calls_df = cna_calls,
#     clinical_df = tcga_clinical,
#     vartype = "cna_ampl",
#     clinical_strata = "site_diagnosis") |>
#     dplyr::filter(!is.na(symbol))
#
#   cohort_stats5 <- cohort_mutation_stats(
#     calls_df = tcga_calls,
#     clinical_df = tcga_clinical,
#     vartype = "snv_indel",
#     clinical_strata = "site",
#     genomic_strata = "gene") |>
#     dplyr::filter(!is.na(symbol))
#
#   cohort_stats7 <- cohort_mutation_stats(
#     calls_df = tcga_calls,
#     clinical_df = tcga_clinical,
#     vartype = "snv_indel",
#     clinical_strata = "site_diagnosis",
#     genomic_strata = "gene") |>
#     dplyr::filter(!is.na(symbol))
#
#
#   tcga_gene_stats <- dplyr::bind_rows(
#     cohort_stats7, cohort_stats5,
#     cohort_stats4, cohort_stats3,
#     cohort_stats2, cohort_stats1)
#
#   saveRDS(
#     tcga_gene_stats,
#     file = tcga_gene_mut_rate_rds
#   )
#
#   return(tcga_gene_stats)
# }
#

#' Calculate TCGA Gene Mutation Rate
#'
#' Function to calculate gene-level mutation and copy number alteration
#' frequencies across TCGA cohorts.
#'
#' @param output_dir Directory to write output files to.
#' @param gdc_release GDC release version.
#' @param gdc_projects Vector of GDC project identifiers (e.g., "TC
#' GA-BRCA").
#' @param overwrite Logical flag indicating whether to overwrite existing
#' results. Default is FALSE.
#'
#' @return Data frame with gene-level mutation and CNA frequencies across TCGA
#' cohorts.
#' @export
#'
calculate_gene_mutation_rate <- function(
    output_dir = NULL,
    gdc_release = NULL,
    gdc_projects = NULL,
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(gdc_release),
    !is.null(gdc_projects),
    !is.null(output_dir),
    msg = "gdc_release or gdc_projects or output_dir is missing"
  )

  if(!dir.exists(
    file.path(
      output_dir, gdc_release, "gene"))){
    dir.create(
      file.path(
        output_dir, gdc_release, "gene"))
  }

  tcga_gene_mut_rate_rds <-
    file.path(
      output_dir,
      gdc_release,
      "gene",
      "tcga_gene_aberration_rate.rds"
    )

  if(file.exists(tcga_gene_mut_rate_rds) &
     overwrite == FALSE){
    tcga_gene_mut_rate <-
      readRDS(file = tcga_gene_mut_rate_rds)
    return(tcga_gene_mut_rate)
  }

  tcga_clinical <- gdc_tcga_clinical(
    gdc_release = gdc_release,
    overwrite = F,
    output_dir = output_dir) |>
    dplyr::select(
      bcr_patient_barcode,
      tumor_sample_barcode,
      tumor,
      sample_type,
      primary_site,
      primary_diagnosis,
      primary_diagnosis_very_simplified
    )

  cna_calls <- data.frame()
  for(project in gdc_projects){
    project2 <- stringr::str_replace(
      project, "TCGA-",""
    )
    cna_calls_project <- readRDS(
      file = file.path(
        output_dir,
        gdc_release,
        "cna",
        glue::glue(
          "tcga_cna_ASCAT3_{project2}.rds"))) |>
      dplyr::filter(
        mut_status %in% c("AMPL","HOMDEL")
      )
    cna_calls <- dplyr::bind_rows(
      cna_calls,
      cna_calls_project
    )

  }

  cna_calls <- cna_calls |>
    dplyr::left_join(
      tcga_clinical, by = "tumor_sample_barcode")

  snv_indel_calls <- gdc_tcga_snv(
    gdc_release = gdc_release,
    gdc_projects = gdc_projects,
    overwrite = F,
    output_dir = output_dir) |>
    dplyr::mutate(
      symbol = Hugo_Symbol
    )

  cohort_stats <- list()
  cohort_stats[['snv_indel']] <- list()
  cohort_stats[['cna_ampl']] <- list()
  cohort_stats[['cna_homdel']] <- list()
  for(t in c('snv_indel','cna_ampl','cna_homdel')){
    for(cs in c('site','site_diagnosis')){
      cohort_stats[[t]][[cs]] <- data.frame()
    }

  }

  cohort_stats[['cna_homdel']][['site']] <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical,
    vartype = "cna_homdel",
    clinical_strata = "site") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats[['cna_homdel']][['site_diagnosis']] <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical,
    vartype = "cna_homdel",
    clinical_strata = "site_diagnosis") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats[['cna_ampl']][['site']] <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical,
    vartype = "cna_ampl",
    clinical_strata = "site") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats[['cna_ampl']][['site_diagnosis']] <- cohort_mutation_stats(
    calls_df = cna_calls,
    clinical_df = tcga_clinical,
    vartype = "cna_ampl",
    clinical_strata = "site_diagnosis") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats[['snv_indel']][['site']] <- cohort_mutation_stats(
    calls_df = snv_indel_calls,
    clinical_df = tcga_clinical,
    vartype = "snv_indel",
    clinical_strata = "site",
    genomic_strata = "gene") |>
    dplyr::filter(!is.na(symbol))

  cohort_stats[['snv_indel']][['site_diagnosis']] <- cohort_mutation_stats(
    calls_df = snv_indel_calls,
    clinical_df = tcga_clinical,
    vartype = "snv_indel",
    clinical_strata = "site_diagnosis",
    genomic_strata = "gene") |>
    dplyr::filter(!is.na(symbol))


  tcga_gene_stats <- dplyr::bind_rows(
    cohort_stats[['snv_indel']][['site_diagnosis']],
    cohort_stats[['snv_indel']][['site']],
    cohort_stats[['cna_ampl']][['site_diagnosis']],
    cohort_stats[['cna_ampl']][['site']],
    cohort_stats[['cna_homdel']][['site_diagnosis']],
    cohort_stats[['cna_homdel']][['site']]
  )

  saveRDS(
    tcga_gene_stats,
    file = tcga_gene_mut_rate_rds
  )

  return(tcga_gene_stats)
}
