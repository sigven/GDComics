
## GISTIC CNA ENCODING
## -2: homozygous deletion (HOMDEL)
## -1: hemizygous deletion (HEMDEL)
##  0: neutral/no change (NONE)
##  1: gain (GAIN)
##  2: high-level amplification (AMPL)

get_cna_calls <- function(
  gdc_projects = NA,
  tcga_clinical_info = NA,
  tcga_release = NA,
  data_raw_dir = NA,
  gOncoX = NA,
  overwrite = F,
  clear_cache = T,
  output_dir = NA){

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
            #dplyr::filter(mut_status == "HOMDEL" | mut_status == "AMPL") |>
            dplyr::mutate(tumor = t) |>
            dplyr::filter(!is.na(entrezgene)) |>
            dplyr::distinct()
        }
        cna_df[[e]] <- cna_calls
      }

      tmp <- cna_df$all |>
        dplyr::select(tumor_sample_barcode,
                      cna_code, entrezgene, bcr_patient_barcode) |>
        dplyr::rename(cna_signal_raw = cna_code)
      #cna_calls <- cna_df

      cna_calls_final <- cna_df$thresholded |>
        dplyr::inner_join(
          tmp,
          by = c("entrezgene",
                 "tumor_sample_barcode",
                 "bcr_patient_barcode")) |>
        dplyr::left_join(
          dplyr::select(tcga_clinical_info$slim, bcr_patient_barcode,
                        primary_site, primary_diagnosis,
                        primary_diagnosis_very_simplified),
          by = "bcr_patient_barcode") |>
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

}


# ### RAW SEGMENTS
#
# all_cna_segments <- data.frame()
# options(scipen = 999)
#
# for(project_id in sort(gdc_projects$project_id)){
#   if(!stringr::str_detect(project_id,"TCGA-")){
#     next
#   }
#   cat(project_id,sep="\n")
#
#   #library(TCGAbiolinks)
#   query_nocnv_primary <- TCGAbiolinks::GDCquery(project = project_id,
#                               data.category = "Copy number variation",
#                               legacy = TRUE,
#                               file.type = "nocnv_hg19.seg",
#                               sample.type = c("Primary Tumor"))
#
#   if(!is.null(query_nocnv_primary)){
#     TCGAbiolinks::GDCdownload(query_nocnv_primary, directory = "data-raw/GDCdata")
#     data <- as.data.frame(
#       TCGAbiolinks::GDCprepare(
#         query_nocnv_primary,
#         directory = "data-raw/GDCdata", save = F,
#         save.filename = paste0("output/cna/segments_SNPArray_",
#                                stringr::str_replace(project_id,"TCGA-",""),
#                                ".", tcga_release,".TP.grch37.rds"))
#     ) |>
#       dplyr::mutate(Tumor_Sample_Barcode = stringr::str_replace(
#         Sample,"-[0-9][0-9][A-Z]-[0-9A-Z]{4}-[0-9]{2}$","")) |>
#       dplyr::mutate(name = paste0(
#         paste(Chromosome,Start,End,sep="_"),":",Tumor_Sample_Barcode,":TP:",Num_Probes,":grch37")) |>
#       dplyr::mutate(score = Segment_Mean) |>
#       dplyr::rename(chrom = Chromosome, chromStart = Start, chromEnd = End) |>
#       dplyr::select(chrom, chromStart, chromEnd, name, score) |>
#       dplyr::distinct()
#     all_cna_segments <- all_cna_segments |>
#       dplyr::bind_rows(data)
#   }
#   query_nocnv_metastatic <- TCGAbiolinks::GDCquery(
#     project = project_id,
#     data.category = "Copy number variation",
#     legacy = TRUE,
#     file.type = "nocnv_hg19.seg",
#     sample.type = c("Metastatic"))
#   if(!is.null(query_nocnv_metastatic)){
#     TCGAbiolinks::GDCdownload(query_nocnv_metastatic, directory = "data-raw/GDCdata")
#     data <- as.data.frame(
#       TCGAbiolinks::GDCprepare(
#         query_nocnv_metastatic,
#         directory = "data-raw/GDCdata", save = F,
#         save.filename = paste0("output/cna/segments_SNPArray_",
#                                stringr::str_replace(project_id,"TCGA-",""),
#                                ".", tcga_release,".TM.grch37.rds"))
#     ) |>
#       dplyr::mutate(Tumor_Sample_Barcode = stringr::str_replace(
#         Sample,"-[0-9][0-9][A-Z]-[0-9A-Z]{4}-[0-9]{2}$","")) |>
#       dplyr::mutate(name = paste0(
#         paste(Chromosome,Start,End,sep="_"),":",Tumor_Sample_Barcode,":TM:",
#         Num_Probes,":grch37")) |>
#       dplyr::mutate(score = Segment_Mean) |>
#       dplyr::rename(chrom = Chromosome, chromStart = Start, chromEnd = End) |>
#       dplyr::select(chrom, chromStart, chromEnd, name, score) |>
#       dplyr::distinct()
#     all_cna_segments <- all_cna_segments |>
#       dplyr::bind_rows(data)
#   }
# }
#
# chrOrder <- c(as.character(c(1:22)),"X","Y")
# all_cna_segments$chrom <- factor(all_cna_segments$chrom, levels=chrOrder)
# all_cna_segments <- all_cna_segments[order(all_cna_segments$chrom),]
# all_cna_segments$chromStart <- as.numeric(all_cna_segments$chromStart)
# all_cna_segments$chromEnd <- as.numeric(all_cna_segments$chromEnd)
# all_cna_segments <- all_cna_segments |>
#   dplyr::filter(chromStart < chromEnd)
#
# all_cna_segments_sorted <- data.frame()
# for(chrom in chrOrder){
#   if(nrow(all_cna_segments[all_cna_segments$chrom == chrom,]) > 0){
#       chrom_regions <- all_cna_segments[all_cna_segments$chrom == chrom,]
#       chrom_cna_segments_sorted <- chrom_regions[with(chrom_regions, order(chromStart, chromEnd)),] |>
#         dplyr::select(chrom, chromStart, chromEnd, name, score) |>
#         dplyr::filter(!is.na(chromStart)) |>
#         dplyr::mutate(chromStart = as.integer(chromStart)) |>
#         dplyr::mutate(chromEnd = as.integer(chromEnd))
#
#       write.table(chrom_cna_segments_sorted,
#                   file = paste0("output/cna/all_cna_segments.",
#                                 chrom,".", tcga_release,".grch37.bed"),
#                   col.names = F, row.names = F, quote = F, sep ="\t")
#       # crossmapr::crossmap_bed(direction = "hg19Tohg38",
#       #                         source_bed = paste0("output/cna/all_cna_segments.", chrom,".", tcga_release,".grch37.bed"),
#       #                         target_bed = paste0("output/cna/all_cna_segments.", chrom,".", tcga_release,".grch38.bed"))
#
#       all_cna_segments_sorted <- all_cna_segments_sorted |>
#         dplyr::bind_rows(chrom_cna_segments_sorted)
#
#   }
#   cat(chrom,'\n')
# }
#
#
#
# write.table(all_cna_segments_sorted, file = paste0("output/cna/all_cna_segments.", chrom,".", tcga_release,".tsv"),
#             col.names = T, row.names = F, quote = F, sep ="\t")
# all_cna_segments_sorted$chrom <- paste0('chr',all_cna_segments_sorted$chrom)
# write.table(all_cna_segments_sorted, file = paste0("output/cna/all_cna_segments.", chrom,".", tcga_release,".grch37.bed"),
#             col.names = F, row.names = F, quote = F, sep ="\t")
# #system(paste0('cat '))


















