
get_bulk_rnaseq <- function(
    gdc_projects = NA,
    tcga_clinical_info = NA,
    tcga_release = NA,
    data_raw_dir = NA,
    overwrite = FALSE,
    gOncoX = NA,
    output_dir = NA){


  if(!dir.exists(
    file.path(
      output_dir, tcga_release, "rnaseq"))){
    dir.create(
      file.path(output_dir, tcga_release, "rnaseq"))
  }

  ## NOTE: As of GDC release v32, GENCODE v36 is used for annotation of RNA-seq genes/transcripts
  ## https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-320
  update_ensembl_gene <- F
  gencode_xref_fname <-
    file.path(data_raw_dir, "gencode_v36_protein_coding_xref.rds")
  gencode_xref <- NULL
  if(!file.exists(gencode_xref_fname) | update_ensembl_gene == T){

    gencode_gtf_fname <- file.path(
      data_raw_dir,
      "gencode",
      "gencode.v36.annotation.gtf.gz")

    gencode_xref <- as.data.frame(
      valr::read_gtf(path = gencode_gtf_fname) |>
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
      ## remove ensembl gene id's that are pointing to the same gene symbol as other
      ## ensembl identifiers
      dplyr::filter(ensembl_gene_id != "ENSG00000145075") |>
      dplyr::filter(ensembl_gene_id != "ENSG00000285437") |>
      dplyr::filter(ensembl_gene_id != "ENSG00000286065") |>
      dplyr::distinct() |>
      dplyr::inner_join(
        dplyr::select(gOncoX$basic$records,
                      symbol, hgnc_id, gene_biotype),
        by = "hgnc_id"
      ) |>
      dplyr::filter(
        gene_biotype == "protein-coding")
    )

    saveRDS(gencode_xref, file = gencode_xref_fname)
  }else{
    gencode_xref <- readRDS(file = gencode_xref_fname)
  }

  for(i in 1:nrow(gdc_projects)){

    rnaseq_query <- NULL
    rnaseq_data <- list()

    project_code <- gdc_projects[i,]$project_id
    tumor <- gdc_projects[i,]$tumor

    cat(i, tumor, sep=" - ")
    cat('\n')

    downloads_raw_gdc <- file.path(
      data_raw_dir, "GDCdata",
      project_code, "harmonized",
      "Transcriptome_Profiling",
      "Gene_Expression_Quantification"
    )

    rds_fname <-
      file.path(
        output_dir, tcga_release, "rnaseq",
        paste0('rnaseq_',tumor,'.rds')
      )

    if(file.exists(rds_fname)){
      next
    }

    rnaseq_data <- list()
    workflow_type <- "STAR - Counts"
    rnaseq_query <-
      TCGAbiolinks::GDCquery(
        project = project_code,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = workflow_type
      )

    TCGAbiolinks::GDCdownload(
      rnaseq_query,
      directory = file.path(
        data_raw_dir, "GDCdata"),
      files.per.chunk = 100)

    rnaseq_data[['tcga_release']] <- tcga_release
    rnaseq_data[['tcga_project']] <- tumor
    rnaseq_data[['sample_metadata']] <- data.frame()
    rnaseq_data[['tcga_biolinks']] <- list()
    rnaseq_data[['tcga_biolinks']][['data.category']] <-
      'Transcriptome Profiling'
    rnaseq_data[['tcga_biolinks']][['data.type']] <-
      'Gene Expression Quantification'
    rnaseq_data[['tcga_biolinks']][['workflow.type']] <-
      workflow_type
    rnaseq_data[['tcga_biolinks']][['version']] <-
      Biobase::package.version("TCGAbiolinks")
    rnaseq_raw_SE <-
      TCGAbiolinks::GDCprepare(
        rnaseq_query, directory =
          file.path(data_raw_dir, "GDCdata"))


    system(paste0('rm -rf ', downloads_raw_gdc, "/*"))
    rnaseq_data[['matrix']] <- list()
    rnaseq_data[['matrix']][['counts']] <- list()
    rnaseq_data[['matrix']][['tpm']] <- list()

    ## Get counts and TPM data matrices from SummarizedExperiment object
    rnaseq_data[['matrix']][['tpm']][['all']] <-
      SummarizedExperiment::assay(rnaseq_raw_SE, "tpm_unstrand")

    rnaseq_data[['matrix']][['counts']][['all']] <-
      SummarizedExperiment::assay(rnaseq_raw_SE)


    ## Strip away barcode information related to Portion/Plate/Center
    ## https://www.biostars.org/p/313063/
    ## Example: TCGA-93-A4JP-01A-11H-A24S-13 --> TCGA-93-A4JP-01A
    ##          "Barcode" --> "Sample"
    ##
    ## Strip away Ensembl gene version, e.g.:
    ## ENSG00000000003.15 --> ENSG00000000003
    for(c in c('tpm', 'counts')){
      colnames(rnaseq_data[['matrix']][[c]][['all']]) <-
        stringr::str_replace(colnames(rnaseq_data[['matrix']][[c]][['all']]),
                             "-[0-9]{2}[A-Z]{1}-[A-Z|0-9]{4}-\\S+$",
                             "")
      rownames(rnaseq_data[['matrix']][[c]][['all']]) <-
        stringr::str_replace(rownames(rnaseq_data[['matrix']][[c]][['all']]),
                             "\\.[0-9]{1,}$",
                             "")

    }

    rnaseq_data[['sample_metadata']] <-
      as.data.frame(SummarizedExperiment::colData(rnaseq_raw_SE)) |>
      dplyr::select(barcode, bcr_patient_barcode, sample_type,
                    name, primary_diagnosis,
                    site_of_resection_or_biopsy)

    rownames(rnaseq_data[['sample_metadata']]) <- NULL

    sample_collection <- list()
    sample_collection[['tumor']] <- rnaseq_data$sample_metadata |>
      dplyr::filter(
        stringr::str_detect(.data$bcr_patient_barcode, "-(0[1-9][A-Z])$"))
    sample_collection[['control']] <- rnaseq_data$sample_metadata |>
      dplyr::filter(
        stringr::str_detect(.data$bcr_patient_barcode, "-(1[0-4][A-Z])$"))

    ## Assign data from tumor and control/normal samples in dedicated TPM/count matrices
    for(exp_measure in c('tpm','counts')){
      for(c in c('tumor', 'control')){

        ## Filter complete sample set against tumor/control sample set
        rnaseq_data[['matrix']][[exp_measure]][[c]] <- as.matrix(
          rnaseq_data[['matrix']][[exp_measure]][['all']][, colnames(rnaseq_data[['matrix']][[exp_measure]][['all']]) %in%
                                                            sample_collection[[c]]$bcr_patient_barcode]
        )


        if(!is.null(rnaseq_data[['matrix']][[exp_measure]][[c]])){
          if(NCOL(rnaseq_data[['matrix']][[exp_measure]][[c]]) > 1){

            ## Exclude duplicated sample names
            rnaseq_data[['matrix']][[exp_measure]][[c]] <-
              rnaseq_data[['matrix']][[exp_measure]][[c]][, !duplicated(
                colnames(rnaseq_data[['matrix']][[exp_measure]][[c]]), fromLast = TRUE)]

            df <- data.frame(
              'ensembl_gene_id' = rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]),
              stringsAsFactors = F) |>
              dplyr::left_join(gencode_xref, by = "ensembl_gene_id")


            rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]) <- df$symbol
            rnaseq_data[['matrix']][[exp_measure]][[c]] <-
              as.matrix(rnaseq_data[['matrix']][[exp_measure]][[c]][!is.na(rownames(rnaseq_data[['matrix']][[exp_measure]][[c]])), ])
            rnaseq_data[['matrix']][[exp_measure]][[c]] <-
              as.matrix(rnaseq_data[['matrix']][[exp_measure]][[c]][rownames(rnaseq_data[['matrix']][[exp_measure]][[c]]) != "NA", ])
          }
        }
      }

    }



    ## remove matrices with all samples (tumor + normal created)
    rnaseq_data$matrix$counts$all <- NULL
    rnaseq_data$matrix$tpm$all <- NULL

    saveRDS(rnaseq_data, rds_fname)

  }
}
