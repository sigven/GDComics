#library(TCGAbiolinks)
#library(BSgenome.Hsapiens.UCSC.hg38)

#source("code/data_preprocessing/data_cleaning_utils.R")

get_tcga_clinical <- function(
    gdc_projects = NA,
    tcga_release = NA,
    output_dir = NA,
    overwrite = FALSE){

  rds_fname = file.path(
    output_dir,
    tcga_release,
    "clinical",
    "tcga_clinical.rds"
  )

  if(file.exists(rds_fname) & overwrite == F){
    tcga_clinical_info <- readRDS(file = rds_fname)
    return(tcga_clinical_info)
  }

  tcga_clinical <- data.frame()

  for (i in 1:nrow(gdc_projects)) {

    project_code <- gdc_projects[i,]$project_id
    tumor <- gdc_projects[i,]$tumor
    cat(i,tumor,project_code,sep = "-")
    cat('\n')

    ## get basic clinical information pr sample for cohort from GDC
    clinical <- TCGAbiolinks::GDCquery_clinic(
      project = project_code, type = "clinical")
    clinical <- clinical |>
      dplyr::select(
        -c(submitter_id,classification_of_tumor, exposure_id,demographic_id,
                       treatments_pharmaceutical_treatment_id, diagnosis_id)) |>
      #plyr::filter(!is.na(updated_datetime)) |>
      dplyr::mutate(tumor = tumor) |>
      dplyr::mutate(project = project_code) |>
      dplyr::distinct()

    if ('NA.' %in% colnames(clinical)) {
      clinical <- clinical |>
        dplyr::select(-c(NA.))
    }


    ## For BRCA and PRAD cohorts, get extra clinical annotations
    ## 1. ER status (BRCA)
    ## 2. PR status (BRCA)
    ## 3. HER2 status (BRCA)
    ## 4. Gleason score (PRAD)
    if (project_code == "TCGA-BRCA" | project_code == "TCGA-PRAD") {
      query <- TCGAbiolinks::GDCquery(project = project_code,
                        data.category = "Clinical",
                        data.type = "Clinical Supplement",
                        data.format = "BCR Biotab")
      TCGAbiolinks::GDCdownload(query)
      clinical_supplement_bcr_biotab <- TCGAbiolinks::GDCprepare(query)

      if (project_code == "TCGA-BRCA") {
        table <- "clinical_patient_brca"
        brca_er_pr_her2_status <- as.data.frame(
          dplyr::select(clinical_supplement_bcr_biotab[[table]],
                                         bcr_patient_barcode,
                                         er_status_by_ihc,
                                         her2_status_by_ihc,
                                         pr_status_by_ihc) |>
          dplyr::filter(startsWith(bcr_patient_barcode,"TCGA-")) |>
            dplyr::rename(ER_status = er_status_by_ihc,
                          HER2_status = her2_status_by_ihc,
                          PR_status = pr_status_by_ihc)

        )
        clinical <- clinical |> dplyr::left_join(brca_er_pr_her2_status)
      }else{
        clinical$ER_status <- as.character(NA)
        clinical$PR_status <- as.character(NA)
        clinical$HER_status <- as.character(NA)
      }

      if (project_code == "TCGA-PRAD") {
        table <- 'clinical_patient_prad'
        prad_gleason <- dplyr::select(
          clinical_supplement_bcr_biotab[[table]],
          bcr_patient_barcode,
          gleason_score) |>
          dplyr::filter(startsWith(bcr_patient_barcode,"TCGA-")) |>
          dplyr::mutate(Gleason_score = as.integer(gleason_score))
        clinical <- clinical |> dplyr::left_join(prad_gleason)

      }else{
        clinical$Gleason_score <- as.integer(NA)
      }
    }

    ## For STAD/READ/COAD/UCEC, get MSI status
    if (project_code == 'TCGA-STAD' | project_code == 'TCGA-READ' |
       project_code == 'TCGA-COAD' | project_code == 'TCGA-UCEC') {
        query <- TCGAbiolinks::GDCquery(
          project = project_code,
          data.category = "Other",
          legacy = TRUE,
          access = "open")
        TCGAbiolinks::GDCdownload(query)
        msi_assay_results <- TCGAbiolinks::GDCprepare_clinic(query, "msi")
        msi_assay_results <- msi_assay_results |>
          dplyr::select(
            bcr_patient_barcode,
            mononucleotide_and_dinucleotide_marker_panel_analysis_status) |>
          dplyr::rename(
            MSI_status = mononucleotide_and_dinucleotide_marker_panel_analysis_status) |>
          dplyr::mutate(MSI_status = as.character(MSI_status)) |>
          dplyr::filter(MSI_status != 'Indeterminate' &
                          !is.na(MSI_status)) |>
          dplyr::distinct() |>
          dplyr::mutate(MSI_status = dplyr::if_else(
            !is.na(MSI_status) & MSI_status == 'MSI-H',
            'MSI.H',as.character(MSI_status)
          )) |>
          dplyr::mutate(MSI_status = dplyr::if_else(
            !is.na(MSI_status) & MSI_status == 'MSI-L',
            'MSI.L',as.character(MSI_status)
          ))

      clinical <- clinical |> dplyr::left_join(msi_assay_results)
    }else{
      clinical$MSI_status <- as.character(NA)
    }

    tcga_clinical <- dplyr::bind_rows(tcga_clinical,clinical)

    i <- i + 1
  }

  ## For all cohorts, get PanCancer subtype information
  tcga_pancan_subtypes <- as.data.frame(TCGAbiolinks::PanCancerAtlas_subtypes()) |>
    janitor::clean_names() |>
    dplyr::rename(bcr_patient_barcode = pan_samples_id,
                  tumor = cancer_type,
                  pancan_subtype_mrna = subtype_m_rna,
                  pancan_subtype_cna = subtype_cna,
                  pancan_subtype_dnameth = subtype_dn_ameth,
                  pancan_subtype_protein = subtype_protein,
                  pancan_subtype_other = subtype_other,
                  pancan_subtype_integrative = subtype_integrative,
                  pancan_subtype_mirna = subtype_mi_rna,
                  pancan_subtype_selected = subtype_selected)

  ## Join uniform clinical annotations with
  ## subtype information, set primary tumor site
  tcga_clinical_final <- tcga_clinical |>
    dplyr::left_join(tcga_pancan_subtypes) |>
    dplyr::mutate(primary_site = "Other/Unknown") |>
    ## Simplify primary sites/tissue
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        primary_diagnosis,
        "(P|p)araganglioma") |
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "adrenal gland"),
      "Adrenal Gland",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "BLCA" & stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "bladder|ureteric"),
      "Bladder/Urinary Tract",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(primary_diagnosis),
        "cholangiocarcinoma") |
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "bile duct|gallbladder"),
      "Biliary Tract",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "bone|bones") &
        !stringr::str_detect(
          tolower(tissue_or_organ_of_origin),"marrow"),
      "Bone",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "bone marrow"),
      "Myeloid",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "brain|cerebrum| lobe"),
      "CNS/Brain",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "breast"),
      "Breast",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "cervix uteri"),
      "Cervix",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "colon|rectum|cecum|rectosigmoid"),
      "Colon/Rectum",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "gastric|stomach|cardia|pylorus|esophagus"),
      "Esophagus/Stomach",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tumor),"sarc|dlbc|hnsc") &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "oral cavity|mouth|larynx|tongue|pharynx|lip|palate|gum|tonsil|cheek|mandible|supraglottis"),"Head and Neck",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "kidney"),"Kidney",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor != "CHOL" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "liver"),"Liver",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "lung|bronchus"),"Lung",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "DLBC" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "lymph nodes"),"Lymphoid",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "OV" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "ovary"),"Ovary/Fallopian Tube",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "PAAD" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "pancreas"),"Pancreas",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "PRAD" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "prostate"),"Prostate",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "SKCM" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "skin"),"Skin",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "thyroid gland"),"Thyroid",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "uterus, nos| uteri|endometrium") &
        !stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "cervix"),"Uterus",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor != "PCPG" & stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "connective"),"Soft Tissue",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "THYM" & stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "thymus|mediastinum"),"Thymus",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(tissue_or_organ_of_origin),
        "testis, nos"),"Testis",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor != "THYM" & tumor != "PCPG" &
        stringr::str_detect(
          tolower(tissue_or_organ_of_origin),
          "mediastinum|pleura"),"Pleura",primary_site)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      tumor == "UVM","Eye",primary_site)) |>
    dplyr::mutate(primary_diagnosis = dplyr::if_else(
      tumor == "READ" & (primary_diagnosis == "Adenocarcinoma, NOS" |
                           primary_diagnosis == "Mucinous adenocarcinoma"),
      paste0("Rectal ",primary_diagnosis),primary_diagnosis)) |>
    dplyr::mutate(primary_diagnosis = dplyr::if_else(
      tumor == "COAD" & (primary_diagnosis == "Adenocarcinoma, NOS" |
                           primary_diagnosis == "Mucinous adenocarcinoma"),
      paste0("Colon ",primary_diagnosis),primary_diagnosis)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        primary_diagnosis,"sarcoma|nerve sheath|Malignant fibrous|Aggressive fibromatosis"),
      "Soft Tissue",as.character(primary_site)))

  ## clean survival data
  tcga_clinical_final <- tcga_clinical_final |>
    dplyr::mutate(days_to_death = dplyr::if_else(stringr::str_detect(
      tolower(vital_status),"alive"),-Inf,as.numeric(days_to_death))) |>
    dplyr::mutate(days_to_last_follow_up = dplyr::if_else(
      stringr::str_detect(tolower(vital_status),"dead"),-Inf,
      as.numeric(days_to_last_follow_up))) |>
    dplyr::mutate(death_event = dplyr::if_else(stringr::str_detect(
      tolower(vital_status),"alive"),as.integer(0),as.integer(1))) |>
    dplyr::mutate(new_death = dplyr::if_else(
      days_to_death == -Inf,days_to_last_follow_up,days_to_death))

  tcga_clinical_counts <- as.data.frame(tcga_clinical_final |>
    dplyr::group_by(tumor,primary_diagnosis,primary_site) |>
    dplyr::summarise(n_samples = dplyr::n()))

  ## Clean small subtypes
  tcga_clinical_counts <- tcga_clinical_counts |>
    dplyr::mutate(
      primary_diagnosis_simplified = dplyr::if_else(
        n_samples <= 5,"Other subtype(s)",as.character(primary_diagnosis))) |>
    dplyr::mutate(
      primary_diagnosis_simplified = dplyr::if_else(
        n_samples <= 5 & stringr::str_detect(
          primary_diagnosis,"^Thymoma"),"Other thymoma subtype(s)",
        as.character(primary_diagnosis_simplified))) |>
    dplyr::mutate(
      primary_diagnosis_simplified = dplyr::if_else(
        n_samples <= 5 & stringr::str_detect(
          primary_diagnosis,"melanoma"),
        "Other melanoma subtype(s)",
        as.character(primary_diagnosis_simplified))) |>
    dplyr::mutate(
      primary_diagnosis_simplified = dplyr::if_else(
        n_samples <= 5 & stringr::str_detect(
          primary_diagnosis,
          "sarcoma|nerve sheath|Malignant fibrous|Aggressive fibromatosis"),
        "Other sarcoma subtype(s)",
        as.character(primary_diagnosis_simplified))) |>
    ### Remove cases in which main tissue/origin of tumor type deviates
    dplyr::mutate(
      primary_site = dplyr::if_else(
        tumor == "DLBC" & primary_site != "Lymphoid",
        as.character(NA),primary_site)) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        tumor == "ESCA" & primary_site != "Esophagus/Stomach",
        as.character(NA),primary_site)) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        tumor == "HNSC" & primary_site != "Head and Neck",
        as.character(NA),primary_site)) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        tumor == "READ" & primary_site != "Colon/Rectum",
        as.character(NA),primary_site)) |>
    dplyr::mutate(
      primary_diagnosis_simplified = dplyr::if_else(
        is.na(primary_site),as.character(NA),
        primary_diagnosis_simplified)) |>
    dplyr::filter(!is.na(primary_diagnosis_simplified)) |>
    ## Clean bigger subtypes
    dplyr::mutate(
      primary_diagnosis_very_simplified = dplyr::if_else(
        n_samples <= 20,"Other subtype(s)",
        as.character(primary_diagnosis_simplified))) |>
    ### Remove cases in which main tissue/origin of tumor type deviates
    dplyr::mutate(primary_diagnosis_very_simplified = dplyr::if_else(
      is.na(primary_site),
      as.character(NA),
      primary_diagnosis_very_simplified)) |>
    dplyr::filter(!is.na(primary_diagnosis_very_simplified)) |>
    dplyr::mutate(tissue_code = paste0(
      toupper(substr(primary_site,0,2)),"_")) |>
    dplyr::mutate(tissue_code = dplyr::if_else(
      primary_site == "Breast","BRCA_",as.character(tissue_code))) |>
    dplyr::mutate(tissue_code = dplyr::if_else(
      primary_site == "Thymus","THYM_",as.character(tissue_code))) |>
    dplyr::mutate(tissue_code = dplyr::if_else(
      primary_site == "Thyroid","THYR_",as.character(tissue_code)))


  tissue_subtype_codes <- as.data.frame(
    tcga_clinical_counts |>
    dplyr::group_by(tissue_code, primary_site,primary_diagnosis_very_simplified) |>
      dplyr::summarise(n2 = sum(n_samples)) |>
      dplyr::arrange(tissue_code,desc(n2)) |>
      dplyr::ungroup() |>
      dplyr::group_by(tissue_code) |>
      dplyr::mutate(tissue_code_2 = paste0(tissue_code,dplyr::row_number()))
    ) |>
    dplyr::select(-c(n2,tissue_code))

  tcga_clinical_complete <- dplyr::select(tcga_clinical_final, -primary_site) |>
    dplyr::left_join(
      dplyr::select(
        tcga_clinical_counts, tumor, primary_diagnosis, primary_site,
        primary_diagnosis_simplified, primary_diagnosis_very_simplified)) |>
    dplyr::filter(tumor != 'PCPG' |
                    (tumor == 'PCPG' & primary_site != 'Other/Unknown')) |>
    dplyr::left_join(tissue_subtype_codes) |>
    dplyr::rename(site_diagnosis_code = tissue_code_2) |>
    dplyr::rename(icd10_code = icd_10_code) |>
    dplyr::mutate(tumor_stage = ajcc_pathologic_stage) |>
    dplyr::mutate(tumor_stage_TNM = paste0(
      ajcc_pathologic_t,"__",ajcc_pathologic_n,"__",ajcc_pathologic_m)) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "LUAD",
      "Lung",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "OV",
      "Ovary/Fallopian Tube",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "GBM",
      "CNS/Brain",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "UCEC",
      "Uterus",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "BRCA",
      "Breast",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "CHOL",
      "Biliary Tract",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) & tumor == "TGCT",
      "Testis",
      as.character(primary_site)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      is.na(primary_site) &
        (tumor == "COAD" | tumor == "READ"),
      "Colon/Rectum",
      as.character(primary_site)
    ))

  ## Make slim version of uniform clinical data frame
  tcga_clinical_complete_slim <- dplyr::select(
    tcga_clinical_complete, bcr_patient_barcode, tumor,
    primary_site, site_diagnosis_code, primary_diagnosis,
    primary_diagnosis_simplified, primary_diagnosis_very_simplified,
    icd10_code, tumor_stage, tumor_stage_TNM,
    days_to_death, days_to_last_follow_up, death_event, new_death,
    prior_malignancy, vital_status, tissue_or_organ_of_origin,
    site_of_resection_or_biopsy, year_of_birth, age_at_index, gender, ethnicity, race,
    cigarettes_per_day, alcohol_history, pack_years_smoked,
    MSI_status, Gleason_score, ER_status, PR_status, HER2_status,
    pancan_subtype_mrna, pancan_subtype_cna, pancan_subtype_dnameth,
    pancan_subtype_other, pancan_subtype_protein, pancan_subtype_integrative,
    pancan_subtype_mirna, pancan_subtype_selected)


  tcga_clinical_info <- list()
  tcga_clinical_info[['complete']] <- tcga_clinical_complete
  tcga_clinical_info[['slim']] <- tcga_clinical_complete_slim

  if(!dir.exists(file.path(output_dir, tcga_release))){
    dir.create(file.path(output_dir, tcga_release))
  }

  if(!dir.exists(file.path(output_dir, tcga_release, "clinical"))){
    dir.create(file.path(output_dir, tcga_release, "clinical"))
  }

  saveRDS(
    tcga_clinical_info,
    file = file.path(
      output_dir,
      tcga_release,
      "clinical",
      "tcga_clinical.rds"
    )
  )

  return(tcga_clinical_info)

}


