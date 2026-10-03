## MAF variant classes used for biomarker matching
truncating_classes <- c(
  "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "Splice_Site")
nonsilent_classes <- c(
  truncating_classes, "Missense_Mutation", "In_Frame_Del",
  "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site")

#' Combine missense calls hitting the same codon in a sample
#'
#' GDC MAFs report multi-base substitutions within a codon as separate SNVs
#' (e.g. BRAF V600K = p.V640E + p.V640M on ENST00000288602). Missense calls
#' at the same protein position, at most 2 bp apart and with tumor VAFs
#' differing by at most max_vaf_diff, get a combined protein change label
#' (sorted, joined by '&', e.g. 'p.V640E&p.V640M') in HGVSp_match. Other
#' calls keep HGVSp_Short.
#'
#' @param muts Data frame of mutation calls (GDComics snv_indel)
#' @param max_vaf_diff Maximum difference in tumor VAF
#'
combine_split_codons <- function(muts = NULL, max_vaf_diff = 0.05){

  muts <- muts |>
    dplyr::mutate(
      vaf = t_alt_count / t_depth,
      aa_pos = as.integer(stringr::str_match(
        HGVSp_Short, "^p\\.[A-Z](\\d+)[A-Z]$")[, 2]))

  split_codons <- muts |>
    dplyr::filter(!is.na(aa_pos)) |>
    dplyr::group_by(tumor_sample_barcode, Hugo_Symbol, aa_pos) |>
    dplyr::summarise(
      n = dplyr::n(),
      span = max(POS) - min(POS),
      vaf_diff = max(vaf) - min(vaf),
      combined = paste(sort(unique(HGVSp_Short)), collapse = "&"),
      .groups = "drop") |>
    dplyr::filter(n > 1, span <= 2, vaf_diff <= max_vaf_diff)

  muts |>
    dplyr::left_join(
      dplyr::select(split_codons, tumor_sample_barcode, Hugo_Symbol,
                    aa_pos, combined),
      by = c("tumor_sample_barcode", "Hugo_Symbol", "aa_pos")) |>
    dplyr::mutate(HGVSp_match = dplyr::coalesce(combined, HGVSp_Short)) |>
    dplyr::distinct(tumor_sample_barcode, Hugo_Symbol,
                    Variant_Classification, HGVSp_match)

}

#' Read curated biomarker definitions for sample set selection
#'
#' @param biomarkers_tsv Path to curated biomarker TSV
#'
#' @export
#'
read_curated_biomarkers <- function(biomarkers_tsv = NULL){

  readr::read_tsv(
    biomarkers_tsv,
    comment = "#",
    col_types = readr::cols(.default = "c"),
    na = c("", "NA")) |>
    dplyr::mutate(
      threshold = as.numeric(threshold))

}

#' Find samples carrying each curated biomarker for a TCGA tumor type
#'
#' @param tumor TCGA tumor code (e.g. "LUAD")
#' @param biomarkers Data frame of curated biomarkers for the tumor type
#' @param release_dir GDComics release output directory
#' @param msi MSI calls (msi/tcga_msi.rds)
#' @param tmb TMB estimates (tmb/tcga_tmb.rds)
#' @param receptor ER/PR/HER2 status (clinical/tcga_er_pr_her2.rds)
#' @param fusions RNA fusions (fusion/tcga_fusions.rds)
#'
#' @return data frame with tumor_sample_barcode, biomarker, evidence
#'
get_biomarker_hits <- function(
    tumor = NULL,
    biomarkers = NULL,
    release_dir = NULL,
    msi = NULL,
    tmb = NULL,
    receptor = NULL,
    fusions = NULL){

  project <- paste0("TCGA-", tumor)
  muts <- readRDS(file.path(
    release_dir, "snv_indel",
    glue::glue("tcga_mutation_{tumor}_grch38.rds"))) |>
    dplyr::filter(Variant_Classification %in% nonsilent_classes) |>
    dplyr::select(tumor_sample_barcode, Hugo_Symbol, Variant_Classification,
                  HGVSp_Short, POS, t_alt_count, t_depth) |>
    combine_split_codons()
  cna <- readRDS(file.path(
    release_dir, "cna",
    glue::glue("tcga_cna_ASCAT3_{tumor}.rds"))) |>
    dplyr::filter(mut_status %in% c("AMPL", "HOMDEL")) |>
    dplyr::select(tumor_sample_barcode, symbol, mut_status, copy_number)
  project_fusions <- fusions |>
    dplyr::filter(project_id == project)

  hits <- purrr::pmap_dfr(biomarkers, function(biomarker, alteration_type, gene,
                                               variant_class, protein_regex,
                                               threshold, ...){
    res <- switch(
      alteration_type,
      "MUT" = {
        classes <- if(is.na(variant_class) || variant_class == "NONSILENT"){
          nonsilent_classes
        }else if(variant_class == "TRUNCATING"){
          truncating_classes
        }else{
          stringr::str_split_1(variant_class, ",")
        }
        m <- muts |>
          dplyr::filter(Hugo_Symbol == gene,
                        Variant_Classification %in% classes)
        if(!is.na(protein_regex)){
          m <- m |>
            dplyr::filter(stringr::str_detect(
              dplyr::coalesce(HGVSp_match, ""), protein_regex))
        }
        m |> dplyr::transmute(
          tumor_sample_barcode,
          evidence = paste(gene, HGVSp_match))
      },
      "AMPL" = , "HOMDEL" = {
        cna |>
          dplyr::filter(symbol == gene, mut_status == alteration_type) |>
          dplyr::transmute(
            tumor_sample_barcode,
            evidence = paste0(gene, " ", alteration_type, " (CN=", copy_number, ")"))
      },
      "FUSION" = {
        partners <- stringr::str_split_1(gene, "--")
        project_fusions |>
          dplyr::filter(
            (partners[1] == "*" | gene5 %in% partners[1]),
            (partners[2] == "*" | gene3 %in% partners[2])) |>
          dplyr::transmute(
            tumor_sample_barcode,
            evidence = paste0(gene5, "--", gene3, " (", split_reads, " split reads)"))
      },
      "MSI_H" = {
        msi |>
          dplyr::filter(msi_status == "MSI-H") |>
          dplyr::transmute(tumor_sample_barcode, evidence = "MSI-H")
      },
      "TMB_HIGH" = {
        tmb |>
          dplyr::filter(TMB_coding_non_silent >= threshold) |>
          dplyr::transmute(
            tumor_sample_barcode,
            evidence = paste0("TMB=", round(TMB_coding_non_silent, 1)))
      },
      "HER2_POSITIVE" = , "TRIPLE_NEGATIVE" = {
        r <- if(alteration_type == "HER2_POSITIVE"){
          dplyr::filter(receptor, her2_status == "Positive")
        }else{
          dplyr::filter(receptor, er_status == "Negative",
                        pr_status == "Negative", her2_status == "Negative")
        }
        r |> dplyr::transmute(
          bcr_patient_barcode, evidence = alteration_type)
      },
      stop("Unknown alteration_type: ", alteration_type)
    )
    dplyr::mutate(res, biomarker = biomarker)
  })

  return(hits)

}

#' Select a set of biomarker-rich TCGA samples per tumor type
#'
#' Candidates are tumors with given sample type codes (default primary (01)
#' and metastatic (06)) with ASCAT QC 'Pass', purity >= min_purity and
#' somatic SNV/InDel calls. Samples are picked greedily to cover as many
#' curated biomarkers as possible (first one sample per biomarker, then a
#' second), rare biomarkers weighted higher, ties broken by purity closest
#' to target_purity (avoiding purity = 1 fits, which may reflect profiles fitted
#' on normal cells), with a cap on hypermutated samples
#' (which otherwise dominate through incidental hits). Remaining slots
#' are given to biomarker-negative control samples.
#'
#' @param tumors Character vector of TCGA tumor codes
#' @param biomarkers_tsv Path to curated biomarker TSV
#' @param output_dir GDComics output directory
#' @param gdc_release TCGA data release version
#' @param sample_list_dir Directory for sample set TSVs (one per tumor type)
#' @param sample_type_codes TCGA sample type codes allowed (01 = primary,
#' 06 = metastatic)
#' @param n_samples Number of samples per tumor type
#' @param n_controls Number of biomarker-negative control samples per set
#' @param min_purity Minimum ASCAT tumor purity
#' @param target_purity Preferred purity for tie-breaks and controls
#' @param tmb_hypermutated TMB (mut/Mb) defining hypermutated samples
#' @param max_hypermutated Maximum number of hypermutated samples per set
#'
#' @export
#'
select_tcga_samples <- function(
    tumors = NULL,
    biomarkers_tsv = NULL,
    output_dir = NULL,
    gdc_release = NULL,
    sample_list_dir = NULL,
    sample_type_codes = c("01", "06"),
    n_samples = 10,
    n_controls = 1,
    min_purity = 0.4,
    target_purity = 0.7,
    tmb_hypermutated = 20,
    max_hypermutated = 2){

  release_dir <- file.path(output_dir, gdc_release)
  if(!dir.exists(sample_list_dir)){
    dir.create(sample_list_dir, recursive = T)
  }

  biomarkers <- read_curated_biomarkers(biomarkers_tsv)
  purity_ploidy <- readRDS(file.path(
    release_dir, "cna", "tcga_cna_purity_ploidy_ASCAT3.rds"))
  tmb <- readRDS(file.path(release_dir, "tmb", "tcga_tmb.rds"))
  msi <- readRDS(file.path(release_dir, "msi", "tcga_msi.rds"))
  receptor <- readRDS(file.path(release_dir, "clinical", "tcga_er_pr_her2.rds"))
  fusions <- readRDS(file.path(release_dir, "fusion", "tcga_fusions.rds"))
  patient_sex <- readRDS(file.path(release_dir, "clinical", "tcga_clinical.rds")) |>
    dplyr::distinct(bcr_patient_barcode, sex = toupper(sex_at_birth)) |>
    dplyr::distinct(bcr_patient_barcode, .keep_all = TRUE)

  sample_sets <- list()
  for(tumor in tumors){

    candidates <- purity_ploidy |>
      dplyr::filter(
        project_id == paste0("TCGA-", tumor),
        substr(tumor_sample_barcode, 14, 15) %in% sample_type_codes,
        ascat_qc == "Pass",
        ascat_purity >= min_purity) |>
      dplyr::inner_join(
        dplyr::select(tmb, tumor_sample_barcode, tmb = TMB_coding_non_silent),
        by = "tumor_sample_barcode") |>
      dplyr::left_join(
        dplyr::select(msi, tumor_sample_barcode, msi_status),
        by = "tumor_sample_barcode") |>
      dplyr::transmute(
        tumor_sample_barcode,
        bcr_patient_barcode = substr(tumor_sample_barcode, 1, 12),
        project_id, purity = ascat_purity, ploidy, tmb = round(tmb, 2),
        msi_status, hypermutated = tmb >= tmb_hypermutated) |>
      dplyr::distinct(tumor_sample_barcode, .keep_all = TRUE)

    tumor_biomarkers <- dplyr::filter(biomarkers, tumor == !!tumor)
    hits_raw <- get_biomarker_hits(
      tumor = tumor, biomarkers = tumor_biomarkers,
      release_dir = release_dir, msi = msi, tmb = tmb,
      receptor = receptor, fusions = fusions)

    ## Harmonize barcodes: receptor status is per patient, fusions per sample
    ## (15-char), other data per sample + vial (16-char)
    for(col in c("bcr_patient_barcode", "tumor_sample_barcode")){
      if(!(col %in% names(hits_raw))) hits_raw[[col]] <- NA_character_
    }
    candidate_ids <- candidates |>
      dplyr::select(bcr_patient_barcode, tumor_sample_barcode) |>
      dplyr::mutate(sample15 = substr(tumor_sample_barcode, 1, 15))
    hits <- dplyr::bind_rows(
      hits_raw |>
        dplyr::filter(!is.na(bcr_patient_barcode)) |>
        dplyr::select(-tumor_sample_barcode) |>
        dplyr::inner_join(
          candidate_ids, by = "bcr_patient_barcode"),
      hits_raw |>
        dplyr::filter(is.na(bcr_patient_barcode)) |>
        dplyr::mutate(sample15 = substr(tumor_sample_barcode, 1, 15)) |>
        dplyr::select(-c("tumor_sample_barcode", "bcr_patient_barcode")) |>
        dplyr::inner_join(
          candidate_ids, by = "sample15")) |>
      dplyr::select(tumor_sample_barcode, biomarker, evidence) |>
      dplyr::distinct()

    ## Greedy selection: maximize newly covered biomarkers (coverage round 1,
    ## then 2), each biomarker weighted by 1/number of carriers (rare first)
    sample_biomarkers <- split(hits$biomarker, hits$tumor_sample_barcode)
    n_carriers <- table(dplyr::distinct(hits, tumor_sample_barcode, biomarker)$biomarker)
    weight <- setNames(as.numeric(1 / n_carriers), names(n_carriers))
    coverage <- setNames(rep(0, length(unique(tumor_biomarkers$biomarker))),
                         unique(tumor_biomarkers$biomarker))
    selected <- character()
    n_biomarker_slots <- n_samples - n_controls
    for(coverage_round in 1:2){
      while(length(selected) < n_biomarker_slots){
        pool <- candidates |>
          dplyr::filter(
            tumor_sample_barcode %in% names(sample_biomarkers),
            !(tumor_sample_barcode %in% selected))
        if(sum(candidates$hypermutated[
          candidates$tumor_sample_barcode %in% selected]) >= max_hypermutated){
          pool <- dplyr::filter(pool, !hypermutated)
        }
        if(nrow(pool) == 0) break
        pool$score <- purrr::map_dbl(
          pool$tumor_sample_barcode,
          function(s){
            b <- unique(sample_biomarkers[[s]])
            sum(weight[b[coverage[b] < coverage_round]])
          })
        best <- pool |>
          dplyr::arrange(dplyr::desc(score), abs(purity - target_purity)) |>
          dplyr::slice(1)
        if(best$score == 0) break
        selected <- c(selected, best$tumor_sample_barcode)
        covered <- unique(sample_biomarkers[[best$tumor_sample_barcode]])
        coverage[covered] <- coverage[covered] + 1
      }
    }

    ## Biomarker-negative controls (not hypermutated, purity closest to target)
    controls <- candidates |>
      dplyr::filter(
        !(tumor_sample_barcode %in% names(sample_biomarkers)),
        !hypermutated) |>
      dplyr::arrange(abs(purity - target_purity)) |>
      dplyr::slice_head(n = n_samples - length(selected)) |>
      dplyr::pull(tumor_sample_barcode)

    sample_set <- candidates |>
      dplyr::filter(tumor_sample_barcode %in% c(selected, controls)) |>
      dplyr::mutate(
        selection_role = dplyr::if_else(
          tumor_sample_barcode %in% selected, "biomarker", "control")) |>
      dplyr::left_join(
        hits |>
          dplyr::group_by(tumor_sample_barcode) |>
          dplyr::summarise(
            biomarkers = paste(sort(unique(biomarker)), collapse = ";"),
            evidence = paste(unique(evidence), collapse = "; "),
            .groups = "drop"),
        by = "tumor_sample_barcode") |>
      dplyr::mutate(
        tumor = tumor,
        bcr_patient_barcode = substr(tumor_sample_barcode, 1, 12)) |>
      dplyr::left_join(patient_sex, by = "bcr_patient_barcode") |>
      dplyr::mutate(sex = dplyr::coalesce(sex, "UNKNOWN")) |>
      dplyr::arrange(
        factor(selection_role, levels = c("biomarker", "control")),
        match(tumor_sample_barcode, selected)) |>
      dplyr::select(
        tumor_sample_barcode, tumor, sex, selection_role, biomarkers,
        evidence, purity, ploidy, tmb, msi_status, hypermutated)

    uncovered <- names(coverage)[coverage == 0]
    available <- unique(hits$biomarker)
    cat(glue::glue(
      "{tumor}: {length(candidates$tumor_sample_barcode)} candidates, ",
      "{nrow(sample_set)} selected, ",
      "{sum(coverage > 0)}/{length(coverage)} biomarkers covered"), sep = "\n")
    if(length(uncovered) > 0){
      cat(paste0(
        "  not covered: ",
        paste(uncovered, ifelse(uncovered %in% available,
                                "(present, crowded out)", "(no candidate)"),
              collapse = ", ")), sep = "\n")
    }

    readr::write_tsv(
      sample_set,
      file = file.path(
        sample_list_dir,
        glue::glue("tcga_{tolower(tumor)}_biomarker_set.tsv")))
    sample_sets[[tumor]] <- sample_set
  }

  all_sets <- dplyr::bind_rows(sample_sets)
  readr::write_tsv(
    all_sets,
    file = file.path(sample_list_dir, "tcga_biomarker_sets_all.tsv"))
  return(all_sets)

}
