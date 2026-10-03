R <- "/Users/sigven/project_data/packages/package__GDComics/GDComics/output/release46_20260810"
suppressMessages(library(dplyr))
cl <- readRDS(file.path(R, "clinical/tcga_clinical.rds"))
meta <- list()
projects <- c("BRCA", "OV", "PRAD", "PAAD")
for (p in projects) {
  seg <- readRDS(file.path(R, "cna", sprintf("tcga_cna_segments_ASCAT3_%s.rds", p)))
  pp  <- readRDS(file.path(R, "cna", sprintf("tcga_cna_purity_ploidy_ASCAT3_%s.rds", p)))
  seg <- seg |> transmute(
    sample = tumor_sample_barcode, Chromosome = sub("^chr", "", chromosome),
    Start = start, End = end, nMajor, nMinor)
  saveRDS(seg, sprintf("segments_%s.rds", p))
  write.table(seg, sprintf("segments_%s.tsv", p), sep = "\t", quote = FALSE, row.names = FALSE)
  mut <- readRDS(file.path(R, "snv_indel", sprintf("tcga_mutation_%s_grch38.rds", p))) |>
    filter(Hugo_Symbol %in% c("BRCA1", "BRCA2"),
           Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del",
                                         "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site")) |>
    group_by(sample = tumor_sample_barcode) |>
    summarise(somatic_lof = paste(unique(paste0(Hugo_Symbol, ":", HGVSp_Short)), collapse = ";"),
              somatic_lof_genes = paste(sort(unique(Hugo_Symbol)), collapse = ","), .groups = "drop")
  m <- pp |> transmute(sample = tumor_sample_barcode, project = p, ploidy_segments = ploidy,
                       ascat_purity, ascat_ploidy, ascat_wgd, ascat_qc) |>
    left_join(cl |> filter(project_id == paste0("TCGA-", p)) |>
                distinct(tumor_sample_barcode, .keep_all = TRUE) |>
                transmute(sample = tumor_sample_barcode, er_status, her2_status, pr_status,
                          subtype_selected, gleason_score, primary_diagnosis_very_simplified),
              by = "sample") |>
    left_join(mut, by = "sample")
  meta[[p]] <- m
}
write.table(bind_rows(meta), "sample_meta.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("done\n")
