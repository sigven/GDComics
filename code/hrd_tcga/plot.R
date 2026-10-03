suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
has_pw <- requireNamespace("patchwork", quietly = TRUE)
d <- read.delim("tcga_hrd_scores.tsv", stringsAsFactors = FALSE)
pj <- c("OV", "BRCA", "PRAD", "PAAD")
n_pj <- table(d$project)[pj]
d$project <- factor(d$project, pj, sprintf("%s (n=%d)", pj, n_pj))
d$brca12_group <- factor(d$brca12_group, c("No BRCA1/2 alteration",
  "BRCA1/2 monoallelic (hemdel/cnLOH/LoF)", "BRCA1/2 biallelic (homdel or LoF+LOH)"),
  c("none", "mono-allelic\n(hemdel/cnLOH/LoF)", "bi-allelic\n(homdel/LoF+LOH)"))
cols <- c("none" = "#7f8c8d", "mono-allelic\n(hemdel/cnLOH/LoF)" = "#e69f00", "bi-allelic\n(homdel/LoF+LOH)" = "#c0392b")
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "none")
thr <- geom_hline(yintercept = 42, linetype = "dashed", colour = "grey30")

pA <- ggplot(d, aes(hrd_sum)) + geom_histogram(binwidth = 4, fill = "#4a6fa5", colour = "white") +
  geom_vline(xintercept = 42, linetype = "dashed") + facet_wrap(~project, ncol = 2, scales = "free_y") +
  labs(x = "HRD score (HRD-LOH + LST + TAI)", y = "Samples", title = "A. Score distribution") + th
pB <- ggplot(d, aes(brca12_group, hrd_sum, fill = brca12_group)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .2, size = .5, alpha = .4) +
  thr + scale_fill_manual(values = cols) + facet_wrap(~project, ncol = 2) +
  labs(x = NULL, y = "HRD score", title = "B. By BRCA1/2 status") + th + theme(axis.text.x = element_text(size = 7))
g <- d |> select(sample, project, hrd_sum, BRCA1 = BRCA1_cn_status, BRCA2 = BRCA2_cn_status) |>
  pivot_longer(c(BRCA1, BRCA2), names_to = "gene", values_to = "cn") |>
  mutate(cn = factor(cn, c("intact", "cnLOH", "HEMDEL", "HOMDEL")))
pC <- ggplot(g, aes(cn, hrd_sum, fill = cn)) + geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .2, size = .5, alpha = .4) + thr +
  scale_fill_manual(values = c(intact = "#7f8c8d", cnLOH = "#9ab973", HEMDEL = "#e69f00", HOMDEL = "#c0392b")) +
  facet_grid(project ~ gene) + labs(x = "Copy number at gene (allele-specific)", y = "HRD score",
  title = "C. By BRCA1/BRCA2 locus CN state") + th
b <- d |> filter(grepl("^BRCA", project), !is.na(subtype_selected), subtype_selected != "")
pr <- d |> filter(grepl("^PRAD", project), !is.na(gleason_score)) |> mutate(gleason_score = factor(gleason_score))
pE <- ggplot(pr, aes(gleason_score, hrd_sum, fill = gleason_score)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .2, size = .5, alpha = .4) + thr +
  labs(x = "Gleason score (PRAD)", y = "HRD score", title = "E. TCGA-PRAD by Gleason score") + th
pD <- ggplot(b, aes(reorder(subtype_selected, -hrd_sum, median), hrd_sum, fill = subtype_selected)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .2, size = .5, alpha = .4) + thr +
  labs(x = "PAM50 subtype (BRCA)", y = "HRD score", title = "D. TCGA-BRCA by subtype") + th
if (has_pw) {
  library(patchwork); p <- (pA | pB) / (pC | (pD / pE)) + plot_layout(heights = c(1, 1.3))
} else { p <- gridExtra::grid.arrange(pA, pB, pC, pD, ncol = 2) }
ggsave("hrd_tcga_cohorts.png", p, width = 14, height = 17, dpi = 150)
cat("saved; patchwork:", has_pw, "\n")
