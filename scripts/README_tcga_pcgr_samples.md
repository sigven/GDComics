# TCGA test samples for PCGR

Scripts to select biomarker-rich TCGA samples and turn them into "PCGR-ready"
input files, based on a GDComics release (`output/<release>/`).

Requires that `scripts/gdc_preprocess.R` has been run for the release (CNA
incl. ASCAT3 segments/purity-ploidy, RNA fusions, mutations, MSI, TMB).

## 1. Select sample sets — `select_tcga_samples.R`

```bash
Rscript scripts/select_tcga_samples.R
```

Picks 10 samples per tumor type (STAD, LUAD, PRAD, BRCA, COAD, PAAD, SKCM,
BLCA, UCEC, HNSC) using `select_tcga_samples()` in `code/sample_selection.R`:

- **Biomarkers** are defined in the curated table
  `data-raw/sample_lists/tcga_biomarkers_curated.tsv` (mutations by protein
  change pattern, loss-of-function, AMPL/HOMDEL, fusions, MSI-H, TMB-high,
  HER2+/triple-negative). Edit this file to change what counts as interesting.
- **Candidates**: primary (01) or metastatic (06) samples, ASCAT QC `Pass`,
  purity ≥ 0.4, with SNV/InDel calls.
- **Selection**: greedy coverage of the biomarker list (rare biomarkers
  first, up to two samples per biomarker), max 2 hypermutated samples
  (TMB ≥ 20), ties broken by purity closest to 0.7, plus ≥ 1
  biomarker-negative control.

Output: `data-raw/sample_lists/<release>/tcga_<tumor>_biomarker_set.tsv`
(one per tumor type, with `selection_role`, `biomarkers`, `evidence`,
purity, ploidy, TMB, MSI) and `tcga_biomarker_sets_all.tsv`.

## 2. Prepare PCGR input — `prepare_tcga_pcgr_batch.py` / `prepare_tcga_pcgr.py`

```bash
# a set of samples (plain list of barcodes, or a TSV with a sample ID column)
scripts/prepare_tcga_pcgr_batch.py data-raw/sample_lists/release46_20260810/tcga_luad_biomarker_set.tsv 46 \
  --output_dir output/release46_20260810/pcgr_input/luad --grch37

# a single sample (patient, sample or sample+vial barcode)
scripts/prepare_tcga_pcgr.py TCGA-29-2429-01 46 --grch37
```

Per sample:

| File | Content |
|---|---|
| `<sample>.grch38.vcf.gz(.tbi)` | Somatic SNVs/InDels (tumor-normal). INFO: `TDP`, `TVAF`, `CDP`, `TAL` |
| `<sample>.grch38.cna.tsv` | ASCAT3 allele-specific segments (`Chromosome Start End nMajor nMinor`) |
| `<sample>.grch38.fusions.tsv` | RNA fusions (`FusionGene LeftBreakpoint RightBreakpoint SplitReads` + `SpanningReads FrameEffect FusionLabel`) |
| `<sample>.grch37.*` | Same, on GRCh37 (with `--grch37`) |
| `<sample>.gene_expression.tsv` | RNA-seq gene expression (`TargetID` = Ensembl gene ID, `TPM`), assembly-independent |
| `<sample>.pcgr_sample_info.tsv` | Project, site, counts, purity/ploidy (ASCAT), QC |

A batch run also writes `pcgr_sample_info.tsv` for all samples.

PCGR tags: `--tumor_dp_tag TDP --tumor_af_tag TVAF --control_dp_tag CDP`;
purity/ploidy for `--tumor_purity`/`--tumor_ploidy` are in the sample info file;
expression via `--input_rna_expression <sample>.gene_expression.tsv`.

Options: `--purity_ploidy` (ASCAT summary file, default
`data-raw/ascat3/summary.ascatv3TCGA.penalty70.hg38.tsv`),
`--mnv_max_vaf_diff` (default 0.05), `--no_mnv_merge`, `--no_expression`.

## Data sources

- **SNVs/InDels**: `output/<release>/snv_indel/tcga_mutation_grch38.tsv.gz` (GDC MAFs)
- **CNA segments**: GDC ASCAT3 allele-specific segments (`data-raw/ascat3/202502/sample_data/`)
- **Purity/ploidy**: [ASCAT TCGA release](https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData/TCGA_SNP6_hg38) (one representative sample per case); ploidy is also computed from the segments
- **Expression**: GDC STAR - Counts TPM, `output/<release>/rnaseq/tcga_rnaseq_TPM_<tumor>.rds`
  (RNA-seq aliquot matched on sample barcode; recorded as `rna_sample_barcode`)
- **Fusions**: TCGA PanCancer Atlas (Gao et al., 2018) via cBioPortal (`output/<release>/fusion/`), original coordinates are GRCh38, lifted to GRCh37 where needed

## Caveats

- **BRAF numbering**: GDC annotates BRAF on ENST00000288602, where V600E is `p.V640E`.
- **Split MNVs**: GDC MAFs report multi-base substitutions as separate SNVs
  (e.g. BRAF V600K = `p.V640E` + `p.V640M`). The PCGR VCFs merge adjacent
  SNVs with similar VAF into MNVs (INFO `MNV_SNVS`); sample selection matches
  such codons as combined changes (`p.V640E&p.V640M`).
- **Somatic only**: no germline variants (e.g. BRCA1/2 carriers are not captured).
- **Fusions**: no fusion file for a sample means either no fusion passed
  the PanCancer Atlas filters or no RNA-seq; the data does not distinguish these.
