#!/usr/bin/env python

"""
Prepare "PCGR-ready" input files for TCGA tumor samples, based on processed
GDComics output (somatic SNVs/InDels) and GDC ASCAT3 allele-specific copy
number segments.

Per sample, the following files are written to the output directory:

  <sample>.grch38.vcf.gz(.tbi)  - somatic SNVs/InDels (tumor-normal), GRCh38
  <sample>.grch38.cna.tsv       - allele-specific CNA segments (PCGR format), GRCh38
  <sample>.grch38.fusions.tsv   - RNA fusions (TCGA PanCancer Atlas, via cBioPortal), GRCh38 (original coordinates)
  <sample>.grch37.fusions.tsv   - RNA fusions, GRCh37 (lifted over; only with --grch37)
  <sample>.grch37.*             - lifted-over versions (only with --grch37)
  <sample>.gene_expression.tsv  - RNA-seq gene expression (TargetID = Ensembl gene ID, TPM),
                                  assembly-independent (GDC STAR - Counts, GDComics rnaseq output)
  <sample>.pcgr_sample_info.tsv - sample metadata (project, site, purity, ploidy, ..)

VCF INFO tags useful for PCGR:
  --tumor_dp_tag TDP --tumor_af_tag TVAF --control_dp_tag CDP
Expression: --input_rna_expression <sample>.gene_expression.tsv

Usage:
  prepare_tcga_pcgr.py TCGA-29-2429-01 release46_20260810
  prepare_tcga_pcgr.py TCGA-29-2429 46 --grch37
"""

import argparse
import csv
import glob
import gzip
import os
import re
import subprocess
import sys
import tempfile

gdcomics_dir = '/Users/sigven/project_data/packages/package__GDComics/GDComics'
crossmap_vcf_script = '/Users/sigven/project_data/packages/package__pcgr/test/scripts/crossmap_vcf.R'
cna_lift_script = '/Users/sigven/project_data/data/data__impress/impress/sample_data/lift_cna_segments.py'

## GDC ASCAT3 raw data (gene-level metadata/downloads + allele-specific segments)
cna_metadata_fname = os.path.join(gdcomics_dir, 'data-raw', 'GDCdata', 'cna', 'cna_gdc.metadata.tsv.gz')
cna_gene_level_dir = os.path.join(gdcomics_dir, 'data-raw', 'GDCdata', 'cna')
ascat_segment_dir = os.path.join(gdcomics_dir, 'data-raw', 'ascat3', '202502', 'sample_data')
## ASCAT TCGA SNP6 release (purity/ploidy), see get_ascat_tcga_release() in code/cna.R
ascat_summary_fname = os.path.join(gdcomics_dir, 'data-raw', 'ascat3', 'summary.ascatv3TCGA.penalty70.hg38.tsv')

chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']

algorithm_codes = {'muse': 1, 'mutect2': 2, 'somaticsniper': 3, 'varscan2': 4, 'pindel': 5}


def resolve_release_dir(release: str) -> str:
    """Accept either full release name (e.g. 'release46_20260810') or just the number ('46')."""
    output_dir = os.path.join(gdcomics_dir, 'output')
    if re.fullmatch(r'\d+', release):
        hits = sorted(glob.glob(os.path.join(output_dir, f'release{release}_*')))
        if not hits:
            sys.exit(f'ERROR: no GDComics output found for release {release} in {output_dir}')
        return hits[-1]
    release_dir = os.path.join(output_dir, release)
    if not os.path.isdir(release_dir):
        sys.exit(f'ERROR: release directory {release_dir} does not exist')
    return release_dir


def validate_barcode(barcode: str) -> str:
    """Allow patient (TCGA-XX-XXXX), sample (TCGA-XX-XXXX-01) or sample+vial (TCGA-XX-XXXX-01A) barcodes."""
    barcode = barcode.strip().upper()
    if not re.fullmatch(r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(-[0-9]{2}[A-Z]?)?', barcode):
        sys.exit(f'ERROR: invalid TCGA barcode: {barcode}')
    return barcode


def select_sample(query: str, candidates: list) -> str | None:
    """
    Pick the tumor sample barcode (16-char, e.g. TCGA-29-2429-01A) that matches the
    query barcode. If several match (e.g. primary + recurrence), prefer primary tumor (01).
    """
    hits = sorted(set(c for c in candidates if c.startswith(query)))
    if not hits:
        return None
    primary = [h for h in hits if h[13:15] == '01']
    chosen = primary[0] if primary else hits[0]
    if len(hits) > 1:
        print(f'  NOTE: multiple tumor samples match {query} ({", ".join(hits)}) - using {chosen}')
    return chosen


def read_mutations(mutation_tsv: str, patient_ids: set) -> dict:
    """
    Read somatic calls for a set of patients from the flat GDComics mutation TSV.
    Rows are pre-filtered on the first column (bcr_patient_barcode) before parsing, which is
    much faster than parsing all rows (and than BSD grep with multiple patterns on macOS).
    Returns dict: tumor_sample_barcode -> list of row dicts
    """
    gz = subprocess.Popen(['gzip', '-dc', mutation_tsv], stdout=subprocess.PIPE, text=True)
    header = gz.stdout.readline()
    matching_lines = (line for line in gz.stdout if line[:12] in patient_ids)

    calls = {}
    reader = csv.DictReader(matching_lines, fieldnames=header.rstrip('\n').split('\t'), delimiter='\t')
    for row in reader:
        calls.setdefault(row['tumor_sample_barcode'], []).append(row)
    gz.wait()
    return calls


def read_cna_sample_map() -> dict:
    """
    Map tumor sample barcode -> dict(project, segment file), using the GDC ASCAT3 gene-level
    metadata (file_id -> barcode) and the downloaded gene-level files, whose names carry the
    same GDC aliquot UUID as the allele-specific segment files.
    """
    sample_map = {}
    with gzip.open(cna_metadata_fname, mode='rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene_level = glob.glob(os.path.join(cna_gene_level_dir, row['project_id'], row['file_id'], '*.ascat3.gene_level_copy_number*'))
            if not gene_level:
                continue
            aliquot_uuid = os.path.basename(gene_level[0]).split('.')[1]
            seg_fname = os.path.join(ascat_segment_dir, f"{row['project_id']}.{aliquot_uuid}.ascat3.allelic_specific.seg.txt")
            if os.path.exists(seg_fname):
                sample_map[row['tumor_sample_barcode']] = {'project': row['project_id'], 'segment_fname': seg_fname}
    return sample_map


def read_purity_ploidy(summary_fname: str) -> dict:
    """
    Read purity/ploidy from the ASCAT TCGA SNP6 release
    (https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData/TCGA_SNP6_hg38,
    summary.ascatv3TCGA.penalty70.hg38.tsv). Keyed by 16-char tumor sample barcode.
    """
    purity_ploidy = {}
    if summary_fname is None or not os.path.exists(summary_fname):
        return purity_ploidy
    with open(summary_fname, mode='rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            purity_ploidy[row['barcodeTumour'][:16]] = {
                'purity': row['purity'], 'ploidy': row['ploidy'], 'ascat_qc': row['QC'],
                'ascat_tumor_aliquot': row['barcodeTumour']}
    return purity_ploidy


def read_fusions(fusion_tsv: str) -> dict:
    """
    Read TCGA PanCancer Atlas RNA fusions (GDComics release, see get_tcga_fusions() in
    code/fusion.R). Keyed by 15-char sample barcode (e.g. TCGA-29-2429-01).
    """
    fusions = {}
    if not os.path.exists(fusion_tsv):
        print(f'WARNING: fusion file {fusion_tsv} not found - no fusion output')
        return fusions
    with gzip.open(fusion_tsv, mode='rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            fusions.setdefault(row['tumor_sample_barcode'], []).append(row)
    return fusions


def write_fusions(fusions: list, fusion_fname: str, assembly: str) -> int:
    """
    Write fusions in PCGR format (FusionGene, LeftBreakpoint, RightBreakpoint, SplitReads),
    with additional annotations (SpanningReads, FrameEffect, FusionLabel). Fusions without
    gene symbols for both partners, or without breakpoints in the given assembly, are skipped.
    Only the entry with the most split reads is kept per gene pair.
    """
    best = {}
    for f in fusions:
        pos5, pos3 = f[f'pos5_{assembly}'], f[f'pos3_{assembly}']
        if 'NA' in (f['gene5'], f['gene3'], pos5, pos3):
            continue
        split_reads = int(f['split_reads']) if f['split_reads'] != 'NA' else 0
        key = (f['gene5'], f['gene3'])
        if key in best and split_reads <= best[key]['SplitReads']:
            continue
        best[key] = {
            'FusionGene': f"{f['gene5']}--{f['gene3']}",
            'LeftBreakpoint': f"{f['chrom5']}:{int(float(pos5))}",
            'RightBreakpoint': f"{f['chrom3']}:{int(float(pos3))}",
            'SplitReads': split_reads,
            'SpanningReads': f['spanning_reads'],
            'FrameEffect': f['frame_effect'],
            'FusionLabel': f['fusion_label']}
    if not best:
        return 0
    with open(fusion_fname, mode='wt') as fh:
        writer = csv.DictWriter(fh, fieldnames=next(iter(best.values())).keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(best.values())
    return len(best)


def to_int(value: str):
    return int(float(value)) if value not in ('', 'NA', None) else None


def merge_adjacent_snvs(records: list, max_vaf_diff: float) -> tuple:
    """
    Merge SNVs at consecutive positions (same chromosome) into MNVs when their tumor
    allelic fractions differ by at most max_vaf_diff, i.e. they are likely one event on
    the same haplotype that was split into separate calls (e.g. BRAF V600K/V600E via
    two-base substitutions, UV-induced CC>TT doublets). Records must be sorted by position.
    Positions with more than one SNV call are left unmerged.
    Returns (records, number of MNVs created).
    """
    pos_count = {}
    for r in records:
        pos_count[(r['chrom'], r['pos'])] = pos_count.get((r['chrom'], r['pos']), 0) + 1

    def mergeable(r):
        return (len(r['ref']) == 1 and len(r['alt']) == 1 and r['vaf'] is not None
                and pos_count[(r['chrom'], r['pos'])] == 1)

    merged, n_mnv, i = [], 0, 0
    while i < len(records):
        chain = [records[i]]
        while (mergeable(chain[-1]) and i + len(chain) < len(records)):
            nxt = records[i + len(chain)]
            if not (mergeable(nxt) and nxt['chrom'] == chain[-1]['chrom']
                    and nxt['pos'] == chain[-1]['pos'] + 1
                    and abs(nxt['vaf'] - chain[-1]['vaf']) <= max_vaf_diff):
                break
            chain.append(nxt)
        if len(chain) == 1:
            merged.append(chain[0])
        else:
            def mean(field):
                values = [r[field] for r in chain if r[field] is not None]
                return round(sum(values) / len(values)) if values else None
            mnv = {
                'chrom': chain[0]['chrom'], 'pos': chain[0]['pos'],
                'ref': ''.join(r['ref'] for r in chain), 'alt': ''.join(r['alt'] for r in chain),
                't_depth': mean('t_depth'), 't_alt': mean('t_alt'), 't_ref': mean('t_ref'),
                'n_depth': mean('n_depth'),
                'algorithms': sorted(set(a for r in chain for a in r['algorithms'])),
                'mnv_snvs': len(chain)}
            mnv['vaf'] = mnv['t_alt'] / mnv['t_depth'] if mnv['t_depth'] else None
            merged.append(mnv)
            n_mnv += 1
        i += len(chain)
    return merged, n_mnv


def write_expression(samples: list, release_dir: str, output_dir: str) -> dict:
    """
    Write gene expression (TPM) per sample in PCGR format (TargetID = Ensembl gene ID, TPM),
    from the GDComics rnaseq output (tcga_rnaseq_TPM_<tumor>.rds; one R call per project to
    extract the relevant sample columns). The RNA-seq aliquot is matched on the full sample
    barcode (e.g. TCGA-29-2429-01A), or else on the sample without vial (TCGA-29-2429-01).
    samples: list of (tumor_sample_barcode, project) tuples.
    Returns dict: tumor_sample_barcode -> (expression file, RNA-seq sample barcode, number of genes)
    """
    expression = {}
    by_project = {}
    for tumor_sample, project in samples:
        by_project.setdefault(project.replace('TCGA-', ''), []).append(tumor_sample)

    r_extract = (
        'args <- commandArgs(trailingOnly = TRUE); x <- readRDS(args[1]); '
        'keep <- substr(colnames(x), 1, 15) %in% strsplit(args[3], ",")[[1]]; '
        'x <- x[!is.na(x$ENSEMBL_GENE_ID), c("ENSEMBL_GENE_ID", colnames(x)[keep]), drop = FALSE]; '
        'write.table(x, args[2], sep = "\\t", quote = FALSE, row.names = FALSE, na = "")')

    for tumor, tumor_samples in by_project.items():
        rds_fname = os.path.join(release_dir, 'rnaseq', f'tcga_rnaseq_TPM_{tumor}.rds')
        if not os.path.exists(rds_fname):
            print(f'WARNING: no RNA-seq data ({rds_fname}) - no expression output for {tumor}')
            continue
        print(f'Extracting RNA-seq expression (TPM) for {len(tumor_samples)} {tumor} sample(s)')
        with tempfile.NamedTemporaryFile(suffix='.tsv', dir=output_dir) as tmp:
            subprocess.run(['Rscript', '-e', r_extract, rds_fname, tmp.name,
                            ','.join(sorted(set(s[:15] for s in tumor_samples)))], check=True)
            with open(tmp.name) as fh:
                reader = csv.reader(fh, delimiter='\t')
                header = next(reader)
                rows = list(reader)
        for tumor_sample in tumor_samples:
            rna_cols = sorted(c for c in header[1:] if c[:15] == tumor_sample[:15])
            if not rna_cols:
                continue
            rna_sample = tumor_sample if tumor_sample in rna_cols else rna_cols[0]
            idx = header.index(rna_sample)
            expression_fname = os.path.join(output_dir, f'{tumor_sample}.gene_expression.tsv')
            n_genes = 0
            with open(expression_fname, mode='wt') as fh:
                fh.write('TargetID\tTPM\n')
                for row in rows:
                    if row[idx] != '':
                        fh.write(f'{row[0]}\t{round(float(row[idx]), 4)}\n')
                        n_genes += 1
            expression[tumor_sample] = (expression_fname, rna_sample, n_genes)
    return expression


def write_vcf(calls: list, sample_barcode: str, project: str, vcf_fname: str,
              mnv_max_vaf_diff: float = 0.05) -> tuple:
    """
    Write somatic calls (tumor-normal) for a sample as a sorted, bgzipped and indexed VCF.
    Adjacent SNVs with similar allelic support are merged into MNVs (disable with
    mnv_max_vaf_diff = None). Returns (number of variants, number of MNVs created).
    """
    rank = {c: i for i, c in enumerate(chrom_order)}
    records = {}
    for c in calls:
        key = (c['CHROM'], int(float(c['POS'])), c['REF'], c['ALT'])
        if key in records:
            continue
        t_depth, t_alt = to_int(c['t_depth']), to_int(c['t_alt_count'])
        records[key] = {
            'chrom': key[0], 'pos': key[1], 'ref': key[2], 'alt': key[3],
            't_depth': t_depth, 't_alt': t_alt, 't_ref': to_int(c['t_ref_count']),
            'n_depth': to_int(c['n_depth']),
            'vaf': t_alt / t_depth if t_depth and t_alt is not None else None,
            'algorithms': c['Algorithms'].split(';'), 'mnv_snvs': 1}
    records = sorted(records.values(), key=lambda r: (rank.get(r['chrom'], 99), r['pos'], r['ref'], r['alt']))
    n_mnv = 0
    if mnv_max_vaf_diff is not None:
        records, n_mnv = merge_adjacent_snvs(records, mnv_max_vaf_diff)

    variants = []
    for r in records:
        info = [f'TCGA_CODE={project.replace("TCGA-", "")}', f'TAL={",".join(r["algorithms"])}']
        if r['t_depth']:
            info.append(f'TDP={r["t_depth"]}')
            if r['t_alt'] is not None:
                info.append(f'TVAF={round(r["t_alt"] / r["t_depth"], 4)}')
        if r['n_depth'] is not None:
            info.append(f'CDP={r["n_depth"]}')
        if r['mnv_snvs'] > 1:
            info.append(f'MNV_SNVS={r["mnv_snvs"]}')
        fmt = ':'.join([
            '0/1',
            str(r['t_depth']) if r['t_depth'] is not None else '.',
            str(r['n_depth']) if r['n_depth'] is not None else '.',
            f'{r["t_ref"] if r["t_ref"] is not None else "."},{r["t_alt"] if r["t_alt"] is not None else "."}'])
        variants.append(f'{r["chrom"]}\t{r["pos"]}\t.\t{r["ref"]}\t{r["alt"]}\t.\tPASS\t'
                        + ';'.join(info) + '\tGT:DPT:DPC:ADT\t' + fmt)

    header = [
        '##fileformat=VCFv4.2',
        '##reference=GRCh38',
        f'##source=GDComics ({os.path.basename(os.path.dirname(os.path.dirname(vcf_fname)))})',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=TCGA_CODE,Number=1,Type=String,Description="TCGA project abbreviation">',
        '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">',
        '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Read depth across variant site in tumor">',
        '##INFO=<ID=CDP,Number=1,Type=Integer,Description="Read depth across variant site in control">',
        '##INFO=<ID=TAL,Number=.,Type=String,Description="Algorithms that called the somatic mutation">',
        '##INFO=<ID=MNV_SNVS,Number=1,Type=Integer,Description="Number of adjacent SNV calls (similar allelic support) merged into this MNV; depths/counts are averages">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DPT,Number=1,Type=Integer,Description="Sequencing depth at variant position (tumor)">',
        '##FORMAT=<ID=DPC,Number=1,Type=Integer,Description="Sequencing depth at variant position (control)">',
        '##FORMAT=<ID=ADT,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles (tumor)">',
    ]
    header += [f'##contig=<ID={c}>' for c in chrom_order]
    header.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_barcode]))

    vcf_plain = vcf_fname.removesuffix('.gz')
    with open(vcf_plain, mode='wt') as fh:
        fh.write('\n'.join(header) + '\n')
        fh.write(''.join(v + '\n' for v in variants))
    subprocess.run(['bgzip', '-f', vcf_plain], check=True)
    subprocess.run(['tabix', '-f', '-p', 'vcf', vcf_fname], check=True)
    return len(variants), n_mnv


def write_cna(segment_fname: str, cna_fname: str) -> tuple:
    """
    Convert GDC ASCAT3 allele-specific segments to PCGR CNA format
    (Chromosome, Start, End, nMajor, nMinor). Also returns the segment-based
    ploidy estimate (length-weighted mean total copy number, as defined in ASCAT).
    """
    n_segments = 0
    total_len = 0
    weighted_cn = 0
    with open(segment_fname, mode='rt') as fh_in, open(cna_fname, mode='wt') as fh_out:
        fh_out.write('Chromosome\tStart\tEnd\tnMajor\tnMinor\n')
        for row in csv.DictReader(fh_in, delimiter='\t'):
            if row['Major_Copy_Number'] in ('', 'NA') or row['Minor_Copy_Number'] in ('', 'NA'):
                continue
            chrom = row['Chromosome'].removeprefix('chr')
            start, end = int(row['Start']), int(row['End'])
            n_major, n_minor = int(row['Major_Copy_Number']), int(row['Minor_Copy_Number'])
            fh_out.write(f'{chrom}\t{start}\t{end}\t{n_major}\t{n_minor}\n')
            n_segments += 1
            seg_len = end - start + 1
            total_len += seg_len
            weighted_cn += seg_len * (n_major + n_minor)
    ploidy = round(weighted_cn / total_len, 4) if total_len > 0 else 'NA'
    return n_segments, ploidy


def liftover_grch37(sample_prefix: str) -> None:
    vcf38 = f'{sample_prefix}.grch38.vcf.gz'
    cna38 = f'{sample_prefix}.grch38.cna.tsv'
    if os.path.exists(vcf38):
        subprocess.run([crossmap_vcf_script, vcf38, f'{sample_prefix}.grch37.vcf', 'hg38Tohg19'], check=False)
    if os.path.exists(cna38):
        subprocess.run([cna_lift_script, cna38, f'{sample_prefix}.grch37.cna.tsv', '--direction', 'hg38Tohg19'], check=False)


def prepare_samples(barcodes: list, release: str, output_dir: str,
                    purity_ploidy_fname: str = None, grch37: bool = False,
                    mnv_max_vaf_diff: float = 0.05, expression: bool = True) -> list:
    """
    Prepare PCGR input files for a set of TCGA barcodes. Mutation data is read in a
    single pass over the release-level mutation file. Returns a list of per-sample info dicts.
    """
    release_dir = resolve_release_dir(release)
    mutation_tsv = os.path.join(release_dir, 'snv_indel', 'tcga_mutation_grch38.tsv.gz')
    fusion_tsv = os.path.join(release_dir, 'fusion', 'tcga_fusions.tsv.gz')
    if not os.path.exists(mutation_tsv):
        sys.exit(f'ERROR: mutation file {mutation_tsv} not found')
    os.makedirs(output_dir, exist_ok=True)

    queries = [validate_barcode(b) for b in barcodes]
    patient_ids = set(q[:12] for q in queries)

    print(f'Reading somatic calls for {len(patient_ids)} patient(s) from {mutation_tsv}')
    calls = read_mutations(mutation_tsv, patient_ids)
    print('Mapping ASCAT3 allele-specific segments to tumor samples')
    cna_map = read_cna_sample_map()
    purity_ploidy = read_purity_ploidy(purity_ploidy_fname)
    fusions = read_fusions(fusion_tsv)

    sample_info = []
    for query in queries:
        tumor_sample = select_sample(query, list(calls.keys()) + list(cna_map.keys()))
        if tumor_sample is None:
            print(f'WARNING: no SNV/InDel or CNA data found for {query} - skipping')
            continue
        print(f'Processing {tumor_sample}')
        prefix = os.path.join(output_dir, tumor_sample)
        sample_calls = calls.get(tumor_sample, [])
        project = cna_map[tumor_sample]['project'] if tumor_sample in cna_map else 'TCGA-' + sample_calls[0]['tumor']

        info = {
            'query_barcode': query,
            'tumor_sample_barcode': tumor_sample,
            'project': project,
            'primary_site': sample_calls[0]['primary_site'] if sample_calls else 'NA',
            'primary_diagnosis': sample_calls[0]['primary_diagnosis'] if sample_calls else 'NA',
            'sample_type': sample_calls[0]['sample_type'] if sample_calls else 'NA',
            'n_snv_indel': 0, 'n_mnv_merged': 0, 'n_cna_segments': 0, 'n_fusions': 0,
            'purity': 'NA', 'ploidy': 'NA', 'ploidy_source': 'NA', 'ploidy_segments': 'NA', 'ascat_qc': 'NA',
            'rna_sample_barcode': 'NA', 'n_expressed_genes': 0,
            'vcf_grch38': 'NA', 'cna_grch38': 'NA', 'fusions_grch38': 'NA', 'gene_expression': 'NA',
        }

        if sample_calls:
            vcf_fname = f'{prefix}.grch38.vcf.gz'
            info['n_snv_indel'], info['n_mnv_merged'] = write_vcf(
                sample_calls, tumor_sample, project, vcf_fname, mnv_max_vaf_diff)
            info['vcf_grch38'] = vcf_fname
        else:
            print(f'  WARNING: no SNV/InDel calls found for {tumor_sample}')

        if tumor_sample in cna_map:
            cna_fname = f'{prefix}.grch38.cna.tsv'
            n_segments, seg_ploidy = write_cna(cna_map[tumor_sample]['segment_fname'], cna_fname)
            info.update({'n_cna_segments': n_segments, 'cna_grch38': cna_fname,
                         'ploidy': seg_ploidy, 'ploidy_segments': seg_ploidy, 'ploidy_source': 'ASCAT3_segments'})
        else:
            print(f'  WARNING: no ASCAT3 segments found for {tumor_sample}')

        ## PanCancer Atlas fusions are reported per sample (no vial letter)
        sample_fusions = fusions.get(tumor_sample[:15], [])
        if sample_fusions:
            fusion_fname = f'{prefix}.grch38.fusions.tsv'
            info['n_fusions'] = write_fusions(sample_fusions, fusion_fname, 'grch38')
            if info['n_fusions'] > 0:
                info['fusions_grch38'] = fusion_fname
            if grch37:
                write_fusions(sample_fusions, f'{prefix}.grch37.fusions.tsv', 'grch37')

        if tumor_sample in purity_ploidy:
            pp = purity_ploidy[tumor_sample]
            info.update({'purity': pp['purity'], 'ploidy': pp['ploidy'],
                         'ploidy_source': 'ASCAT_TCGA_SNP6_release', 'ascat_qc': pp['ascat_qc']})

        if grch37:
            liftover_grch37(prefix)

        sample_info.append(info)

    if expression and sample_info:
        sample_expression = write_expression(
            [(i['tumor_sample_barcode'], i['project']) for i in sample_info], release_dir, output_dir)
        for info in sample_info:
            if info['tumor_sample_barcode'] in sample_expression:
                fname, rna_sample, n_genes = sample_expression[info['tumor_sample_barcode']]
                info.update({'gene_expression': fname, 'rna_sample_barcode': rna_sample,
                             'n_expressed_genes': n_genes})
            else:
                print(f"  WARNING: no RNA-seq expression found for {info['tumor_sample_barcode']}")

    for info in sample_info:
        with open(os.path.join(output_dir, f"{info['tumor_sample_barcode']}.pcgr_sample_info.tsv"), mode='wt') as fh:
            writer = csv.DictWriter(fh, fieldnames=info.keys(), delimiter='\t')
            writer.writeheader()
            writer.writerow(info)

    return sample_info


def main():
    parser = argparse.ArgumentParser(
        description='Prepare PCGR-ready input files (SNV/InDel VCF + CNA segments) for a TCGA sample',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('barcode', help='TCGA barcode (patient, e.g. TCGA-29-2429, or sample, e.g. TCGA-29-2429-01A)')
    parser.add_argument('release', help="GDComics release, e.g. 'release46_20260810' or '46'")
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: <GDComics>/output/<release>/pcgr_input)')
    parser.add_argument('--purity_ploidy', default=ascat_summary_fname,
                        help='ASCAT TCGA summary file (summary.ascatv3TCGA.penalty70.hg38.tsv) with purity/ploidy')
    parser.add_argument('--grch37', action='store_true', help='Also produce GRCh37 (lifted-over) files')
    parser.add_argument('--mnv_max_vaf_diff', type=float, default=0.05,
                        help='Merge adjacent SNVs into MNVs if tumor VAFs differ by at most this (default: 0.05)')
    parser.add_argument('--no_mnv_merge', action='store_true', help='Do not merge adjacent SNVs into MNVs')
    parser.add_argument('--no_expression', action='store_true', help='Do not write RNA-seq gene expression files')
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.join(resolve_release_dir(args.release), 'pcgr_input')
    prepare_samples([args.barcode], args.release, output_dir, args.purity_ploidy, args.grch37,
                    None if args.no_mnv_merge else args.mnv_max_vaf_diff, not args.no_expression)


if __name__ == '__main__':
    main()
