#!/usr/bin/env python

"""
Prepare "PCGR-ready" input files for a set of TCGA samples (see prepare_tcga_pcgr.py).

The sample file is either a plain list of barcodes (one per line) or a TSV with a
header, where barcodes are taken from the first column named 'Sample ID',
'sample_id', 'tumor_sample_barcode' or 'Patient ID' (e.g. a cBioPortal table export).

A combined sample sheet (pcgr_sample_info.tsv) is written to the output directory.

Usage:
  prepare_tcga_pcgr_batch.py ov_brca_cases.txt release46_20260810 --output_dir OV_BRCA
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_tcga_pcgr import prepare_samples, resolve_release_dir, ascat_summary_fname

id_columns = ['Sample ID', 'sample_id', 'tumor_sample_barcode', 'Patient ID', 'patient_id']


def read_barcodes(fname: str) -> list:
    with open(fname, mode='rt') as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip() and not l.startswith('#')]
    header = lines[0].split('\t')
    id_col = next((c for c in id_columns if c in header), None)
    if id_col is None:
        return [l.split('\t')[0].strip() for l in lines]
    idx = header.index(id_col)
    return [l.split('\t')[idx].strip() for l in lines[1:]]


def main():
    parser = argparse.ArgumentParser(description='Prepare PCGR-ready input files for a set of TCGA samples')
    parser.add_argument('sample_file', help='File with TCGA barcodes (plain list or TSV with a sample ID column)')
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

    barcodes = list(dict.fromkeys(read_barcodes(args.sample_file)))
    print(f'Read {len(barcodes)} barcodes from {args.sample_file}')
    output_dir = args.output_dir or os.path.join(resolve_release_dir(args.release), 'pcgr_input')

    sample_info = prepare_samples(barcodes, args.release, output_dir, args.purity_ploidy, args.grch37,
                                  None if args.no_mnv_merge else args.mnv_max_vaf_diff, not args.no_expression)
    if not sample_info:
        sys.exit('ERROR: no samples processed')

    sheet_fname = os.path.join(output_dir, 'pcgr_sample_info.tsv')
    with open(sheet_fname, mode='wt') as fh:
        writer = csv.DictWriter(fh, fieldnames=sample_info[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(sample_info)
    print(f'Processed {len(sample_info)}/{len(barcodes)} samples - sample sheet written to {sheet_fname}')


if __name__ == '__main__':
    main()
