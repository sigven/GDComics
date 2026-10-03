import sys, time, logging
sys.path.insert(0, "/Users/sigven/project_data/packages/package__pcgr/pcgr")
import pandas as pd
from multiprocessing import Pool
from pcgr import hrd

CHROMSIZES = "/Users/sigven/project_data/data/data__pcgrdb/bundle_output/20261001/data/grch38/chromsize.grch38.tsv"
GENES = {"BRCA1": ("17", 43044295, 43170245), "BRCA2": ("13", 32315086, 32400268)}  # GRCh38
log = logging.getLogger("hrd"); log.setLevel(logging.ERROR)
arms = hrd.scarhrd_chrom_arms(hrd.read_chrom_arms(CHROMSIZES), "grch38")   # scarHRD centromere coordinates (identical scores to scarHRD)

def gene_status(seg, chrom, gs, ge):
    s = seg[(seg.Chromosome.astype(str) == chrom) & (seg.Start <= ge) & (seg.End >= gs)].copy()
    if s.empty:
        return "NA", None, None
    s["ov"] = s[["End"]].min(axis=1).clip(upper=ge) - s["Start"].clip(lower=gs)
    r = s.sort_values("ov", ascending=False).iloc[0]
    tot, minor = int(r.nMajor + r.nMinor), int(r.nMinor)
    if tot == 0: st = "HOMDEL"
    elif minor == 0 and tot == 1: st = "HEMDEL"
    elif minor == 0: st = "cnLOH"
    else: st = "intact"
    return st, tot, minor

def one(args):
    sample, seg = args
    seg = seg.drop(columns="sample").reset_index(drop=True)
    row = {"sample": sample}
    try:
        row["hrd_loh"] = hrd.calc_hrd_loh(seg, logger=log)
        row["lst"] = hrd.calc_lst(seg, arms, logger=log)
        row["tai"] = hrd.calc_tai(seg, arms, logger=log)
        row["hrd_sum"] = row["hrd_loh"] + row["lst"] + row["tai"]
        row["fga"] = hrd.calc_fraction_genome_altered(seg, chrom_arms=arms, logger=log)
        w = hrd.calc_genome_doubling(seg, chrom_arms=arms, logger=log)
        row["wgd_fraction"], row["genome_doubled"] = (w["wgd_fraction"], w["genome_doubled"]) if w else (None, None)
    except Exception as e:
        row["error"] = repr(e)
    for g, (c, a, b) in GENES.items():
        st, tot, mi = gene_status(seg, c, a, b)
        row[f"{g}_cn_status"], row[f"{g}_total_cn"], row[f"{g}_minor_cn"] = st, tot, mi
    return row

if __name__ == "__main__":
    import os
    projects = sys.argv[1].split(",") if len(sys.argv) > 1 else ["BRCA", "OV", "PRAD", "PAAD"]
    for p in projects:
        out_fname = f"hrd_scores_raw_{p}.tsv"
        if os.path.exists(out_fname):
            print(p, "cached:", out_fname); continue
        seg = pd.read_csv(f"segments_{p}.tsv", sep="\t", dtype={"Chromosome": str})
        groups = list(seg.groupby("sample"))
        t = time.time()
        with Pool(8) as pool:
            rows = pool.map(one, groups, chunksize=8)
        out = pd.DataFrame(rows)
        out.to_csv(out_fname, sep="\t", index=False)
        print(p, len(groups), f"{time.time()-t:.0f}s", out.get("error", pd.Series()).notna().sum(), "errors", flush=True)
