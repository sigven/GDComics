import pandas as pd, numpy as np
import glob
s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob("hrd_scores_raw_*.tsv"))])
m = pd.read_csv("sample_meta.tsv", sep="\t")
d = m.merge(s, on="sample", how="inner")
def brca_status(r):
    hits = []
    for g in ["BRCA1", "BRCA2"]:
        cn = r[f"{g}_cn_status"]
        lof = isinstance(r.somatic_lof_genes, str) and g in r.somatic_lof_genes.split(",")
        loh = cn in ("HEMDEL", "cnLOH")
        if cn == "HOMDEL": hits.append(f"{g}:homdel")
        elif lof and (loh): hits.append(f"{g}:LoF+LOH")
        elif lof: hits.append(f"{g}:LoF(mono)")
        elif cn == "HEMDEL": hits.append(f"{g}:hemdel")
        elif cn == "cnLOH": hits.append(f"{g}:cnLOH")
    return ";".join(hits) if hits else "none"
d["brca12_status"] = d.apply(brca_status, axis=1)
def group(r):
    b = r.brca12_status
    if "homdel" in b or "LoF+LOH" in b: return "BRCA1/2 biallelic (homdel or LoF+LOH)"
    if b == "none": return "No BRCA1/2 alteration"
    return "BRCA1/2 monoallelic (hemdel/cnLOH/LoF)"
d["brca12_group"] = d.apply(group, axis=1)
d["hrd_ge42"] = d.hrd_sum >= 42
d.to_csv("tcga_hrd_scores.tsv", sep="\t", index=False)
print(len(d), d.groupby("project").size().to_dict())
for p, x in d.groupby("project"):
    print(f"\n== {p}: HRD sum quantiles", x.hrd_sum.quantile([.05,.25,.5,.75,.95]).round(0).tolist(), f" >=42: {x.hrd_ge42.mean():.0%}")
    print(x.groupby("brca12_group").hrd_sum.agg(["count","median","mean", lambda v:(v>=42).mean()]).round(2))
    print(x.groupby(["BRCA1_cn_status"]).size().to_dict(), x.groupby(["BRCA2_cn_status"]).size().to_dict())
    print(x.brca12_status.value_counts().to_dict())
pc = d[d.project=="PRAD"]
print("\nPRAD by Gleason:"); print(pc.groupby("gleason_score").hrd_sum.agg(["count","median",lambda v:(v>=42).mean()]).round(2))
bc = d[d.project=="BRCA"]
print("\nBRCA by subtype:"); print(bc.groupby("subtype_selected").hrd_sum.agg(["count","median",lambda v:(v>=42).mean()]).round(2))
print("\nBRCA by ER:"); print(bc.groupby("er_status").hrd_sum.agg(["count","median",lambda v:(v>=42).mean()]).round(2))
print("\nBy WGD:"); print(d.groupby(["project","genome_doubled"]).hrd_sum.agg(["count","median"]))
print(d[["hrd_loh","lst","tai","hrd_sum","fga"]].corr().round(2))
