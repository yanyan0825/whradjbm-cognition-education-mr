# -*- coding: utf-8 -*-
"""Table S2: APOE 19q13 区域剔除前后 IVW 对比（Barton + Pulit 双数据集）"""
import json, gzip, math
import numpy as np
from scipy import stats

APOE_WIN = (44_400_000, 46_500_000)  # GRCh37/38 兼容窗口（chr19 中心 45.4Mb）


def ivw(rows):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    r = by / bx; w = 1 / (sy / np.abs(bx)) ** 2
    b = float((w * r).sum() / w.sum()); se_f = float(math.sqrt(1 / w.sum()))
    k = len(rows); Q = float((w * (r - b) ** 2).sum()); se = se_f * math.sqrt(max(1, Q / (k - 1)))
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))), k=k)


def is_apoe(r):
    return str(r.get("chr")) == "19" and r.get("pos") and APOE_WIN[0] <= int(r["pos"]) <= APOE_WIN[1]


rows_out = []
# A. Barton
H = json.load(open("whradjbmi_harmonised.json"))
names = {"ebi-a-GCST006572": "认知表现(复合)", "ukb-b-15625": "加工速度", "ukb-b-20140": "执行功能(Trail#2)",
         "finn-b-F5_DEMENTIA": "全因痴呆", "finn-b-G6_ALZHEIMER": "AD(FinnGen)", "ieu-a-1239": "教育年数"}
for oid, name in names.items():
    rows = H[oid]
    apoe_snps = [r["rsid"] for r in rows if is_apoe(r)]
    kept = [r for r in rows if not is_apoe(r)]
    m0, m1 = ivw(rows), ivw(kept)
    print("%-14s k=%d->%d  b %+0.4f->%+0.4f  p %.4f->%.4f  (剔除: %s)"
          % (name, m0["k"], m1["k"], m0["b"], m1["b"], m0["p"], m1["p"], ",".join(apoe_snps)))
    rows_out.append(dict(dataset="Barton 2021", outcome=name, k_full=m0["k"], k_no=m1["k"],
                         apoe_removed=",".join(apoe_snps),
                         b_full=round(m0["b"], 5), p_full=round(m0["p"], 5),
                         b_no=round(m1["b"], 5), p_no=round(m1["p"], 5)))

# B. Pulit（主结局）
exp = json.load(open("whradjbmi_instruments.json"))
cache = json.load(open("matrix_cache.json"))["assoc"]
om = cache.get("ebi-a-GCST006572", {})
G = {}
with gzip.open("giant/whradjbmi_pulit2019_hapmap2.txt.gz", "rt") as f:
    hdr = f.readline().split()
    ix = {k: i for i, k in enumerate(hdr)}
    rsids = {e["rsid"] for e in exp}
    for line in f:
        p = line.split()
        snp = p[ix["SNP"]].split(":")[0]
        if snp in rsids:
            G[snp] = dict(ta=p[ix["Tested_Allele"]].upper(), oa=p[ix["Other_Allele"]].upper(),
                          beta=float(p[ix["BETA"]]), se=float(p[ix["SE"]]))
rows_p = []
for e in exp:
    g = G.get(e["rsid"]); o = om.get(e["rsid"])
    if not g or not o:
        continue
    ea, oa = e["ea"].upper(), e["nea"].upper()
    if {ea, oa} in ({"A", "T"}, {"C", "G"}):
        continue
    bx = g["beta"]
    if g["ta"] == ea and g["oa"] == oa:
        pass
    elif g["ta"] == oa and g["oa"] == ea:
        bx = -bx
    else:
        continue
    bo, oo = (o.get("ea") or "").upper(), (o.get("nea") or "").upper()
    by = o["beta"]
    if ea == bo and oa == oo:
        pass
    elif ea == oo and oa == bo:
        by = -by
    else:
        continue
    rows_p.append(dict(rsid=e["rsid"], b_exp=bx, se_exp=g["se"], b_out=by, se_out=o["se"]))
kept_p = [r for r in rows_p if r["rsid"] not in ("rs429358", "rs7412")]
m0, m1 = ivw(rows_p), ivw(kept_p)
print("Pulit 主结局     k=%d->%d  b %+0.4f->%+0.4f  p %.4f->%.4f (剔除 rs429358/rs7412 中命中者)"
      % (m0["k"], m1["k"], m0["b"], m1["b"], m0["p"], m1["p"]))
rows_out.append(dict(dataset="Pulit 2019", outcome="认知表现(复合)", k_full=m0["k"], k_no=m1["k"],
                     apoe_removed=",".join(sorted({"rs429358", "rs7412"} & {r["rsid"] for r in rows_p})),
                     b_full=round(m0["b"], 5), p_full=round(m0["p"], 5),
                     b_no=round(m1["b"], 5), p_no=round(m1["p"], 5)))

json.dump(rows_out, open("apoe_exclusion_results.json", "w"), ensure_ascii=False, indent=1)
print("\n写出 apoe_exclusion_results.json")
