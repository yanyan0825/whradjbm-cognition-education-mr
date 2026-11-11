# -*- coding: utf-8 -*-
"""
Radial MR (Bowden et al. 2018; 对齐 RadialMR 包 ivw_radial 实现)
  - ratio 尺度加权回归: r_i = by_i/bx_i, W_i = (bx_i/se_out_i)^2, theta = sum(W r)/sum(W)
  - SE(fixed, first-order) = sqrt(1/sum(W))
  - 每位点 Q_j = W_j (r_j - theta)^2 ~ chi2(1); 离群判定: p < alpha (包默认名义 0.05, 同时报告 Bonferroni)
  - 迭代剔除最极端离群点后重拟合
数据源:
  A. Barton 2021 暴露 (whradjbmi_harmonised.json)
  B. GIANT Pulit 2019 暴露复现 (本地文件 + matrix_cache.json)
"""
import json, gzip, math
import numpy as np
from scipy import stats

OUT_LABEL = {
    "ebi-a-GCST006572": "认知表现(复合)",
    "ukb-b-2084": "视觉记忆",
    "ukb-b-15625": "加工速度",
    "ukb-b-20140": "执行功能(Trail#2)",
    "ukb-b-13823": "反应时",
    "ieu-a-1239": "教育年数(阴性对照)",
}


def radial_fit(rows):
    """ratio 尺度 radial IVW（对齐 RadialMR::ivw_radial, first-order weights）"""
    bx = np.array([r["b_exp"] for r in rows])
    by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    W = (bx / sy) ** 2                      # = 1/se_ratio^2
    r = by / bx
    theta = float((W * r).sum() / W.sum())
    se = float(math.sqrt(1.0 / W.sum()))    # lm 摘要等价 (FE, first-order)
    Qj = W * (r - theta) ** 2
    Q = float(Qj.sum())
    k = len(rows)
    Qp = float(stats.chi2.sf(Q, k - 1))
    return dict(theta=theta, se=se, Qj=Qj, Q=Q, Qp=Qp, k=k,
                p=float(2 * stats.norm.sf(abs(theta / se))),
                W=W, r=r, rsid=[x["rsid"] for x in rows])


def radial_outliers(rows, alpha=0.05):
    """迭代剔除：每轮挑 Qj 贡献最大且 p<alpha 的位点（与包一致用名义 alpha；剔除上限 10%）"""
    cur = list(rows)
    removed = []
    max_rm = max(5, int(0.10 * len(rows)))
    while True:
        m = radial_fit(cur)
        pj = stats.chi2.sf(m["Qj"], 1)
        bad = np.where(pj < alpha)[0]
        if len(bad) == 0 or len(cur) <= 5 or len(removed) >= max_rm:
            return cur, removed, m
        worst = int(bad[np.argmax(m["Qj"][bad])])
        removed.append(cur[worst]["rsid"])
        cur.pop(worst)


def report(rows, label):
    m0 = radial_fit(rows)
    cur, rm, m1 = radial_outliers(rows)
    print("\n[%s]  k=%d" % (label, m0["k"]))
    print("  Radial IVW 全样本 : b=%+.5f se=%.5f p=%.5f  Q=%.1f p=%.3g"
          % (m0["theta"], m0["se"], m0["p"], m0["Q"], m0["Qp"]))
    if rm:
        print("  剔除离群 %d 个: %s" % (len(rm), ",".join(rm)))
        print("  Radial IVW 剔除后 : b=%+.5f se=%.5f p=%.5f (95%%CI %+.4f..%+.4f)  Q=%.1f p=%.3f k=%d"
              % (m1["theta"], m1["se"], m1["p"],
                 m1["theta"] - 1.96 * m1["se"], m1["theta"] + 1.96 * m1["se"],
                 m1["Q"], m1["Qp"], m1["k"]))
    else:
        print("  无离群点 (chi2(1) 名义 p<0.05)")
    return dict(label=label, k0=m0["k"], b0=m0["theta"], se0=m0["se"], p0=m0["p"],
                Q0=m0["Q"], Qp0=m0["Qp"], removed=rm,
                k1=len(cur), b1=m1["theta"], se1=m1["se"], p1=m1["p"], Q1=m1["Q"], Qp1=m1["Qp"])


def main():
    results = []
    # A. Barton 2021 暴露
    H = json.load(open("whradjbmi_harmonised.json"))
    for oid in ["ebi-a-GCST006572", "ukb-b-2084", "ukb-b-15625", "ukb-b-20140", "ukb-b-13823", "ieu-a-1239"]:
        if oid in H:
            results.append(report(H[oid], "Barton2021 -> " + OUT_LABEL[oid]))

    # B. GIANT Pulit 2019 暴露复现（主结局）
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
    results.append(report(rows_p, "Pulit2019 -> 认知表现(复合)"))

    json.dump(results, open("radial_mr_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n写出 radial_mr_results.json")


if __name__ == "__main__":
    main()
