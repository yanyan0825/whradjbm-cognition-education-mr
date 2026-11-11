# -*- coding: utf-8 -*-
"""
复现 Wang et al. 2022 (Int J Obes 46:1487-1492, PMID 35538205) 欧洲人群 MR 结果
  - 暴露: ieu-a-79 WHRadjBMI (Shungin 2015, n=210,082, European, BMI-adjusted)  [原文 n=210,088]
          ieu-a-2  BMI       (Locke 2015,   n=339,224, Mixed)                   [原文 n=322,154]
  - 结局: 本地 Lee 2018 UKB 认知表现 GCST006572 (n=257,841)                     [原文 n=257,828 UKB+COGENT]
  - 流程对齐 TwoSampleMR: tophits(p<5e-8, LD clump r2<0.001) -> harmonise -> IVW(主)/WM/Egger/mode/LOO
  - 离群检测: RadialMR 式 radial IVW (替代 MRPRESSO)
原文目标: WHRadjBMI->cog  beta=-0.07 [-0.12,-0.02] p=0.006
          BMI->cog        beta= 0.00 [-0.07, 0.07] p=0.91
"""
import json, gzip, math, time, sys, os
import numpy as np
from scipy import stats
from mr_rerun import (hdr, api_post, unwrap, get_tophits, harmonize, ratios,
                      ivw, weighted_median, mr_egger, weighted_mode, leave_one_out)
from radial_mr import radial_fit, radial_outliers


def get_tophits_cached(tok, gid):
    cache = "wang2022_tophits_%s.json" % gid
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    hits = get_tophits(tok, gid)
    json.dump(hits, open(cache, "w", encoding="utf-8"))
    return hits

COG_GZ = "ldsc/sumstats/cognition_GCST006572.h.tsv.gz"
EXPOSURES = [
    ("ieu-a-79", "WHRadjBMI (Shungin 2015)"),
    ("ieu-a-2",  "BMI (Locke 2015)"),
]
P_THRESH = 5e-8


def load_outcome(rsids):
    """从本地谐和后的 Lee 2018 认知 sumstats 提取目标 SNP（hm_* 列）"""
    want = set(rsids)
    out = {}
    with gzip.open(COG_GZ, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rs = parts[ix["hm_rsid"]]
            if rs not in want:
                continue
            beta = parts[ix["hm_beta"]]
            se = parts[ix["standard_error"]]
            ea = parts[ix["hm_effect_allele"]]
            nea = parts[ix["hm_other_allele"]]
            eaf = parts[ix["hm_effect_allele_frequency"]]
            if beta in ("", "NA") or se in ("", "NA") or ea in ("", "NA"):
                continue
            out[rs] = dict(rsid=rs, beta=float(beta), se=float(se),
                           ea=ea, nea=nea,
                           eaf=None if eaf in ("", "NA") else float(eaf))
    print("    本地结局命中 %d / %d 个 SNP" % (len(out), len(rsids)))
    return list(out.values())


def fmt_est(m, key="b", se_key="se"):
    return "%.4f [%.4f, %.4f] p=%.2e" % (m[key], m[key] - 1.96 * m[se_key],
                                         m[key] + 1.96 * m[se_key], m["p"])


def run_one(tok, gid, label):
    print("\n" + "=" * 72)
    print("暴露: %s (%s)" % (gid, label))
    print("=" * 72)
    hits = get_tophits_cached(tok, gid)
    rsids = [h["rsid"] for h in hits]
    print("    clumped 显著 SNP: %d" % len(rsids))
    t0 = time.time()
    out_recs = load_outcome(rsids)
    print("    本地结局提取耗时 %.0fs" % (time.time() - t0))
    rows = harmonize(hits, out_recs)
    if len(rows) < 3:
        sys.exit("可用 SNP 不足 3 个")
    r, sr = ratios(rows)
    res_ivw = ivw(r, sr)
    res_wm = weighted_median(r, sr)
    res_eg = mr_egger(np.array([d["b_exp"] for d in rows]),
                      np.array([d["b_out"] for d in rows]),
                      np.array([d["se_out"] for d in rows]))
    res_mo = weighted_mode(r, sr)

    print("\n--- 结果 (n_SNP=%d) ---" % len(rows))
    print("IVW   fixed: %s" % fmt_est(dict(b=res_ivw["b"], se=res_ivw["se_fixed"], p=res_ivw["p_fixed"])))
    print("IVW   random(=TwoSampleMR 默认): %s" % fmt_est(dict(b=res_ivw["b"], se=res_ivw["se_re"], p=res_ivw["p_random"])))
    print("IVW   Cochran Q=%.1f (df=%d) p=%.2e" % (res_ivw["Q"], res_ivw["Q_df"], res_ivw["Q_p"]))
    print("WM          : %s" % fmt_est(res_wm))
    print("Egger       : %.4f [%.4f, %.4f] p=%.2e   intercept=%+.4f p=%.3f" % (
        res_eg["slope"], res_eg["slope"] - 1.96 * res_eg["se"], res_eg["slope"] + 1.96 * res_eg["se"],
        res_eg["p"], res_eg["intercept"], res_eg["p_int"]))
    print("WeightedMode: %s" % fmt_est(res_mo))

    # radial IVW 离群检测（替代 MRPRESSO）
    try:
        cur, removed, m = radial_outliers(list(rows))
        print("Radial IVW  : b=%.4f (FE)  离群点=%d/%d (剔除后 Q_p=%.2e)" % (
            m["theta"], len(removed), len(rows), m["Qp"]))
    except Exception as e:
        print("Radial: failed:", e)

    loo = leave_one_out(rows)
    bs = [x[1] for x in loo]
    print("LOO 范围    : [%.4f, %.4f]" % (min(bs), max(bs)))

    return dict(exposure=gid, label=label, n_snp=len(rows),
                ivw=dict(b=res_ivw["b"], se_fixed=res_ivw["se_fixed"], se_random=res_ivw["se_re"],
                         p_fixed=res_ivw["p_fixed"], p_random=res_ivw["p_random"],
                         Q=res_ivw["Q"], Q_df=res_ivw["Q_df"], Q_p=res_ivw["Q_p"]),
                wm=res_wm, egger=res_eg, mode=res_mo,
                loo=[dict(rsid=x[0], b=x[1], se=x[2], p=x[3]) for x in loo],
                snps=[dict(rsid=d["rsid"], b_exp=d["b_exp"], b_out=d["b_out"],
                           se_out=d["se_out"]) for d in rows])


def main():
    token = open(".opengwas_jwt", encoding="utf-8").read().strip()
    results = {}
    for gid, label in EXPOSURES:
        results[gid] = run_one(token, gid, label)
        time.sleep(1)
    with open("wang2022_replication.json", "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)

    print("\n" + "=" * 72)
    print("与原文对照 (欧洲人群)")
    print("=" * 72)
    w = results["ieu-a-79"]
    b = results["ieu-a-2"]
    print("WHRadjBMI->认知: 原文 -0.07 [-0.12,-0.02] p=0.006 | 本复现 random %.3f p=%.4f" % (
        w["ivw"]["b"], w["ivw"]["p_random"]))
    print("BMI->认知      : 原文  0.00 [-0.07, 0.07] p=0.91  | 本复现 random %.3f p=%.4f" % (
        b["ivw"]["b"], b["ivw"]["p_random"]))


if __name__ == "__main__":
    main()
