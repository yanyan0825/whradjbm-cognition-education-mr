# -*- coding: utf-8 -*-
"""
反向 MR + Steiger 方向性检验
反向：认知表现 / 受教育年限 -> WHRadjBMI, 腰围, BMI
Steiger：对正向(脂肪分布->认知)与反向(认知->脂肪分布)两套 harmonised 数据分别计算
        r2 = b^2 / (b^2 + se^2 * N)
"""
import json, math, time
import numpy as np
import requests
from scipy import stats
import mr_matrix as M

API = M.API
HD = M.HD
CACHE = M.CACHE

REV_EXPOSURES = [
    ("ebi-a-GCST006572", "Cognitive performance", 5e-8),
    ("ieu-a-1239", "Educational attainment", 5e-10),
]
REV_OUTCOMES = [
    ("ebi-a-GCST90025996", "WHRadjBMI"),
    ("ieu-a-66", "Waist circumference"),
    ("ieu-b-40", "BMI"),
]


def sample_sizes(ids):
    out = {}
    for i in range(0, len(ids), 10):
        d = M.post("/gwasinfo", {"id": ids[i:i + 10]})
        for v in (d if isinstance(d, list) else d.values()):
            ss = v.get("sample_size") or v.get("n")
            try:
                out[v["id"]] = float(ss)
            except Exception:
                pass
    return out


def steiger(rows, n_exp, n_out):
    """返回 r2_exp > r2_out 的比例，以及中位 r2"""
    re_, ro_ = [], []
    for r in rows:
        r2e = r["b_exp"] ** 2 / (r["b_exp"] ** 2 + r["se_exp"] ** 2 * n_exp)
        r2o = r["b_out"] ** 2 / (r["b_out"] ** 2 + r["se_out"] ** 2 * n_out)
        re_.append(r2e); ro_.append(r2o)
    re_ = np.array(re_); ro_ = np.array(ro_)
    return float((re_ > ro_).mean()), float(np.median(re_)), float(np.median(ro_)), re_, ro_


def main():
    cache = M.cache_load()
    nsamp = sample_sizes(["ebi-a-GCST006572", "ieu-a-1239", "ebi-a-GCST90025996",
                          "ieu-a-66", "ieu-b-40", "ukb-b-2084", "ukb-b-15625",
                          "ukb-b-20140", "ukb-b-13823"])
    print("样本量:", {k: int(v) for k, v in nsamp.items()}, flush=True)

    out_rows = []
    print("\n===== 反向 MR =====", flush=True)
    for eid, ename, pthr in REV_EXPOSURES:
        key = "tophits_rev_%s_%g" % (eid, pthr)
        if key in cache:
            exp = cache[key]
        else:
            exp = M.post("/tophits", {"id": eid, "pval": pthr, "clump": 1, "r2": 0.001, "kb": 10000})
            cache[key] = exp
            json.dump(cache, open(CACHE, "w"))
        rs = [x["rsid"] for x in exp]
        F = [(x["beta"] / x["se"]) ** 2 for x in exp]
        print("\n### %s (%s)  %d SNP, F median %.1f min %.1f"
              % (ename, eid, len(exp), np.median(F), min(F)), flush=True)
        print("  %-24s %4s %8s %8s %9s %8s" % ("outcome", "k", "beta", "se", "p", "WM_p"), flush=True)
        for oid, oname in REV_OUTCOMES:
            om = M.get_assoc(cache, oid, rs)
            rows = M.harmonise(exp, om)
            if len(rows) < 8:
                print("  %-24s 可用SNP过少(%d)" % (oname, len(rows)), flush=True)
                continue
            m = M.ivw(rows); w = M.wm(rows); eg = M.egger(rows)
            frac, mre, mro, _, _ = steiger(rows, nsamp[eid], nsamp[oid])
            print("  %-24s %4d %+8.4f %8.4f %9.5f %8.4f   Steiger方向正确 %.0f%% (r2 %.2e vs %.2e)"
                  % (oname, m["k"], m["b"], m["se"], m["p"], w["p"], frac * 100, mre, mro), flush=True)
            out_rows.append(dict(direction="reverse", exposure=eid, exposure_name=ename,
                                 outcome=oid, outcome_name=oname, k=m["k"],
                                 b=round(m["b"], 5), se=round(m["se"], 5), p=round(m["p"], 5),
                                 lo=round(m["lo"], 5), hi=round(m["hi"], 5),
                                 wm_b=round(w["b"], 5), wm_p=round(w["p"], 5),
                                 egger_ic_p=round(eg["ic_p"], 3),
                                 Q=round(m["Q"], 1), Qp=float("%.3g" % m["Qp"]),
                                 steiger_frac=round(frac, 3),
                                 mdes=round(2.8 * m["se"], 4)))

    # ---- Steiger 正向复核：WHRadjBMI -> 5 认知结局 ----
    print("\n===== Steiger 正向复核 (WHRadjBMI -> 认知) =====", flush=True)
    H = json.load(open("whradjbmi_harmonised.json"))
    for oid, oname in [("ebi-a-GCST006572", "General cognitive performance"),
                       ("ukb-b-2084", "Visual memory"), ("ukb-b-15625", "Processing speed"),
                       ("ukb-b-20140", "Executive function"), ("ukb-b-13823", "Reaction time")]:
        rows = H.get(oid)
        if not rows:
            continue
        frac, mre, mro, _, _ = steiger(rows, nsamp["ebi-a-GCST90025996"], nsamp[oid])
        print("  %-30s k=%3d  Steiger方向正确 %.0f%%   median r2_exp=%.2e  r2_out=%.2e"
              % (oname, len(rows), frac * 100, mre, mro), flush=True)
        out_rows.append(dict(direction="forward_steiger", exposure="ebi-a-GCST90025996",
                             exposure_name="WHRadjBMI", outcome=oid, outcome_name=oname,
                             k=len(rows), b="", se="", p="", lo="", hi="",
                             wm_b="", wm_p="", egger_ic_p="", Q="", Qp="",
                             steiger_frac=round(frac, 3), mdes=""))

    with open("mr_reverse_results.json", "w", encoding="utf-8") as f:
        json.dump(out_rows, f, ensure_ascii=False, indent=1)
    print("\n写出 mr_reverse_results.json (%d 行)" % len(out_rows), flush=True)


if __name__ == "__main__":
    main()
