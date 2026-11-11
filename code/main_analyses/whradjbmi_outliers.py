# -*- coding: utf-8 -*-
"""WHRadjBMI -> 认知：离群位点校正 + 逐位点影响力诊断"""
import json, math
import numpy as np
from scipy import stats

H = json.load(open("whradjbmi_harmonised.json"))


def ivw(rows_return_full=False, rows=None):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    r = by / bx; se_r = sy / np.abs(bx); w = 1 / se_r ** 2
    b = float((w * r).sum() / w.sum()); se_f = float(math.sqrt(1 / w.sum()))
    k = len(rows); Q = float((w * (r - b) ** 2).sum()); phi = max(1.0, Q / (k - 1))
    se = se_f * math.sqrt(phi)
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))), Q=Q, k=k,
                Qp=float(stats.chi2.sf(Q, k - 1)), lo=b - 1.96 * se, hi=b + 1.96 * se)


def prune_outliers(rows, alpha=0.05, max_frac=0.10, label=""):
    cur = list(rows)
    removed = []
    while True:
        m = ivw(rows=cur)
        if m["Qp"] >= alpha or len(cur) <= 5 or len(removed) >= int(max_frac * len(rows)):
            break
        # 逐位点贡献
        contrib = []
        for i in range(len(cur)):
            mm = ivw(rows=cur[:i] + cur[i + 1:])
            contrib.append((mm["Q"], cur[i]["rsid"]))
        contrib.sort()
        worst = contrib[0][1]
        removed.append(worst)
        cur = [r for r in cur if r["rsid"] != worst]
    return cur, removed


print("=== 离群位点迭代剔除 (MR-PRESSO 式) ===")
for oid in ["ebi-a-GCST006572", "ukb-b-20140", "ukb-b-15625"]:
    rows = H[oid]
    cur, rm = prune_outliers(rows, label=oid)
    m0 = ivw(rows=rows); m1 = ivw(rows=cur)
    print("%-18s 原: b=%+.4f p=%.4f Q=%.1f (p=%.1e, k=%d)" % (oid, m0["b"], m0["p"], m0["Q"], m0["Qp"], m0["k"]))
    print("%-18s 剔除 %d 个离群: b=%+.4f (%.3f,%.3f) p=%.4f Q=%.1f (p=%.3f, k=%d)  剔除=%s" %
          ("", len(rm), m1["b"], m1["lo"], m1["hi"], m1["p"], m1["Q"], m1["Qp"], m1["k"], ",".join(rm[:8])))

# --- 单 SNP 杠杆诊断：每个 SNP 的 Cook 式影响 ---
print("\n=== 影响力最大的 8 个位点 (认知表现) ===")
rows = H["ebi-a-GCST006572"]
full = ivw(rows=rows)
inf = []
for i in range(len(rows)):
    m = ivw(rows=rows[:i] + rows[i + 1:])
    inf.append((abs(m["b"] - full["b"]) / full["se"], rows[i]["rsid"], round(rows[i]["b_exp"], 4), round(rows[i]["b_out"] / rows[i]["b_exp"], 2)))
inf.sort(reverse=True)
for x in inf[:8]:
    print("   %-12s b_exp=%+.4f ratio=%+.2f  影响=%.2f SE" % (x[1], x[2], x[3], x[0]))
