# -*- coding: utf-8 -*-
"""
MVMR: WHRadjBMI + Educational attainment -> 认知结局
目的：检验 WHRadjBMI 的认知关联在校正受教育年限后是否消失（即教育多效性解释）
"""
import json, math
import numpy as np
from scipy import stats
import mr_matrix as M

OUT = [("ebi-a-GCST006572", "General cognitive performance"),
       ("ukb-b-2084", "Visual memory"),
       ("ukb-b-15625", "Processing speed"),
       ("ukb-b-20140", "Executive function"),
       ("ukb-b-13823", "Reaction time")]


def mvmr(X, y, se):
    w = 1 / se ** 2
    XtW = X.T * w
    cov = np.linalg.inv(XtW @ X)
    b = cov @ (XtW @ y)
    res = y - X @ b
    k, p = X.shape
    s2 = float((w * res ** 2).sum() / (k - p))
    cov = cov * max(1.0, s2)
    se_b = np.sqrt(np.diag(cov))
    return b, se_b


def cond_f(xj, sej, others):
    """Sanderson 条件F：把 xj 对其他暴露回归，权重用暴露自身 SE"""
    A = np.column_stack([np.ones(len(xj))] + others)
    w = 1.0 / sej ** 2
    AtW = A.T * w
    beta = np.linalg.solve(AtW @ A, AtW @ xj)
    r = xj - A @ beta
    Q = float((w * r ** 2).sum())
    df = len(xj) - A.shape[1]
    return Q / df


def main():
    cache = M.cache_load()
    wer = cache["tophits_ebi-a-GCST90025996"]
    edu = cache["tophits_rev_ieu-a-1239_5e-10"]
    union = {x["rsid"]: x for x in edu}
    for x in wer:
        union.setdefault(x["rsid"], x)
    rs = list(union)
    print("并集 SNP: %d (WHRadj %d, EA %d)" % (len(rs), len(wer), len(edu)), flush=True)

    ea_eff = M.get_assoc(cache, "ieu-a-1239", rs)
    print("EA 效应可用: %d" % len(ea_eff), flush=True)

    res = []
    for oid, oname in OUT:
        om = M.get_assoc(cache, oid, rs)
        X, y, se = [], [], []
        se1, se2 = [], []
        for r in rs:
            w = union[r]; a = ea_eff.get(r); o = om.get(r)
            if not (a and o):
                continue
            e1, n1 = w["ea"].upper(), w["nea"].upper()
            if {e1, n1} in ({"A", "T"}, {"C", "G"}):
                continue
            # EA 对齐到 WHRadjBMI 的效应等位基因
            a1, a2 = (a.get("ea") or "").upper(), (a.get("nea") or "").upper()
            b_ea = a["beta"]
            if a1 == e1 and a2 == n1:
                pass
            elif a1 == n1 and a2 == e1:
                b_ea = -b_ea
            else:
                continue
            # 结局对齐
            o1, o2 = (o.get("ea") or "").upper(), (o.get("nea") or "").upper()
            b_y = o["beta"]
            if o1 == e1 and o2 == n1:
                pass
            elif o1 == n1 and o2 == e1:
                b_y = -b_y
            else:
                continue
            X.append([w["beta"], b_ea]); y.append(b_y); se.append(o["se"])
            se1.append(w["se"]); se2.append(a["se"])
        if len(y) < 10:
            print("  %-30s 可用SNP过少 %d" % (oname, len(y)), flush=True)
            continue
        X = np.array(X); y = np.array(y); se = np.array(se)
        se1 = np.array(se1); se2 = np.array(se2)
        b, sb = mvmr(X, y, se)
        p1 = float(2 * stats.norm.sf(abs(b[0] / sb[0])))
        p2 = float(2 * stats.norm.sf(abs(b[1] / sb[1])))
        cf1 = cond_f(X[:, 0], se1, [X[:, 1]])
        cf2 = cond_f(X[:, 1], se2, [X[:, 0]])
        print("  %-30s k=%3d  WHRadj %+.4f (%.4f) p=%.4f | EA %+.4f (%.4f) p=%.3g | condF %.1f / %.1f"
              % (oname, len(y), b[0], sb[0], p1, b[1], sb[1], p2, cf1, cf2), flush=True)
        res.append(dict(outcome=oid, name=oname, k=len(y),
                        whr_b=round(float(b[0]), 4), whr_se=round(float(sb[0]), 4), whr_p=float("%.4g" % p1),
                        ea_b=round(float(b[1]), 4), ea_se=round(float(sb[1]), 4), ea_p=float("%.4g" % p2),
                        condF_whr=round(cf1, 1), condF_ea=round(cf2, 1)))
    json.dump(res, open("mvmr_whr_edu_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n写出 mvmr_whr_edu_results.json", flush=True)


if __name__ == "__main__":
    main()
