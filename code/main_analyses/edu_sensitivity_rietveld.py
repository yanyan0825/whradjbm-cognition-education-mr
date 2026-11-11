# -*- coding: utf-8 -*-
"""
非 UKB 教育 GWAS 敏感性分析（Rietveld 2014 PNAS, GCST002598, n=106,736, 无 UK Biobank）
对 241 个 WHRadjBMI 工具变量：
  1) 取 Rietveld 教育效应 z 值（对齐等位基因）
  2) 重跑教育分层 IVW（全部 / 剔除|z|>3 / 剔除|z|>2 / 仅|z|<1），5 个认知结局
  3) 重跑 MVMR（WHRadjBMI + Rietveld 教育 -> 认知），5 个认知结局
输出：edu_rietveld_strat.json + mvmr_rietveld_results.json
"""
import json, math, gzip
import numpy as np
from scipy import stats

cache = json.load(open("matrix_cache.json"))
mc = json.load(open("whradjbmi_assoc_cache.json"))
exp = cache["tophits_ebi-a-GCST90025996"]          # 241 WHRadjBMI tophits
lee = mc["ieu-a-1239"]                             # Lee 2018 教育效应（主分析）
OUT = [("ebi-a-GCST006572", "General cognitive performance"),
       ("ukb-b-2084", "Visual memory"),
       ("ukb-b-15625", "Processing speed"),
       ("ukb-b-20140", "Executive function"),
       ("ukb-b-13823", "Reaction time")]
RIET = "ldsc/sumstats/education_rietveld2014_GCST002598.h.tsv.gz"


def load_rietveld():
    """流式读取 harmonised 文件，取 241 个 rsid 的教育效应"""
    want = {e["rsid"] for e in exp}
    got = {}
    with gzip.open(RIET, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rs = p[ix["hm_rsid"]]
            if rs not in want or rs in got:
                continue
            b, pv = p[ix["hm_beta"]], p[ix["p_value"]]
            if not b or not pv or b in ("NA", "") or pv in ("NA", ""):
                continue
            # harmonised 文件缺 SE 列，由 beta 与 p 反推：se = |beta| / z
            pval = float(pv)
            if pval <= 0 or pval >= 1:
                continue
            z = abs(stats.norm.ppf(pval / 2.0))
            se = abs(float(b)) / z if z > 0 else None
            if se is None or se == 0:
                continue
            got[rs] = dict(ea=(p[ix["hm_effect_allele"]] or "").upper(),
                           nea=(p[ix["hm_other_allele"]] or "").upper(),
                           beta=float(b), se=float(se))
    return got


def align_beta(ref_ea, ref_nea, cand):
    """把候选效应按参考 (ea, nea) 对齐；无法对齐返回 None"""
    a1, a2 = cand["ea"], cand["nea"]
    if a1 == ref_ea and a2 == ref_nea:
        return cand["beta"], cand["se"]
    if a1 == ref_nea and a2 == ref_ea:
        return -cand["beta"], cand["se"]
    return None, None


def ivw(rows):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    r = by / bx; w = 1 / (sy / np.abs(bx)) ** 2
    b = float((w * r).sum() / w.sum()); se_f = float(math.sqrt(1 / w.sum()))
    k = len(rows); Q = float((w * (r - b) ** 2).sum()); se = se_f * math.sqrt(max(1, Q / (k - 1)))
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))), k=k)


def mvmr(X, y, se):
    w = 1 / se ** 2
    XtW = X.T * w
    cov = np.linalg.inv(XtW @ X)
    b = cov @ (XtW @ y)
    res = y - X @ b
    k, p = X.shape
    s2 = float((w * res ** 2).sum() / (k - p))
    cov = cov * max(1.0, s2)
    return b, np.sqrt(np.diag(cov))


def cond_f(xj, sej, others):
    A = np.column_stack([np.ones(len(xj))] + others)
    w = 1.0 / sej ** 2
    AtW = A.T * w
    beta = np.linalg.solve(AtW @ A, AtW @ xj)
    r = xj - A @ beta
    Q = float((w * r ** 2).sum())
    return Q / (len(xj) - A.shape[1])


def build(oid, riet):
    om = cache["assoc"].get(oid, {})
    rows = []
    for e in exp:
        o = om.get(e["rsid"])
        if not o:
            continue
        be, oe = e["ea"].upper(), e["nea"].upper()
        # 结局对齐
        bo, oo = (o.get("ea") or "").upper(), (o.get("nea") or "").upper()
        by = o["beta"]
        if be == bo and oe == oo:
            pass
        elif be == oo and oe == bo:
            by = -by
        else:
            continue
        if {be, oe} in ({"A", "T"}, {"C", "G"}):
            eaf = e.get("eaf")
            if eaf is None or 0.42 < eaf < 0.58:
                continue
        # Rietveld 教育效应（对齐 WHRadjBMI 效应等位基因）
        b_r, se_r = (None, None)
        if e["rsid"] in riet:
            b_r, se_r = align_beta(be, oe, riet[e["rsid"]])
        z_r = b_r / se_r if b_r is not None and se_r else None
        # Lee 教育效应（对齐同侧）
        z_l = None
        if e["rsid"] in lee:
            z_l = lee[e["rsid"]]["beta"] / lee[e["rsid"]]["se"]
        rows.append(dict(rsid=e["rsid"], b_exp=float(e["beta"]), se_exp=float(e["se"]),
                         b_out=float(by), se_out=float(o["se"]), z_riet=z_r, z_lee=z_l))
    return rows


def main():
    riet = load_rietveld()
    print("Rietveld 命中工具变量: %d / 241" % len(riet), flush=True)

    rows0 = build("ebi-a-GCST006572", riet)
    have = [(r["z_riet"], r["z_lee"]) for r in rows0 if r["z_riet"] is not None and r["z_lee"] is not None]
    if len(have) > 10:
        a, b = np.array(have).T
        rho = float(stats.pearsonr(a, b)[0])
        print("z 分布: |z|<=1 %.0f%%, <=2 %.0f%%, <=3 %.0f%%; 与 Lee z 相关 r=%.3f"
              % (100*np.mean(abs(a)<=1), 100*np.mean(abs(a)<=2), 100*np.mean(abs(a)<=3), rho), flush=True)
    n_noz = sum(1 for r in rows0 if r["z_riet"] is None)
    print("无教育信息的工具变量: %d" % n_noz, flush=True)

    # --- 分层 ---
    strat_out = []
    print("\n== Rietveld 教育分层 IVW ==")
    for oid, name in OUT:
        rows = build(oid, riet)
        rec = dict(outcome=oid, name=name)
        line = [name]
        for lab, flt in [("all", None), ("z3", 3), ("z2", 2), ("z1", 1)]:
            if flt is None:
                sub = rows
            elif lab == "z1":
                sub = [r for r in rows if r["z_riet"] is not None and abs(r["z_riet"]) < 1]
            else:
                sub = [r for r in rows if r["z_riet"] is not None and abs(r["z_riet"]) <= flt]
            if len(sub) < 10:
                line.append("k<10"); rec[lab] = None; continue
            m = ivw(sub)
            line.append("%+.3f (p=%.3f, k=%d)" % (m["b"], m["p"], m["k"]))
            rec[lab] = dict(b=round(m["b"], 4), se=round(m["se"], 4), p=round(m["p"], 4), k=m["k"])
        print("  " + " | ".join(line))
        strat_out.append(rec)
    json.dump(strat_out, open("edu_rietveld_strat.json", "w"), ensure_ascii=False, indent=1)

    # --- MVMR ---
    print("\n== Rietveld MVMR ==")
    mvmr_out = []
    for oid, name in OUT:
        om = cache["assoc"].get(oid, {})
        X, y, se, se1, se2 = [], [], [], [], []
        for e in exp:
            o = om.get(e["rsid"])
            a = riet.get(e["rsid"])
            if not (o and a):
                continue
            be, ne = e["ea"].upper(), e["nea"].upper()
            if {be, ne} in ({"A", "T"}, {"C", "G"}):
                continue
            b_r, se_r = align_beta(be, ne, a)
            if b_r is None:
                continue
            bo, oo = (o.get("ea") or "").upper(), (o.get("nea") or "").upper()
            b_y = o["beta"]
            if bo == be and oo == ne:
                pass
            elif bo == ne and oo == be:
                b_y = -b_y
            else:
                continue
            X.append([e["beta"], b_r]); y.append(b_y); se.append(o["se"])
            se1.append(e["se"]); se2.append(se_r)
        if len(y) < 10:
            print("  %-30s 可用SNP过少 %d" % (name, len(y)))
            continue
        X = np.array(X); y = np.array(y); se = np.array(se)
        se1 = np.array(se1); se2 = np.array(se2)
        b, sb = mvmr(X, y, se)
        p1 = float(2 * stats.norm.sf(abs(b[0] / sb[0])))
        p2 = float(2 * stats.norm.sf(abs(b[1] / sb[1])))
        cf1 = cond_f(X[:, 0], se1, [X[:, 1]])
        cf2 = cond_f(X[:, 1], se2, [X[:, 0]])
        print("  %-30s k=%3d  WHRadj %+.4f (%.4f) p=%.4f | Rietveld-edu %+.4f (%.4f) p=%.3g | condF %.1f / %.1f"
              % (name, len(y), b[0], sb[0], p1, b[1], sb[1], p2, cf1, cf2))
        mvmr_out.append(dict(outcome=oid, name=name, k=len(y),
                             whr_b=round(float(b[0]), 4), whr_se=round(float(sb[0]), 4), whr_p=float("%.4g" % p1),
                             edu_b=round(float(b[1]), 4), edu_se=round(float(sb[1]), 4), edu_p=float("%.4g" % p2),
                             condF_whr=round(cf1, 1), condF_edu=round(cf2, 1)))
    json.dump(mvmr_out, open("mvmr_rietveld_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n写出 edu_rietveld_strat.json / mvmr_rietveld_results.json")


if __name__ == "__main__":
    main()
