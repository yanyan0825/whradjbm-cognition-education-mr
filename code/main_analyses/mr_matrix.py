# -*- coding: utf-8 -*-
"""
主分析矩阵：3 个暴露 × (1 主要结局 + 4 认知域)
暴露：腰围 WC (ieu-a-66) / WHRadjBMI (ebi-a-GCST90025996) / BMI (ieu-b-40)
结局：一般认知表现 / 视觉记忆 / 加工速度 / 执行功能 / 反应时
输出：mr_matrix_results.csv + 控制台表格
"""
import json, math, os, time, csv
import numpy as np
import requests
from scipy import stats

API = "https://api.opengwas.io/api"
TOK = open(".opengwas_jwt").read().strip()
HD = {"Authorization": "Bearer " + TOK, "Content-Type": "application/json"}
CACHE = "matrix_cache.json"
BATCH = 45

EXPOSURES = [("ieu-a-66", "Waist circumference"), ("ebi-a-GCST90025996", "WHRadjBMI"), ("ieu-b-40", "BMI")]
OUTCOMES = [
    ("ebi-a-GCST006572", "General cognitive performance", "primary"),
    ("ukb-b-2084", "Visual memory (correct matches)", "domain"),
    ("ukb-b-15625", "Processing speed (symbol digit)", "domain"),
    ("ukb-b-20140", "Executive function (trail #2 time)", "domain"),
    ("ukb-b-13823", "Reaction time", "domain"),
]
HIGHER_IS_WORSE = {"ukb-b-20140", "ukb-b-13823"}


def cache_load():
    return json.load(open(CACHE)) if os.path.exists(CACHE) else {"assoc": {}}


def post(path, body, retry=6):
    for k in range(retry):
        try:
            r = requests.post(API + path, headers=HD, json=body, timeout=150)
            if r.status_code == 200:
                return r.json()
            print("   HTTP", r.status_code, "retry", k, flush=True)
        except Exception as e:
            print("   net", str(e)[:60], flush=True)
        time.sleep(6 * (k + 1))
    raise RuntimeError(path)


def get_assoc(cache, gid, rsids):
    key = gid
    got = cache["assoc"].get(key, {})
    missing = [r for r in rsids if r not in got]
    for i in range(0, len(missing), BATCH):
        chunk = missing[i:i + BATCH]
        rows = post("/associations", {"variant": chunk, "id": [gid]})
        for x in rows:
            got[x["rsid"]] = x
        time.sleep(0.7)
    if missing:
        cache["assoc"][key] = got
        json.dump(cache, open(CACHE, "w"))
    return got


def tophits(cache, gid):
    k = "tophits_" + gid
    if k in cache:
        return cache[k]
    h = post("/tophits", {"id": gid, "pval": 5e-8, "clump": 1, "r2": 0.001, "kb": 10000})
    cache[k] = h
    json.dump(cache, open(CACHE, "w"))
    return h


def harmonise(exp, outmap):
    rows = []
    for e in exp:
        o = outmap.get(e["rsid"])
        if not o:
            continue
        be, oe = e["ea"].upper(), e["nea"].upper()
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
        rows.append(dict(rsid=e["rsid"], b_exp=float(e["beta"]), se_exp=float(e["se"]),
                         b_out=float(by), se_out=float(o["se"]), p_out=float(o["p"]),
                         n_out=o.get("n")))
    return rows


def ivw(rows):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    r = by / bx; w = 1.0 / (sy / np.abs(bx)) ** 2
    b = float((w * r).sum() / w.sum()); se_f = float(math.sqrt(1 / w.sum()))
    k = len(rows); Q = float((w * (r - b) ** 2).sum()); phi = max(1.0, Q / (k - 1))
    se = se_f * math.sqrt(phi)
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))), lo=b - 1.96 * se, hi=b + 1.96 * se,
                Q=Q, Qdf=k - 1, Qp=float(stats.chi2.sf(Q, k - 1)), k=k)


def wm(rows, nboot=2000, seed=11):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows]); w = (bx / sy) ** 2; r = by / bx

    def est(rr, ww):
        o = np.argsort(rr); rr2, ww2 = rr[o], ww[o]
        cw = np.cumsum(ww2) - 0.5 * ww2; cw /= cw[-1]
        i = np.where(cw < 0.5)[0][-1]
        return float(rr2[i] + (rr2[i + 1] - rr2[i]) * (0.5 - cw[i]) / (cw[i + 1] - cw[i]))
    b = est(r, w)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(nboot):
        i = rng.integers(0, len(r), len(r))
        try:
            bs.append(est(r[i], w[i]))
        except Exception:
            pass
    se = float(np.std(bs, ddof=1))
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))))


def egger(rows):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows]); w = 1 / sy ** 2
    X = np.column_stack([np.ones(len(bx)), bx]); XtW = X.T * w
    cov = np.linalg.inv(XtW @ X); beta = cov @ (XtW @ by)
    res = by - X @ beta; k = len(bx)
    s2 = float((w * res ** 2).sum() / (k - 2)); cov = cov * max(1.0, s2)
    se = np.sqrt(np.diag(cov))
    return dict(b=float(beta[1]), p=float(2 * stats.norm.sf(abs(beta[1] / se[1]))),
                ic=float(beta[0]), ic_p=float(2 * stats.norm.sf(abs(beta[0] / se[0]))))


def trim_ivw(rows, cut=5.0):
    r = np.array([x["b_out"] / x["b_exp"] for x in rows])
    keep = [rows[i] for i in range(len(rows)) if abs(r[i]) <= cut]
    if len(keep) == len(rows) or len(keep) < 5:
        return None, 0
    return ivw(keep), len(rows) - len(keep)


def main():
    cache = cache_load()
    results = []
    for eid, ename in EXPOSURES:
        exp = tophits(cache, eid)
        rs = [x["rsid"] for x in exp]
        F = [(x["beta"] / x["se"]) ** 2 for x in exp]
        print("\n### %s (%s)  %d SNPs, F median %.1f" % (ename, eid, len(exp), np.median(F)), flush=True)
        print("  %-38s %4s %8s %8s %8s %8s" % ("outcome", "k", "beta", "se", "p", "beta_trim"))
        for oid, oname, tier in OUTCOMES:
            om = get_assoc(cache, oid, rs)
            rows = harmonise(exp, om)
            if len(rows) < 10:
                print("  %-38s 可用SNP过少(%d)" % (oname[:38], len(rows))); continue
            m = ivw(rows); mt, nt = trim_ivw(rows)
            m_wm = wm(rows); m_eg = egger(rows)
            b_trim = round(mt["b"], 4) if mt else ""
            p_trim = round(mt["p"], 4) if mt else ""
            print("  %-38s %4d %+8.4f %8.4f %8.4f %8s" % (oname[:38], m["k"], m["b"], m["se"], m["p"], b_trim), flush=True)
            results.append(dict(exposure=eid, exposure_name=ename, outcome=oid, outcome_name=oname, tier=tier,
                                higher_is_worse=int(oid in HIGHER_IS_WORSE), k=m["k"],
                                b=round(m["b"], 5), se=round(m["se"], 5), p=round(m["p"], 5),
                                lo=round(m["lo"], 5), hi=round(m["hi"], 5),
                                wm_b=round(m_wm["b"], 5), wm_p=round(m_wm["p"], 5),
                                egger_ic_p=round(m_eg["ic_p"], 3),
                                Q=round(m["Q"], 1), Qp=float("%.3g" % m["Qp"]),
                                trim_n=nt, trim_b=b_trim, trim_p=p_trim,
                                mdes=round(2.8 * m["se"], 4)))
    with open("mr_matrix_results.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader(); w.writerows(results)
    print("\n写出 mr_matrix_results.csv  (%d 行)" % len(results))


if __name__ == "__main__":
    main()
