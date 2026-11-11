# -*- coding: utf-8 -*-
"""
独立重跑: 腰围 -> 全因痴呆 两样本MR (不依赖任何旧结果)
数据源 (OpenGWAS API):
  暴露: ieu-a-65 (GIANT 腰围, Shungin 2015)
  结局: finn-b-KRA_PSY_DEMENTIA_EXMORE (FinnGen R10 全因痴呆, 排除更多对照)
流程: 全基因组显著SNP(P<5e-8) -> 服务端LD clumping -> 效应值提取 -> harmonization
      -> IVW(固定/随机) / 加权中位数 / MR-Egger / 加权众数 / Cochran's Q / 留一法
用法: python mr_rerun.py --token <JWT>   (或环境变量 OPENGWAS_JWT)
"""
import argparse, json, math, random, sys, time
import requests
import numpy as np
from scipy import stats

API = "https://api.opengwas.io/api"
EXPOSURE = "ieu-a-65"
OUTCOME = "finn-b-KRA_PSY_DEMENTIA_EXMORE"
P_THRESH = 5e-8

rng = random.Random(20261109)

def hdr(tok):
    return {"Authorization": "Bearer " + tok, "Content-Type": "application/json"}

def api_get(tok, path):
    r = requests.get(API + path, headers=hdr(tok), timeout=60)
    if r.status_code == 429:
        print("  限流, 等待30s..."); time.sleep(30)
        return api_get(tok, path)
    r.raise_for_status()
    return r.json()

def api_post(tok, path, body, retry=3):
    for i in range(retry):
        r = requests.post(API + path, headers=hdr(tok), data=json.dumps(body), timeout=120)
        if r.status_code == 429:
            print("  限流, 等待30s..."); time.sleep(30); continue
        if r.status_code >= 400:
            print("  HTTP %d: %s" % (r.status_code, r.text[:300]))
        r.raise_for_status()
        return r.json()
    raise RuntimeError("API重试耗尽: " + path)

def unwrap(resp):
    """兼容 API 返回 {data:...} / {associations:...} / 直接列表"""
    if isinstance(resp, list):
        return resp
    for k in ("data", "associations", "results", "tophits"):
        if k in resp:
            return resp[k]
    return resp

# ---------- 数据获取 ----------
def get_tophits(tok, gid):
    print("[1] 获取 %s 全基因组显著SNP并LD-clumping (P<%g)..." % (gid, P_THRESH))
    resp = api_post(tok, "/tophits", {"id": gid, "pval": P_THRESH, "clump": 1})
    hits = unwrap(resp)
    if hits and isinstance(hits, dict):
        hits = list(hits.values())[0]
    print("    独立显著SNP数: %d" % len(hits))
    return hits

def get_outcome(tok, gid, rsids):
    print("[2] 提取结局 %s 在 %d 个SNP上的效应值..." % (gid, len(rsids)))
    out = []
    for i in range(0, len(rsids), 100):
        batch = rsids[i:i+100]
        resp = api_post(tok, "/associations", {"variant": batch, "id": [gid]})
        unw = unwrap(resp)
        if isinstance(unw, dict):
            # {rsid: [records]} 形态
            for rsid, recs in unw.items():
                for rec in recs:
                    if rec.get("id") == gid:
                        rec["_rsid"] = rsid; out.append(rec)
        elif isinstance(unw, list):
            for rec in unw:
                if rec.get("id") == gid:
                    out.append(rec)
        time.sleep(0.6)
    print("    拿到结局记录: %d" % len(out))
    return out

# ---------- Harmonization ----------
def harmonize(hits, out_recs):
    print("[3] 等位基因 harmonization...")
    omap = {r.get("_rsid") or r.get("rsid"): r for r in out_recs}
    rows, dropped = [], []
    for h in hits:
        rs = h.get("rsid") or h.get("_rsid")
        o = omap.get(rs)
        if not o or o.get("beta") is None or o.get("ea") in (None, ""):
            dropped.append((rs, "结局无数据")); continue
        ea_x, nea_x = h["ea"].upper(), h["nea"].upper()
        ea_y, nea_y = o["ea"].upper(), o["nea"].upper()
        b, se, eaf_y = float(o["beta"]), float(o["se"]), o.get("eaf")
        palindrome = {ea_x, nea_x} in ({"A","T"}, {"C","G"})
        if {ea_x, nea_x} == {ea_y, nea_y}:
            if ea_x == ea_y:  # 方向一致
                pass
            else:             # 等位基因互换 -> 翻转
                b = -b
                if eaf_y is not None: eaf_y = 1 - float(eaf_y)
        else:
            dropped.append((rs, "等位基因不匹配")); continue
        eaf_x = float(h["eaf"]) if h.get("eaf") not in (None, "") else None
        if palindrome:
            if eaf_x is None or eaf_y is None:
                dropped.append((rs, "回文SNP无EAF")); continue
            diff = abs(eaf_x - eaf_y)
            if min(diff, abs(diff - 1)) > 0.42 and False:
                pass
            if diff > 0.42:
                dropped.append((rs, "回文SNP EAF无法判定(差%.2f)" % diff)); continue
        rows.append(dict(rsid=rs, b_exp=float(h["beta"]), se_exp=float(h["se"]),
                         b_out=b, se_out=se, eaf_exp=eaf_x,
                         eaf_out=float(eaf_y) if eaf_y not in (None, "") else None))
    print("    保留 %d 个SNP, 剔除 %d 个: %s" % (len(rows), len(dropped), dropped or "无"))
    return rows

# ---------- MR 分析 ----------
def ratios(rows):
    r = np.array([d["b_out"]/d["b_exp"] for d in rows])
    sr = np.array([d["se_out"]/abs(d["b_exp"]) for d in rows])
    return r, sr

def ivw(r, sr):
    w = 1.0/sr**2
    b = float((w*r).sum()/w.sum())
    se_fixed = float(math.sqrt(1.0/w.sum()))
    k = len(r)
    Q = float((w*(r-b)**2).sum())
    pQ = float(stats.chi2.sf(Q, k-1)) if k > 1 else float("nan")
    se_re = se_fixed*math.sqrt(max(1.0, Q/(k-1))) if k > 1 else float("nan")
    p_fix = 2*stats.norm.sf(abs(b/se_fixed))
    p_re = 2*stats.norm.sf(abs(b/se_re)) if k > 1 else float("nan")
    return dict(b=b, se_fixed=se_fixed, se_re=se_re, p_fixed=p_fix, p_random=p_re,
                Q=Q, Q_df=k-1, Q_p=pQ, or_=math.exp(b),
                lci=math.exp(b-1.96*se_fixed), uci=math.exp(b+1.96*se_fixed))

def weighted_median(r, sr, nboot=1000):
    w = 1.0/sr**2
    idx = np.argsort(r)
    r_s, w_s = r[idx], w[idx]
    cum = np.cumsum(w_s)/w_s.sum()
    b = float(r_s[np.searchsorted(cum, 0.5)])
    boots = []
    for _ in range(nboot):
        samp = np.random.RandomState(rng.randrange(1<<30)).choice(len(r), size=len(r), replace=True,
                                                                  p=w/w.sum())
        boots.append(np.median(r[samp]))
    se = float(np.std(boots, ddof=1))
    return dict(b=b, se=se, p=2*stats.norm.sf(abs(b/se)), or_=math.exp(b),
                lci=math.exp(b-1.96*se), uci=math.exp(b+1.96*se))

def mr_egger(b_exp, b_out, se_out):
    X = np.column_stack([np.ones(len(b_exp)), b_exp])
    W = np.diag(1.0/se_out**2)
    XtW = X.T @ W
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ XtW @ b_out
    yfit = X @ beta
    k = len(b_exp)
    sigma2 = float(((1.0/se_out**2)*(b_out-yfit)**2).sum()/(k-2)) if k > 2 else 1.0
    cov = cov * sigma2
    se = np.sqrt(np.diag(cov))
    I = float(stats.norm.sf(abs(beta[1]/se[1])))
    Iint = float(stats.norm.sf(abs(beta[0]/se[0])))
    return dict(slope=float(beta[1]), se=float(se[1]), p=I,
                intercept=float(beta[0]), se_int=float(se[0]), p_int=Iint)

def weighted_mode(r, sr, nboot=1000):
    bw = 0.5*float(np.mean(sr))
    grid = np.linspace(r.min()-3*bw, r.max()+3*bw, 2000)
    def dens(rvals):
        d = np.zeros_like(grid)
        for v, s in zip(rvals, 0.5*np.mean(sr)*np.ones(len(rvals))):
            d += np.exp(-0.5*((grid-v)/s)**2)/s
        return d
    d = dens(r)
    b = float(grid[d.argmax()])
    boots = []
    for _ in range(nboot):
        samp = np.random.RandomState(rng.randrange(1<<30)).choice(len(r), size=len(r), replace=True)
        ds = dens(r[samp])
        boots.append(grid[ds.argmax()])
    se = float(np.std(boots, ddof=1))
    return dict(b=b, se=se, p=2*stats.norm.sf(abs(b/se)), or_=math.exp(b),
                lci=math.exp(b-1.96*se), uci=math.exp(b+1.96*se))

def leave_one_out(rows):
    print("[5] 留一法分析...")
    res = []
    for i in range(len(rows)):
        sub = rows[:i] + rows[i+1:]
        r, sr = ratios(sub)
        m = ivw(r, sr)
        res.append((rows[i]["rsid"], m["b"], m["se_fixed"], m["p_fixed"]))
    return res

def fmt(m):
    return "b=%.7f se=%.7f p=%.7f OR=%.4f (95%%CI %.4f-%.4f)" % (
        m.get("b", float("nan")), m.get("se_fixed", m.get("se", float("nan"))),
        m.get("p_fixed", m.get("p", float("nan"))), m["or_"], m["lci"], m["uci"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", default=None)
    ap.add_argument("--token-file", default=".opengwas_jwt")
    ap.add_argument("--exposure", default=EXPOSURE)
    ap.add_argument("--outcome", default=OUTCOME)
    args = ap.parse_args()
    import os
    tok = args.token or os.environ.get("OPENGWAS_JWT")
    if not tok and os.path.exists(args.token_file):
        tok = open(args.token_file, encoding="utf-8").read().strip()
    if not tok:
        sys.exit("缺少 OpenGWAS JWT token: 用 --token/--token-file 或环境变量 OPENGWAS_JWT 提供")

    print("== 数据集元信息 ==")
    for gid in (args.exposure, args.outcome):
        try:
            meta = api_get(tok, "/gwas/" + gid)
            meta = unwrap(meta) if isinstance(meta, dict) else meta
            if isinstance(meta, list): meta = meta[0]
            print("  %s: %s (N=%s, %s)" % (gid, meta.get("trait"), meta.get("n"), meta.get("population")))
        except Exception as e:
            print("  [警告] 元信息获取失败(%s), 继续执行..." % e)

    hits = get_tophits(tok, args.exposure)
    if not hits:
        sys.exit("  无显著SNP")
    rsids = [h.get("rsid") for h in hits]
    out_recs = get_outcome(tok, args.outcome, rsids)
    rows = harmonize(hits, out_recs)
    if len(rows) < 3:
        sys.exit("  可用SNP不足3个, 无法做MR")

    with open("harmonized_snps.csv", "w", encoding="utf-8") as f:
        f.write("rsid,b_exp,se_exp,b_out,se_out,eaf_exp,eaf_out\n")
        for d in rows:
            f.write("%s,%s,%s,%s,%s,%s,%s\n" % (d["rsid"], d["b_exp"], d["se_exp"],
                                                d["b_out"], d["se_out"], d["eaf_exp"], d["eaf_out"]))

    r, sr = ratios(rows)
    m_ivw = ivw(r, sr)
    m_wm  = weighted_median(r, sr)
    m_eg  = mr_egger(np.array([d["b_exp"] for d in rows]),
                     np.array([d["b_out"] for d in rows]),
                     np.array([d["se_out"] for d in rows]))
    m_wmo = weighted_mode(r, sr)
    loo   = leave_one_out(rows)

    print("\n===== 独立重跑结果 (k=%d SNPs) =====" % len(rows))
    print("IVW (固定效应)      :", fmt(m_ivw))
    print("IVW (乘法随机效应)  : b=%.7f se=%.7f p=%.7f OR=%.4f" % (
        m_ivw["b"], m_ivw["se_re"], m_ivw["p_random"], math.exp(m_ivw["b"])))
    print("Cochran's Q         : Q=%.3f, df=%d, P=%.4f" % (m_ivw["Q"], m_ivw["Q_df"], m_ivw["Q_p"]))
    print("加权中位数          :", fmt(m_wm))
    print("MR-Egger 斜率       : b=%.7f se=%.7f p=%.7f" % (m_eg["slope"], m_eg["se"], m_eg["p"]))
    print("MR-Egger 截距       : b=%.7f se=%.7f p=%.7f" % (m_eg["intercept"], m_eg["se_int"], m_eg["p_int"]))
    print("加权众数            :", fmt(m_wmo))
    print("\n留一法 IVW (逐一剔除):")
    for rsid, b, se, p in loo:
        print("  -%-16s b=%.4f se=%.4f p=%.4f" % (rsid, b, se, p))

    print("\n===== 与原报告对照 =====")
    print("原报告(损坏Excel) : IVW b=0.077333 se=0.159200 p=0.027137 OR=1.2804")
    print("原报告(中介TE)    : IVW b=0.077333 se=0.159200 p=0.627137 OR=1.0804")
    print("本次独立重跑      :", fmt(m_ivw))
    print("\nharmonized_snps.csv 已保存到当前目录")

if __name__ == "__main__":
    main()
