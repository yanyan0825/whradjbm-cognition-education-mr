# -*- coding: utf-8 -*-
"""WHRadjBMI 各认知结局的教育多效性分层（主分析用表）"""
import json, math
import numpy as np
from scipy import stats

cache = json.load(open("matrix_cache.json"))
mc = json.load(open("whradjbmi_assoc_cache.json"))          # 含 ieu-a-1239 教育效应
exp = cache["tophits_ebi-a-GCST90025996"]
edu = mc["ieu-a-1239"]
OUT = [("ebi-a-GCST006572", "一般认知表现"), ("ukb-b-2084", "视觉记忆"),
       ("ukb-b-15625", "加工速度"), ("ukb-b-20140", "执行功能"), ("ukb-b-13823", "反应时")]


def ivw(rows):
    bx = np.array([r["b_exp"] for r in rows]); by = np.array([r["b_out"] for r in rows])
    sy = np.array([r["se_out"] for r in rows])
    r = by / bx; w = 1 / (sy / np.abs(bx)) ** 2
    b = float((w * r).sum() / w.sum()); se_f = float(math.sqrt(1 / w.sum()))
    k = len(rows); Q = float((w * (r - b) ** 2).sum()); se = se_f * math.sqrt(max(1, Q / (k - 1)))
    return dict(b=b, se=se, p=float(2 * stats.norm.sf(abs(b / se))), k=k)


def build(oid):
    om = cache["assoc"].get(oid, {})
    rows = []
    for e in exp:
        o = om.get(e["rsid"])
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
        z = None
        if e["rsid"] in edu:
            z = edu[e["rsid"]]["beta"] / edu[e["rsid"]]["se"]
        rows.append(dict(rsid=e["rsid"], b_exp=float(e["beta"]), se_exp=float(e["se"]),
                         b_out=float(by), se_out=float(o["se"]), z_edu=z))
    return rows


print("%-14s %-26s %-26s %-26s %-26s" % ("结局", "全部", "剔除|z|>3", "剔除|z|>2", "仅|z|<1"))
out = []
for oid, name in OUT:
    rows = build(oid)
    line = [name]
    rec = dict(outcome=oid, name=name)
    for lab, flt in [("all", None), ("z3", 3), ("z2", 2), ("z1", 1)]:
        if flt is None:
            sub = rows
        elif lab == "z1":
            sub = [r for r in rows if r["z_edu"] is not None and abs(r["z_edu"]) < 1]
        else:
            sub = [r for r in rows if r["z_edu"] is not None and abs(r["z_edu"]) <= flt]
        if len(sub) < 10:
            line.append("k<10"); rec[lab] = None; continue
        m = ivw(sub)
        line.append("%+.3f (p=%.3f, k=%d)" % (m["b"], m["p"], m["k"]))
        rec[lab] = dict(b=round(m["b"], 4), se=round(m["se"], 4), p=round(m["p"], 4), k=m["k"])
    print("%-14s " % name + " ".join("%-26s" % x for x in line[1:]))
    out.append(rec)
json.dump(out, open("edu_stratification.json", "w"), ensure_ascii=False, indent=1)
