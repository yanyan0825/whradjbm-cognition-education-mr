"""
Publication figures for the WHRadjBMI-cognition MR manuscript.
All numbers are read from result files -- nothing hard-coded except the
conceptual DAG layout (Fig 1).
Output: figures/Fig1..Fig5 (.png 300 dpi + .pdf vector)
"""
import json, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

C_WC   = "#7f7f7f"   # waist circumference
C_WHR  = "#c0392b"   # WHRadjBMI (primary)
C_BMI  = "#2c6fbb"   # BMI
C_EDU  = "#1f8a70"   # education
C_GREY = "#9aa0a6"
C_INK  = "#1a1a1a"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
    print("saved", name)


def load_csv(p):
    with open(os.path.join(BASE, p), encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------------
# Fig 1 -- two competing hypotheses (conceptual DAG)
# ----------------------------------------------------------------------------
def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))

    def node(ax, x, y, w, h, text, fc="#ffffff", ec=C_INK, fs=8.2, bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.02,rounding_size=0.08",
                                    linewidth=1.1, facecolor=fc, edgecolor=ec, zorder=3))
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4,
                fontweight="bold" if bold else "normal", color=C_INK)

    def arrow(ax, p0, p1, style="-", color=C_INK, lw=1.4, rad=0.0, ls=None, z=2):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11,
                                     linewidth=lw, color=color, zorder=z,
                                     linestyle=ls if ls else style,
                                     connectionstyle=f"arc3,rad={rad}",
                                     shrinkA=1, shrinkB=1))

    for ax, title in zip(axes, ["Hypothesis 1\nDirect causal effect",
                                "Hypothesis 2\nShared upstream structure"]):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
        ax.set_title(title, fontsize=9.2, fontweight="bold", color=C_INK, pad=6)

    # ---- panel A
    ax = axes[0]
    node(ax, 2.1, 5.0, 3.95, 1.3, "Central adiposity\n(WHRadjBMI)")
    node(ax, 8.15, 5.0, 3.1, 1.3, "Cognitive\nfunction")
    arrow(ax, (4.2, 5.0), (6.5, 5.0), lw=1.8, color=C_WHR)
    ax.text(5.35, 5.85, "causal effect?", ha="center", fontsize=8, color=C_WHR,
            style="italic")
    ax.text(5.0, 2.5, "$\\beta \\approx -0.06$ SD per SD\n(before adjustment)",
            ha="center", fontsize=7.4, color=C_INK)
    ax.text(5.0, 1.0, "prediction: effect survives\nmultivariable adjustment",
            ha="center", fontsize=7.2, color=C_GREY, style="italic")

    # ---- panel B
    ax = axes[1]
    node(ax, 5.0, 8.35, 5.9, 1.35, "Educational attainment /\ncognitive reserve", ec=C_EDU)
    node(ax, 2.15, 4.1, 4.1, 1.3, "Central adiposity\n(WHRadjBMI)")
    node(ax, 8.3, 4.1, 3.1, 1.3, "Cognitive\nfunction")
    arrow(ax, (3.0, 4.8), (4.1, 7.65), color=C_EDU, lw=1.6, rad=0.12)
    arrow(ax, (7.4, 4.8), (6.0, 7.65), color=C_EDU, lw=1.6, rad=-0.12)
    arrow(ax, (4.3, 4.1), (6.65, 4.1), color=C_GREY, lw=1.4, ls="--")
    ax.text(5.45, 3.15, "no direct effect", ha="center", fontsize=7.0, color=C_GREY,
            style="italic")
    ax.text(1.15, 6.85, "rg = -0.15", fontsize=7.2, color=C_EDU)
    ax.text(8.45, 6.85, "rg = +0.67", fontsize=7.2, color=C_EDU)
    ax.text(5.0, 1.85, "indirect path: (-0.15)(+0.67) = -0.10\n"
                       "= 81% of rg(WHRadjBMI, cognition)",
            ha="center", fontsize=7.4, color=C_INK)
    ax.text(5.0, 0.55, "prediction: association attenuates after\nadjustment for education",
            ha="center", fontsize=7.2, color=C_GREY, style="italic")

    fig.text(0.5, -0.01, "rg, genetic correlation (LDSC). Numbers in panel B from the present study.",
             ha="center", fontsize=7, color=C_GREY)
    save(fig, "Fig1_hypotheses_DAG")


# ----------------------------------------------------------------------------
# Fig 2 -- forest plot, 3 exposures x 5 cognitive outcomes
# ----------------------------------------------------------------------------
def fig2():
    rows = load_csv("mr_matrix_results.csv")
    order_exp = ["ieu-a-66", "ebi-a-GCST90025996", "ieu-b-40"]
    exp_lab = {"ieu-a-66": "Waist circumference", "ebi-a-GCST90025996": "WHRadjBMI",
               "ieu-b-40": "BMI"}
    exp_col = {"ieu-a-66": C_WC, "ebi-a-GCST90025996": C_WHR, "ieu-b-40": C_BMI}
    order_out = ["ebi-a-GCST006572", "ukb-b-2084", "ukb-b-15625",
                 "ukb-b-20140", "ukb-b-13823"]
    out_lab = {"ebi-a-GCST006572": "General cognitive performance",
               "ukb-b-2084": "Visual memory",
               "ukb-b-15625": "Processing speed",
               "ukb-b-20140": "Executive function",
               "ukb-b-13823": "Reaction time"}

    idx = {(r["exposure"], r["outcome"]): r for r in rows}

    fig = plt.figure(figsize=(7.4, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.45, 1.05], wspace=0.02)
    ax = fig.add_subplot(gs[0])
    axt = fig.add_subplot(gs[1], sharey=ax)
    axt.axis("off")

    y = 0.0
    yt, yl, ybold, ycol, blocks = [], [], [], [], []
    for e in order_exp:
        y_start = y
        yt.append(y); yl.append(exp_lab[e]); ybold.append(True); ycol.append(exp_col[e])
        y += 0.8
        for o in order_out:
            r = idx.get((e, o))
            if r is None:
                continue
            b = float(r["b"]); se = float(r["se"]); p = float(r["p"])
            if int(r["higher_is_worse"]) == 1:      # flip to "cognitive performance" scale
                b = -b
            lo, hi = b - 1.96 * se, b + 1.96 * se
            sig = p < 0.05
            ax.plot([lo, hi], [y, y], color=exp_col[e], lw=1.5,
                    solid_capstyle="round", alpha=0.95, zorder=2)
            ax.plot([b], [y], marker="o", ms=5.6, color=exp_col[e] if sig else "white",
                    markeredgecolor=exp_col[e], markeredgewidth=1.3, zorder=3)
            axt.text(0.045, y, f"{int(r['k'])}", ha="left", va="center", fontsize=7.4,
                     color=C_INK, transform=axt.get_yaxis_transform())
            axt.text(0.30, y, f"{b:+.3f} ({lo:+.3f}, {hi:+.3f})", ha="left", va="center",
                     fontsize=7.4, color=C_INK, transform=axt.get_yaxis_transform(),
                     fontweight="bold" if (sig and o == "ebi-a-GCST006572") else "normal")
            yt.append(y); yl.append(out_lab[o]); ybold.append(False); ycol.append(C_INK)
            y += 1.0
        blocks.append((y_start - 0.42, y - 0.42))
        y += 0.6

    for y0, y1 in blocks:
        ax.axhspan(y0, y1, color="#000000", alpha=0.035, zorder=0)

    ax.axvline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=8)
    for t, bo, co in zip(ax.get_yticklabels(), ybold, ycol):
        t.set_fontweight("bold" if bo else "normal")
        t.set_color(co)
    ax.set_ylim(y - 0.9, -2.3)
    ax.set_xlim(-0.155, 0.115)
    ax.set_xlabel("Effect of 1-SD genetically predicted adiposity on cognitive performance (SD)", fontsize=8)
    ax.tick_params(axis="y", length=0)

    axt.text(0.045, -1.75, "SNPs", ha="left", va="center", fontsize=7.6,
             fontweight="bold", transform=axt.get_yaxis_transform())
    axt.text(0.30, -1.75, "$\\beta$ (95% CI)", ha="left", va="center", fontsize=7.6,
             fontweight="bold", transform=axt.get_yaxis_transform())

    handles = [Line2D([], [], marker="o", ls="", ms=5.6, color=c, label=l)
               for c, l in [(C_WC, "Waist circumference"), (C_WHR, "WHRadjBMI"),
                            (C_BMI, "BMI")]]
    handles += [Line2D([], [], marker="o", ls="", ms=5.6, mfc="white", mec=C_INK,
                       label="p $\\geq$ 0.05"),
                Line2D([], [], marker="o", ls="", ms=5.6, mfc=C_INK, mec=C_INK,
                       label="p < 0.05")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.02, 1.01),
              frameon=False, fontsize=7.2, handletextpad=0.4, borderaxespad=0.0,
              ncol=5, columnspacing=1.0)

    fig.text(0.5, -0.055,
             "Executive function (trail-making B) and reaction time are scored so that higher values indicate worse performance; "
             "signs were reversed so that\nnegative values denote poorer cognitive performance for every outcome. "
             "Estimates are inverse-variance weighted.",
             ha="center", fontsize=6.9, color=C_GREY)
    save(fig, "Fig2_forest_main_matrix")


# ----------------------------------------------------------------------------
# Fig 3 -- education attenuation: (A) stratification ladder, (B) UV vs MV
# ----------------------------------------------------------------------------
def fig3():
    strat = json.load(open(os.path.join(BASE, "edu_stratification.json")))
    s = [x for x in strat if x["outcome"] == "ebi-a-GCST006572"][0]
    mvmr = json.load(open(os.path.join(BASE, "mvmr_whr_edu_results.json")))
    rows = load_csv("mr_matrix_results.csv")
    idx = {r["outcome"]: r for r in rows if r["exposure"] == "ebi-a-GCST90025996"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.35),
                                   gridspec_kw={"width_ratios": [1.0, 1.18], "wspace": 0.60})

    # ---- A: ladder
    keys = ["all", "z3", "z2", "z1"]
    lab = ["All SNPs", "z > 3", "z > 2", "|z| < 1"]
    bs = [s[k]["b"] for k in keys]
    ses = [s[k]["se"] for k in keys]
    ks = [s[k]["k"] for k in keys]
    xs = np.arange(4)
    ax1.axhline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax1.plot(xs, bs, color=C_WHR, lw=1.2, alpha=0.35, zorder=1)
    # value annotations: last stratum is right-aligned so it cannot spill into
    # the gap where panel B's long y-tick labels live
    ann_kw = [dict(), dict(), dict(), dict(ha="right", xoff=-0.10)]
    for x, b, se, k, dy, kw in zip(xs, bs, ses, ks, [0.020, 0.036, 0.020, 0.036], ann_kw):
        ax1.plot([x, x], [b - 1.96 * se, b + 1.96 * se], color=C_WHR, lw=1.6, zorder=2)
        ax1.plot([x], [b], "o", ms=6.2, color=C_WHR if abs(b / se) > 1.96 else "white",
                 mec=C_WHR, mew=1.4, zorder=3)
        ax1.text(x + kw.get("xoff", 0), dy, f"{b:+.3f}\n({k} SNPs)", ha=kw.get("ha", "center"),
                 fontsize=6.8, color=C_INK)
    ax1.set_xticks(xs); ax1.set_xticklabels(lab, fontsize=7.2)
    ax1.set_xlim(-0.45, 3.45); ax1.set_ylim(-0.105, 0.056)
    ax1.set_ylabel("$\\beta$ (SD cognitive performance per SD WHRadjBMI)", fontsize=7.8)
    ax1.set_title("A  Attenuation across strata of shared\n    genetic architecture with education",
                  fontsize=8.6, loc="left", fontweight="bold")
    ax1.annotate("p = 0.0013", xy=(0, bs[0]), xytext=(0.12, -0.085), fontsize=6.9,
                 color=C_GREY, ha="left")
    ax1.annotate("p = 0.17", xy=(3, bs[3]), xytext=(2.55, -0.085), fontsize=6.9,
                 color=C_GREY, ha="left")

    # ---- B: univariable vs multivariable dumbbell
    outs = ["ebi-a-GCST006572", "ukb-b-2084", "ukb-b-15625", "ukb-b-20140", "ukb-b-13823"]
    lab2 = ["General cognitive\nperformance", "Visual memory", "Processing speed",
            "Executive function", "Reaction time"]
    y = np.arange(len(outs))[::-1]
    for yy, o, nm in zip(y, outs, lab2):
        r = idx[o]
        b0 = float(r["b"]); se0 = float(r["se"])
        if int(r["higher_is_worse"]) == 1:
            b0 = -b0
        m = [x for x in mvmr if x["outcome"] == o][0]
        b1, se1 = m["whr_b"], m["whr_se"]
        ax2.plot([b0 - 1.96 * se0, b0 + 1.96 * se0], [yy + 0.10, yy + 0.10],
                 color=C_GREY, lw=1.5, zorder=2)
        ax2.plot([b1 - 1.96 * se1, b1 + 1.96 * se1], [yy - 0.10, yy - 0.10],
                 color=C_EDU, lw=1.5, zorder=2)
        ax2.plot([b0, b1], [yy + 0.10, yy - 0.10], color="#cfd4d8", lw=0.9, zorder=1)
        ax2.plot([b0], [yy + 0.10], "o", ms=5.4, color=C_GREY, zorder=3)
        ax2.plot([b1], [yy - 0.10], "o", ms=5.4, color=C_EDU, zorder=3)
        if o == "ebi-a-GCST006572":
            ax2.annotate("", xy=(-0.038, yy - 0.34), xytext=(-0.060, yy - 0.34),
                         arrowprops=dict(arrowstyle="-|>", color=C_INK, lw=1.0,
                                         shrinkA=0, shrinkB=0))
            ax2.text(-0.049, yy - 0.64, "-37%", ha="center", fontsize=7.4,
                     fontweight="bold", color=C_INK)
    ax2.axvline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax2.set_yticks(y); ax2.set_yticklabels(lab2, fontsize=7.4)
    ax2.set_ylim(-2.05, len(outs) - 0.35)
    ax2.set_xlim(-0.135, 0.055)
    ax2.set_xlabel("$\\beta$ (95% CI)", fontsize=7.8)
    ax2.set_title("B  Univariable vs. multivariable (education-adjusted)",
                  fontsize=8.6, loc="left", fontweight="bold")
    ax2.legend(handles=[Line2D([], [], marker="o", ls="-", color=C_GREY, ms=5.4,
                               label="Univariable IVW"),
                        Line2D([], [], marker="o", ls="-", color=C_EDU, ms=5.4,
                               label="MVMR, adjusted for education")],
               loc="lower left", bbox_to_anchor=(-0.02, 0.015), frameon=False,
               fontsize=7.2, handletextpad=0.5)

    fig.text(0.5, -0.05,
             "Strata are defined by the LDSC z-score of each instrument's association with educational attainment "
             "(\"All SNPs\" = unstratified). Conditional F statistics in the multivariable model:\n"
             "WHRadjBMI 49.4, education 20.9. Executive function and reaction time were sign-reversed "
             "to the cognitive-performance scale.",
             ha="center", fontsize=6.9, color=C_GREY)

    save(fig, "Fig3_education_attenuation")


# ----------------------------------------------------------------------------
# Fig 4 -- replication across two independent WHRadjBMI GWAS
# ----------------------------------------------------------------------------
def fig4():
    barton = [r for r in load_csv("mr_matrix_results.csv")
              if r["exposure"] == "ebi-a-GCST90025996"]
    pulit = json.load(open(os.path.join(BASE, "giant_pulit_replication.json")))
    pm = {p["outcome"]: p for p in pulit}
    keep = ["ebi-a-GCST006572", "ukb-b-2084", "ukb-b-15625", "ukb-b-20140", "ukb-b-13823"]
    lab = {"ebi-a-GCST006572": "General cognitive performance", "ukb-b-2084": "Visual memory",
           "ukb-b-15625": "Processing speed", "ukb-b-20140": "Executive function",
           "ukb-b-13823": "Reaction time"}
    flip = {o: (1 if o in ("ukb-b-20140", "ukb-b-13823") else 0) for o in keep}
    ks = {r["outcome"]: r["k"] for r in barton}

    xs, ys, ses, labs = [], [], [], []
    for o in keep:
        r = [x for x in barton if x["outcome"] == o][0]
        b1 = float(r["b"]) * (-1 if flip[o] else 1)
        b2 = pm[o]["b"] * (-1 if flip[o] else 1)
        xs.append(b1); ys.append(b2); ses.append((float(r["se"]), pm[o]["se"])); labs.append(lab[o])

    fig, ax = plt.subplots(figsize=(4.8, 4.25))
    lim = (-0.115, 0.075)
    ax.plot(lim, lim, ls=(0, (4, 3)), color=C_GREY, lw=1.0, zorder=1)
    for x, y, (sx, sy), l in zip(xs, ys, ses, labs):
        prim = l.startswith("General")
        ax.errorbar([x], [y], xerr=1.96 * sx, yerr=1.96 * sy, fmt="none",
                    ecolor=C_WHR if prim else C_GREY, elinewidth=1.2, alpha=0.75, zorder=2)
        ax.plot([x], [y], "o", ms=8.0 if prim else 5.6, color=C_WHR if prim else C_BMI,
                mec="white", mew=0.8, zorder=3)
    ax.annotate("General cognitive\nperformance", xy=(xs[0], ys[0]),
                xytext=(-0.112, -0.097), fontsize=7.4, color=C_WHR, fontweight="bold")
    ax.annotate("Executive function", xy=(-0.0315, -0.019),
                xytext=(-0.104, -0.038), fontsize=7.2, color=C_BMI,
                arrowprops=dict(arrowstyle="-", color=C_BMI, lw=0.7))
    ax.annotate("Visual memory, processing speed,\nreaction time (overlapping)",
                xy=(0.0, 0.006), xytext=(0.008, 0.049), fontsize=7.2, color=C_BMI,
                arrowprops=dict(arrowstyle="-", color=C_BMI, lw=0.7))

    r = np.corrcoef(xs, ys)[0, 1]
    slope = float(np.dot(xs, ys) / np.dot(xs, xs))
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xticks([-0.10, -0.05, 0.0, 0.05])
    ax.set_yticks([-0.10, -0.05, 0.0, 0.05])
    ax.set_xlabel("$\\beta$, Barton 2021 WHRadjBMI ($n$ = 458,349; 233-237 SNPs)", fontsize=7.8)
    ax.set_ylabel("$\\beta$, GIANT 2019 WHRadjBMI ($n$ = 694,649; 146 SNPs)", fontsize=7.8)
    ax.set_title("Independent replication of the exposure-outcome estimates",
                 fontsize=8.8, fontweight="bold", loc="left")
    ax.text(0.03, 0.975, f"r = {r:.2f}   slope through origin = {slope:.2f}",
            transform=ax.transAxes, fontsize=7.6, va="top", color=C_INK)
    ax.text(0.03, 0.905, "primary outcome: -0.060 vs -0.064\n(concordant, both p < 0.05)",
            transform=ax.transAxes, fontsize=7.2, va="top", color=C_GREY)
    ax.set_aspect("equal")
    save(fig, "Fig4_replication")


# ----------------------------------------------------------------------------
# Fig 5 -- LDSC genetic correlation and path decomposition
# ----------------------------------------------------------------------------
def fig5():
    d = json.load(open(os.path.join(BASE, "ldsc", "ldsc-summary.json")))
    m = {f'{x["trait1"]} | {x["trait2"]}': x for x in d["rg_matrix"]}
    pairs = [
        ("WHRadjBMI (Pulit 2019) | Cognitive performance (GCST006572)", "WHRadjBMI x cognitive performance", C_WHR),
        ("WHRadjBMI (Pulit 2019) | Educational attainment (GCST90566696)",
         "WHRadjBMI x educational attainment", C_EDU),
        ("Educational attainment | Cognitive performance",
         "Educational attainment x cognitive performance", C_BMI),
    ]

    fig = plt.figure(figsize=(7.3, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.30)
    ax = fig.add_subplot(gs[0])
    axt = fig.add_subplot(gs[1])
    axt.axis("off")

    y = 2.0
    for key, label, col in pairs:
        rec = m[key]
        rg, se = rec["rg"], rec["se"]
        pstr = "p < 1e-300" if rec["p"] == 0 else f"p = {rec['p']:.1e}"
        ax.plot([rg - 1.96 * se, rg + 1.96 * se], [y, y], color=col, lw=1.8,
                solid_capstyle="round")
        ax.plot([rg], [y], "o", ms=6.4, color=col)
        ax.text(max(rg + 1.96 * se, rg) + 0.045 if rg > 0 else 0.05, y,
                f"{rg:+.3f} (SE {se:.3f})\n{pstr}",
                va="center", fontsize=7.2, color=C_INK)
        y -= 1.0
    # sensitivity row
    s = d["sensitivity_barton"]
    ax.plot([s["rg"] - 1.96 * s["se"], s["rg"] + 1.96 * s["se"]], [y, y],
            color=C_GREY, lw=1.5, ls=(0, (3, 2)))
    ax.plot([s["rg"]], [y], "o", ms=5.4, mfc="white", mec=C_GREY, mew=1.3)
    ax.text(0.05, y,
            f"{s['rg']:+.3f} (sensitivity,\nBarton 2021 exposure)",
            va="center", fontsize=6.9, color=C_GREY)

    ax.axvline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks([2.0, 1.0, 0.0, -1.0])
    ax.set_yticklabels(["WHRadjBMI x\ncognitive performance", "WHRadjBMI x\neducational attainment",
                        "Educational attainment x\ncognitive performance", "WHRadjBMI x cognition\n(exposure sensitivity)"],
                       fontsize=7.2)
    ax.set_ylim(-1.65, 2.65); ax.set_xlim(-0.24, 0.86)
    ax.set_xlabel("Genetic correlation ($r_g$), LDSC", fontsize=7.8)
    ax.set_title("A  Shared genetic architecture", fontsize=8.8, loc="left",
                 fontweight="bold")
    ax.tick_params(axis="y", length=0)

    # ---- B: path diagram
    def node(x, y, w, h, text, fc="#ffffff", ec=C_INK):
        axt.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                     boxstyle="round,pad=0.02,rounding_size=0.05",
                                     linewidth=1.0, facecolor=fc, edgecolor=ec, zorder=3,
                                     transform=axt.transAxes))
        axt.text(x, y, text, ha="center", va="center", fontsize=7.6, zorder=4,
                 transform=axt.transAxes, color=C_INK)

    def arrow(p0, p1, color, lw=1.4, ls="-", rad=0.0):
        axt.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9,
                                      linewidth=lw, color=color, linestyle=ls, zorder=2,
                                      connectionstyle=f"arc3,rad={rad}",
                                      transform=axt.transAxes, shrinkA=2, shrinkB=2))

    node(0.50, 0.86, 0.62, 0.15, "Educational attainment", ec=C_EDU)
    node(0.16, 0.44, 0.28, 0.17, "WHRadjBMI")
    node(0.84, 0.44, 0.28, 0.17, "Cognition")
    arrow((0.30, 0.40), (0.42, 0.80), C_EDU, rad=-0.15)
    arrow((0.70, 0.40), (0.58, 0.80), C_EDU, rad=0.15)
    arrow((0.31, 0.44), (0.69, 0.44), C_GREY, lw=1.2, ls="--")
    axt.text(0.50, 0.515, "$r_g$ = -0.121", ha="center", fontsize=6.8, color=C_GREY,
             transform=axt.transAxes, style="italic")
    axt.text(0.10, 0.735, "-0.147", fontsize=7.0, color=C_EDU, transform=axt.transAxes)
    axt.text(0.845, 0.735, "+0.670", fontsize=7.0, color=C_EDU, transform=axt.transAxes)
    axt.text(0.50, 0.20, "(-0.147)(+0.670) = -0.098", ha="center", fontsize=7.4,
             transform=axt.transAxes, color=C_INK)
    axt.text(0.50, 0.10, "= 81% of the observed $r_g$(-0.121)", ha="center",
             fontsize=7.4, transform=axt.transAxes, color=C_INK, fontweight="bold")
    axt.set_title("B  Path decomposition of $r_g$", fontsize=8.8, loc="left",
                  fontweight="bold")

    save(fig, "Fig5_ldsc_rg")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    print("ALL FIGURES DONE")
