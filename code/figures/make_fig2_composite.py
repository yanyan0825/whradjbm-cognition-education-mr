"""
Figure 2 four-panel composite for the IJO Technical Report:
  A  LDSC genetic correlations (forest)          [from Fig5 A]
  B  Path decomposition of rg                    [from Fig5 B]
  C  Education-stratification ladder             [from Fig3 A]
  D  Univariable vs multivariable dumbbell       [from Fig3 B]
Counts as ONE display item. Times New Roman, 300 dpi PNG + vector PDF.
"""
import json, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")

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

C_WHR  = "#c0392b"
C_EDU  = "#1f8a70"
C_BMI  = "#2c6fbb"
C_GREY = "#9aa0a6"
C_INK  = "#1a1a1a"


def load_csv(p):
    with open(os.path.join(BASE, p), encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def panel_A(ax):
    """LDSC genetic correlations, forest with value labels."""
    d = json.load(open(os.path.join(BASE, "ldsc", "ldsc-summary.json")))
    m = {f'{x["trait1"]} | {x["trait2"]}': x for x in d["rg_matrix"]}
    # Joint multivariate LDSC (Genomic SEM) estimates: these are the values used
    # in the text and in Table S4, so the figure must match them exactly.
    JOINT = {
        "WHRadjBMI (Pulit 2019) | Cognitive performance (GCST006572)":
            (-0.1225, 0.0155, "p = 2.7e-15", "-0.122 (SE 0.016)"),
        "WHRadjBMI (Pulit 2019) | Educational attainment (GCST90566696)":
            (-0.1471, 0.0135, "p = 1.2e-27", "-0.147 (SE 0.014)"),
        "Educational attainment | Cognitive performance":
            (0.6756, 0.0133, "p < 1e-300", "+0.676 (SE 0.013)"),
    }
    pairs = [
        ("WHRadjBMI (Pulit 2019) | Cognitive performance (GCST006572)", C_WHR),
        ("WHRadjBMI (Pulit 2019) | Educational attainment (GCST90566696)", C_EDU),
        ("Educational attainment | Cognitive performance", C_BMI),
    ]
    y = 2.0
    for key, col in pairs:
        rec = m[key]
        rg, se, pstr, disp = JOINT[key]
        ax.plot([rg - 1.96 * se, rg + 1.96 * se], [y, y], color=col, lw=1.8,
                solid_capstyle="round")
        ax.plot([rg], [y], "o", ms=6.0, color=col)
        ax.text(0.05, y, f"{disp}\n{pstr}",
                va="center", fontsize=7.0, color=C_INK)
        y -= 1.0
    s = d["sensitivity_barton"]
    ax.plot([s["rg"] - 1.96 * s["se"], s["rg"] + 1.96 * s["se"]], [y, y],
            color=C_GREY, lw=1.5, ls=(0, (3, 2)))
    ax.plot([s["rg"]], [y], "o", ms=5.2, mfc="white", mec=C_GREY, mew=1.3)
    ax.text(0.05, y, f"{s['rg']:+.3f} (sensitivity,\nBarton 2021 exposure)",
            va="center", fontsize=6.8, color=C_GREY)
    ax.axvline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks([2.0, 1.0, 0.0, -1.0])
    ax.set_yticklabels(["WHRadjBMI x\ncognitive performance",
                        "WHRadjBMI x\neducational attainment",
                        "Educational attainment x\ncognitive performance",
                        "WHRadjBMI x cognition\n(exposure sensitivity)"],
                       fontsize=6.9)
    ax.set_ylim(-1.65, 2.65)
    ax.set_xlim(-0.24, 0.86)
    ax.set_xlabel("Genetic correlation ($r_g$), LDSC", fontsize=7.8)
    ax.set_title("A  Shared genetic architecture", fontsize=8.8, loc="left",
                 fontweight="bold")
    ax.tick_params(axis="y", length=0)


def panel_B(ax):
    """Path decomposition diagram (axis off)."""
    ax.axis("off")

    def node(x, y, w, h, text, ec=C_INK):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.02,rounding_size=0.05",
                                    linewidth=1.0, facecolor="#ffffff", edgecolor=ec,
                                    zorder=3, transform=ax.transAxes))
        ax.text(x, y, text, ha="center", va="center", fontsize=7.6, zorder=4,
                transform=ax.transAxes, color=C_INK)

    def arrow(p0, p1, color, lw=1.4, ls="-", rad=0.0):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9,
                                     linewidth=lw, color=color, linestyle=ls, zorder=2,
                                     connectionstyle=f"arc3,rad={rad}",
                                     transform=ax.transAxes, shrinkA=2, shrinkB=2))

    node(0.50, 0.86, 0.62, 0.15, "Educational attainment", ec=C_EDU)
    node(0.16, 0.44, 0.28, 0.17, "WHRadjBMI")
    node(0.84, 0.44, 0.28, 0.17, "Cognition")
    arrow((0.30, 0.40), (0.42, 0.80), C_EDU, rad=-0.15)
    arrow((0.70, 0.40), (0.58, 0.80), C_EDU, rad=0.15)
    arrow((0.31, 0.44), (0.69, 0.44), C_GREY, lw=1.2, ls="--")
    ax.text(0.50, 0.515, "$r_g$ = -0.122", ha="center", fontsize=6.8, color=C_GREY,
            transform=ax.transAxes, style="italic")
    ax.text(0.02, 0.68, "-0.147", fontsize=7.0, color=C_EDU,
            transform=ax.transAxes)
    ax.text(0.91, 0.68, "+0.676", fontsize=7.0, color=C_EDU,
            transform=ax.transAxes)
    ax.text(0.50, 0.285, "Product check: (-0.147)(+0.676) = -0.099  =  81%",
            ha="center", fontsize=7.0, transform=ax.transAxes, color=C_INK)
    ax.text(0.50, 0.155, "Conditional on education: $r_g$ = -0.031  =  74% attenuation",
            ha="center", fontsize=7.0, transform=ax.transAxes, color=C_EDU)
    ax.set_title("B  Decomposition of $r_g$", fontsize=8.8, loc="left",
                 fontweight="bold")


def panel_C(ax):
    """Education-stratification ladder."""
    strat = json.load(open(os.path.join(BASE, "edu_stratification.json")))
    s = [x for x in strat if x["outcome"] == "ebi-a-GCST006572"][0]
    keys = ["all", "z3", "z2", "z1"]
    lab = ["All SNPs", "z > 3", "z > 2", "|z| < 1"]
    bs = [s[k]["b"] for k in keys]
    ses = [s[k]["se"] for k in keys]
    ks = [s[k]["k"] for k in keys]
    xs = np.arange(4)
    ax.axhline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax.plot(xs, bs, color=C_WHR, lw=1.2, alpha=0.35, zorder=1)
    ann_kw = [dict(), dict(), dict(), dict(ha="right", xoff=-0.10)]
    for x, b, se, k, dy, kw in zip(xs, bs, ses, ks,
                                   [0.020, 0.036, 0.020, 0.036], ann_kw):
        ax.plot([x, x], [b - 1.96 * se, b + 1.96 * se], color=C_WHR, lw=1.6, zorder=2)
        ax.plot([x], [b], "o", ms=6.0,
                color=C_WHR if abs(b / se) > 1.96 else "white",
                mec=C_WHR, mew=1.4, zorder=3)
        ax.text(x + kw.get("xoff", 0), dy, f"{b:+.3f}\n({k} SNPs)",
                ha=kw.get("ha", "center"), fontsize=6.8, color=C_INK)
    ax.set_xticks(xs)
    ax.set_xticklabels(lab, fontsize=7.0)
    ax.set_xlim(-0.45, 3.45)
    ax.set_ylim(-0.105, 0.056)
    ax.set_ylabel("$\\beta$ (SD per SD WHRadjBMI)", fontsize=7.8)
    ax.set_title("C  Attenuation across strata of shared\n    architecture with education",
                 fontsize=8.8, loc="left", fontweight="bold")
    ax.annotate("p = 0.0013", xy=(0, bs[0]), xytext=(0.12, -0.085), fontsize=6.9,
                color=C_GREY, ha="left")
    ax.annotate("p = 0.17", xy=(3, bs[3]), xytext=(2.55, -0.085), fontsize=6.9,
                color=C_GREY, ha="left")


def panel_D(ax):
    """Univariable vs multivariable dumbbell."""
    mvmr = json.load(open(os.path.join(BASE, "mvmr_whr_edu_results.json")))
    idx = {r["outcome"]: r for r in load_csv("mr_matrix_results.csv")
           if r["exposure"] == "ebi-a-GCST90025996"}
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
        ax.plot([b0 - 1.96 * se0, b0 + 1.96 * se0], [yy + 0.10, yy + 0.10],
                color=C_GREY, lw=1.5, zorder=2)
        ax.plot([b1 - 1.96 * se1, b1 + 1.96 * se1], [yy - 0.10, yy - 0.10],
                color=C_EDU, lw=1.5, zorder=2)
        ax.plot([b0, b1], [yy + 0.10, yy - 0.10], color="#cfd4d8", lw=0.9, zorder=1)
        ax.plot([b0], [yy + 0.10], "o", ms=5.2, color=C_GREY, zorder=3)
        ax.plot([b1], [yy - 0.10], "o", ms=5.2, color=C_EDU, zorder=3)
        if o == "ebi-a-GCST006572":
            ax.annotate("", xy=(-0.038, yy - 0.34), xytext=(-0.060, yy - 0.34),
                        arrowprops=dict(arrowstyle="-|>", color=C_INK, lw=1.0,
                                        shrinkA=0, shrinkB=0))
            ax.text(-0.049, yy - 0.72, "-37%", ha="center", fontsize=7.2,
                    fontweight="bold", color=C_INK)
    ax.axvline(0, color=C_INK, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels(lab2, fontsize=7.0)
    ax.set_ylim(-2.05, len(outs) - 0.35)
    ax.set_xlim(-0.135, 0.055)
    ax.set_xlabel("$\\beta$ (95% CI)", fontsize=7.8)
    ax.set_title("D  Univariable vs. multivariable\n    (education-adjusted)",
                 fontsize=8.8, loc="left", fontweight="bold")
    ax.legend(handles=[Line2D([], [], marker="o", ls="-", color=C_GREY, ms=5.2,
                              label="Univariable IVW"),
                       Line2D([], [], marker="o", ls="-", color=C_EDU, ms=5.2,
                              label="MVMR, adjusted for education")],
              loc="lower left", bbox_to_anchor=(-0.02, 0.015), frameon=False,
              fontsize=6.9, handletextpad=0.5)


def main():
    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.08, 1.0],
                          height_ratios=[1.0, 1.0],
                          wspace=0.52, hspace=0.55,
                          left=0.16, right=0.97, top=0.90, bottom=0.14)
    panel_A(fig.add_subplot(gs[0, 0]))
    panel_B(fig.add_subplot(gs[0, 1]))
    panel_C(fig.add_subplot(gs[1, 0]))
    panel_D(fig.add_subplot(gs[1, 1]))

    fig.text(0.5, 0.015,
             "Strata in panel C are defined by the LDSC z-score of each instrument's association with "
             "educational attainment (\"All SNPs\" = unstratified).\nExecutive function and reaction time in "
             "panel D were sign-reversed to the cognitive-performance scale; conditional F statistics in the "
             "multivariable model:\nWHRadjBMI 49.4, education 20.9.",
             ha="center", fontsize=6.8, color=C_GREY)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"Fig2_fourpanel_composite.{ext}"), dpi=300)
    plt.close(fig)
    print("saved Fig2_fourpanel_composite")


if __name__ == "__main__":
    main()
