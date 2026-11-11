# Genetic correlation between central adiposity and cognitive function is shared with educational attainment

Code, results and figures supporting a Mendelian randomization (MR) re-analysis of the
adiposity–cognition association. The manuscript shows that the genetic correlation between
central adiposity (WHRadjBMI) and cognitive performance is largely carried by an axis that
also carries educational attainment: conditioning on education reduces the genetic
correlation by **74%** (−0.122 → −0.031; multivariate LD score regression / Genomic SEM),
and the MR estimate attenuates monotonically as education-linked variants are removed.

> Manuscript under review (International Journal of Obesity, Technical Report format).
> This repository contains the analysis code and key outputs; it does **not** redistribute
> GWAS summary statistics (all data are publicly available from the sources listed below).

## Repository structure

```
├── code/
│   ├── reproduce_wang2022.py      # exact replication of Wang et al. 2022 (Int J Obes 46:1487–92), see results/wang2022_replication/
│   ├── utils/                     # shared MR machinery (IVW, weighted median/mode, MR-Egger, leave-one-out, radial MR)
│   ├── main_analyses/             # primary MR matrix, MVMR, education stratification, reverse MR, outlier & APOE sensitivity
│   ├── genomic_sem/               # multivariate LDSC (GenomicSEM) + Cholesky-type genetic-correlation decomposition (R)
│   └── figures/                   # figure-generation scripts (matplotlib, Times New Roman)
├── results/                       # key outputs (JSON/CSV) for every analysis, incl. Genomic SEM S/V matrix
├── figures/                       # final figures (PNG/PDF)
└── README.md
```

## Data sources (all public)

| Trait | Source | ID / URL |
|---|---|---|
| Cognitive performance (n = 257,841) | Lee et al. 2018, UK Biobank | GWAS Catalog **GCST006572** (`ftp.ebi.ac.uk/.../GCST006572/harmonised/`) |
| WHRadjBMI (n = 458,349) | Barton et al. 2021, UKB | GWAS Catalog **GCST90025996** |
| WHRadjBMI (n = 210,082) | Shungin et al. 2015, GIANT | OpenGWAS **ieu-a-79** |
| BMI (n = 339,224) | Locke et al. 2015, GIANT | OpenGWAS **ieu-a-2** |
| Educational attainment (n = 753,152) | Xia et al. 2025 | GWAS Catalog **GCST90566696** |
| Educational attainment, UKB-free (n = 106,736) | Rietveld et al. 2014 | GWAS Catalog **GCST002598** |

LD score reference: 1000 Genomes Phase 3 EUR (Broad/MRC-IEU LDSC files, `eur_w_ld_chr` + `w_hm3.snplist`).

## Requirements

- Python ≥ 3.10: `numpy`, `scipy`, `requests` (+ `pandas`, `matplotlib` for figures)
- R ≥ 4.2 with [`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM) (for `code/genomic_sem/`)
- [OpenGWAS](https://api.opengwas.io) account: put your JWT in a file named `.opengwas_jwt`
  next to the scripts (never commit it).

## Reproducing the analyses

```bash
# 1. Replication of Wang et al. 2022 (uses local GCST006572 + OpenGWAS instruments)
python code/reproduce_wang2022.py
# -> matches published IVW beta = -0.07 [-0.12, -0.02] (our: -0.089 [-0.143, -0.035])

# 2. Main MR matrix (all exposure × outcome cells, primary + sensitivity estimators)
python code/main_analyses/mr_matrix.py

# 3. Multivariable MR (WHRadjBMI + education), education stratification, reverse MR
python code/main_analyses/mvmr_whr_edu.py
python code/main_analyses/edu_strat.py
python code/main_analyses/mr_reverse.py

# 4. UKB-free education sensitivity (Rietveld 2014)
python code/main_analyses/edu_sensitivity_rietveld.py

# 5. Multivariate LDSC + genetic-correlation decomposition (R)
Rscript code/genomic_sem/run_ldsc3.R     # S and V matrices
Rscript code/genomic_sem/decompose.R     # education-conditional rg: -0.031 (74% attenuation)

# 6. Figures
python code/figures/make_figures.py
python code/figures/make_fig2_composite.py
```

Pre-computed outputs for steps 1–5 are in `results/`.

## Headline results

| Quantity | Estimate | 95% CI |
|---|---|---|
| rg(WHRadjBMI, cognition) | −0.122 | −0.153 to −0.092 |
| rg(WHRadjBMI, cognition \| education) | **−0.031** | −0.071 to +0.009 |
| Attenuation | **74%** | — |
| Univariable MR (IVW, composite) | −0.060 | −0.096 to −0.025 |
| MVMR (education-adjusted) | −0.038 | −0.061 to −0.015 |

## Citation

If you use this code, please cite the manuscript (citation to be added upon acceptance)
and the original GWAS/method papers listed in the manuscript reference list.

## License

MIT — see [LICENSE](LICENSE).
