# An Efficient Rotation-Sign-Permutation Algorithm to Solve Rotational Indeterminacy in Bayesian Exploratory Factor Analysis

**Authors:** Ricardo Rey-Sáez & Javier Revuelta

**Affiliations:** Universidad Autónoma de Madrid

---

## Overview

This repository contains all data, analysis scripts, simulation code, and manuscript sources needed to reproduce the findings reported in the paper. The work introduces an **Efficient RSP (E-RSP) algorithm** that overcomes the limitations of existing alignment methods for Bayesian Exploratory Factor Analysis (EFA). The E-RSP algorithm makes exact signed-permutation alignment scalable to a large number of factors by reducing the alignment step to a single **Linear Assignment Problem (LAP)**, solved with the Hungarian algorithm—eliminating the exponential enumeration of $2^M$ sign configurations required by the original RSP method.

The algorithm is implemented in C++ and available in the R package [`BayesEFA`](https://github.com/RicardoReySaez/BayesEFA).

> **Open Science:** This repository accompanies the OSF project at https://osf.io/4nyd7/.

---

## Repository Structure

```
ERSP/
│
├── ERSP.Rproj                  # RStudio project file (sets working directory)
├── README.md                   # This file
│
├── Data/
│   └── df_chem.rds                # NHANES 2015-2016 chemical exposure dataset
│                                    (used in the empirical illustration)
│
├── Figures/
│   ├── RSP_agreement.png          # Scatterplot: perfect agreement in absolute fit
│   │                                errors between exact RSP and E-RSP (Simulation 1)
│   └── factor_loading_matrix.png  # Posterior means and 99% CIs of factor loadings
│                                    aligned via E-RSP (Lopes & West illustration)
│
├── R scripts/
│   ├── R functions/
│   │   └── helper_fcts.R              # Auxiliary functions from Poworoznek et al. (2025)
│   │
│   ├── Simulation study/
│   │   ├── RSP_comparison_sim.R       # Simulation Study 1: run E-RSP vs. exact RSP
│   │   ├── RSP_comparison_analysis.R  # Simulation Study 1: analyze results
│   │   ├── MA_comparison_sim.R        # Simulation Study 2: run E-RSP vs. MatchAlign
│   │   └── MA_comparison_analysis.R   # Simulation Study 2: analyze results
│   │
│   └── Empirical analysis/
│       ├── NHANES_illustration.R      # Empirical: high-dimensional NHANES analysis
│       │                                (E-RSP vs. MatchAlign, M = 7 factors)
│       └── Lopes_illustration.R       # Empirical: Lopes & West exchange-rate data
│                                        (exact RSP vs. MatchAlign vs. E-RSP)
│
├── Results/
│   ├── Simulation study/
│   │   ├── RSP comparison/
│   │   │   ├── rsp_comparison_simresults.rds      # Raw simulation output (Study 1)
│   │   │   ├── rsp_comparison_simresults.rdata    # Same as above (.RData format)
│   │   │   └── simres1_execution_time.rds         # Summarized execution times
│   │   │
│   │   └── MatchAlign comparison/
│   │       ├── MatchAlign_comp_simres.rds         # Raw simulation output (Study 2)
│   │       ├── MatchAlign_comp_simres.rdata       # Same as above (.RData format)
│   │       └── simres2_marginal_means.rds         # Summarized marginal means
│   │
│   └── Empirical analysis/
│       ├── NHANES_fits.rds            # Fitted MCMC objects for the NHANES analysis
│       │                                (~2.5 GB; NOT included in this repo due to
│       │                                 size — available on OSF: https://osf.io/4nyd7/)
│       ├── NHANES_results_table.rds   # Summary table of NHANES results
│       └── LW_results_table.rds       # Summary table of Lopes & West results
│
└── Manuscript/
    ├── Manuscript.qmd           # Quarto source: narrative + embedded R code for
    │                              dynamic table generation
    ├── Manuscript.tex           # LaTeX output rendered from the .qmd
    ├── Manuscript.pdf           # Compiled PDF of the manuscript
    ├── bibliography.bib         # BibTeX references
    └── _extensions/
        └── wjschne/apaquarto/   # Quarto extension for APA 7th ed. formatting
```

---

## Simulation Studies

| Study | Comparison | Purpose | Script |
|-------|-----------|---------|--------|
| **1** | E-RSP vs. exact RSP | Verify that E-RSP yields identical solutions in low-dimensional settings ($M \leq 10$) | `RSP_comparison_sim.R` |
| **2** | E-RSP vs. MatchAlign | Benchmark alignment accuracy and runtime as dimensionality increases ($M$ up to 50) | `MA_comparison_sim.R` |

## Empirical Illustrations

| Example | Dataset | Factors | Script |
|---------|---------|---------|--------|
| **NHANES** | Chemical exposure biomarkers (2015-2016) | $M = 25$ and $M = 50$ | `NHANES_illustration.R` |
| **Lopes & West** | Exchange-rate data | $M = 5$ | `Lopes_illustration.R` |

---

## Software & Dependencies

The analyses were conducted in **R (≥ 4.5.0)**. Key packages:

| Package | Role |
|---------|------|
| [`BayesEFA`](https://github.com/RicardoReySaez/BayesEFA) | E-RSP algorithm (C++ implementation) |
| `factor.switching` | Exact RSP algorithm of Papastamoulis & Ntzoufras (2022) |
| `infinitefactor` | MatchAlign algorithm, `linearDL` and `linearMGSP` samplers |
| `SimDesign` | Reproducible simulation framework |
| `mvnfast` | Fast multivariate normal data generation |
| `tidyverse` | Data manipulation and visualization |

---

## Reproducing the Results

1. **Clone** this repository (or download from OSF) preserving the directory structure.
2. **Open** `ERSP.Rproj` in RStudio to set the working directory.
3. **Run simulations**: execute scripts in `R scripts/Simulation study/`.
4. **Run empirical analyses**: execute scripts in `R scripts/Empirical analysis/`.
5. **Inspect outputs**: tables and intermediate results are saved to `Results/`; figures to `Figures/`.
6. **Render the manuscript**: run `quarto render Manuscript/Manuscript.qmd` to produce the PDF with embedded results.

---

## Citation

If you use this work, please cite:

> Rey-Sáez, R. & Revuelta, J. (2026). An Efficient Rotation-Sign-Permutation Algorithm to Solve Rotational Indeterminacy in Bayesian Exploratory Factor Analysis. *Manuscript submitted for publication*.

---

## Contact

Ricardo Rey-Sáez — [ricardoreysaez95@gmail.com](mailto:ricardoreysaez95@gmail.com)  
Department of Basic Psychology, Universidad Autónoma de Madrid
