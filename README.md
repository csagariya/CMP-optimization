# Controlled mass pollination optimization

This repository supports the manuscript on controlled mass pollination (CMP)
optimization using general combining ability (GCA), specific combining ability
(SCA), and an explicit status-number constraint.

## Quantitative-genetic conventions used by the code

The repository uses **standard diallel GCA effects** `g_i`, defined by

`family mean = mu + g_i + g_j + s_ij`.

Accordingly, the family-level optimization coefficient is

`g_i + g_j + s_ij`.

Do not use `0.5*(g_i+g_j)+s_ij` unless the input values called “GCA” are
actually additive breeding values `BV_i = 2*g_i`.

For diversity, input relationship matrices are additive relationship matrices
`A`. The R code forms the coancestry matrix `K = A/2`, computes
`Theta = p' K p`, and imposes `Theta <= 1/(2*Ns)`. Equivalently, using `A`
directly gives `p' A p <= 1/Ns`.

## Repository structure

```text
CMP-optimization-main/
├── README.md
├── LICENSE
├── 01_CMP_optimization/
│   ├── 001_R/
│   │   ├── CMP_optimization.R
│   │   └── input.xlsx
│   └── 002_Excel/
│       └── CMP_solver.xlsx
└── 02_manuscript_reproduction/
    └── CMP_manuscript_simulation_optimization.R
```

## 01_CMP_optimization: user-facing R workflow

Run from the repository root:

```bash
Rscript 01_CMP_optimization/001_R/CMP_optimization.R \
  01_CMP_optimization/001_R/input.xlsx 5 1000
```

Arguments are: (1) Excel input file, (2) target status number `Ns`, and
(3) requested total number of operational crosses.

The input workbook contains:

| Sheet | Required content |
|---|---|
| 1 | Standard GCA vector `g`, one value per parent |
| 2 | Symmetric SCA matrix |
| 3 | Additive relationship matrix `A` |
| 4 | Optional cross-specific upper-limit matrix |

The script solves the **continuous** family-allocation QCP. It writes the
continuous parent/family results plus an integer operational translation back
to the workbook. The integer translation is a deployment aid; any strict
integer, sex-specific, or direction-specific operational constraints should be
modeled explicitly if they must be guaranteed by the optimizer.

Required R packages are `Matrix`, `readxl`, `openxlsx`, and the Gurobi R API.
Gurobi requires a separate installation and license.

## 02_manuscript_reproduction: simulation and study scenarios

Run:

```bash
Rscript 02_manuscript_reproduction/CMP_manuscript_simulation_optimization.R
```

The corrected reproduction workflow uses:

- 50 independent stochastic population replicates;
- `h2 = 0.2, 0.5`;
- `Vd/Va = 0.25, 0.5, 1.0`;
- ordinary dominance QTL (`n.dominant`), not overdominance;
- target status numbers `Ns = 2, 5, 10, 15, 20`;
- one deterministic optimization per stochastic replicate × parameter setting × scenario × `Ns`.

Study scenarios:

| Scenario | Decision information used in optimization |
|---|---|
| `sc1_true` | True standard GCA + true SCA |
| `sc2_hs` | Half-sib GCA only |
| `sc3_30` | Full-sib GCA + SCA estimated from 30 progeny per cross |
| `sc4_100` | Full-sib GCA + SCA estimated from 100 progeny per cross |

The script records realized variance-component and heritability diagnostics so
the simulated architecture can be checked rather than inferred only from input
settings.

### Dependencies

CRAN/package dependencies include `dplyr`, `tidyr`, `Matrix`, `openxlsx`, and
`ggplot2`, plus MoBPS and its dependencies. ASReml-R and Gurobi are separately
licensed dependencies.

## Reproducibility note

Results, confidence intervals, tables, and figures in the manuscript should be
regenerated after changing GCA scaling, dominance architecture, residual
variance, or the replicate structure. The previous numerical results should not
be mixed with outputs from the corrected implementation.

## Contact
Milan Lstiburek: lstiburek@fld.czu.cz
