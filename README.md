# A novel convex optimization framework for controlled mass pollination

This repository supports the article **“A novel convex optimization framework for controlled mass pollination: capturing additive and non-additive genetic gain in seed orchards.”**

The repository is organized into two main workflows:

- **`01_CMP_optimization/`**: user-facing controlled mass pollination (CMP) optimization tools in R and Excel.
- **`02_manuscript_reproduction/`**: manuscript-scale simulation and optimization workflow used to reproduce the study scenarios.

---

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
    └── Simulation_optimization_script
```

---

## 01_CMP_optimization

This folder contains the main CMP optimization tools for users.

### R implementation

Main script:

```text
01_CMP_optimization/001_R/CMP_optimization.R
```

Example input workbook:

```text
01_CMP_optimization/001_R/input.xlsx
```

Run the script from the repository root with the bundled input workbook:

```bash
Rscript 01_CMP_optimization/001_R/CMP_optimization.R 01_CMP_optimization/001_R/input.xlsx 5 1000
```

Command-line arguments:

```text
1. Input Excel workbook path
2. Target status number, Ns
3. Total number of operational crosses
```

The script can also be run without command-line arguments. In that case, it looks for `input.xlsx` in the current working directory and uses the default settings `Ns = 5` and `number of crosses = 1000`:

```bash
cd 01_CMP_optimization/001_R
Rscript CMP_optimization.R
```

Expected input workbook structure:

| Sheet | Required content |
|---|---|
| Sheet 1 | GCA vector, one column, `n` parents × 1 |
| Sheet 2 | SCA matrix, `n` parents × `n` parents |
| Sheet 3 | Genomic relationship matrix, `n` parents × `n` parents |
| Sheet 4 | Optional cross-limit matrix; `0` = banned cross, values between `0` and `1` = upper bound, blank = no custom limit |

The script writes three output sheets back into the same workbook:

| Output sheet | Description |
|---|---|
| `Output_Parents_p` | Optimized total parental contribution for each parent |
| `Output_Families_Y` | Optimized family proportions and operational number of crosses |
| `Operational_Plan_Simple` | Editable operational plan for female and male contributions |

### Excel implementation

Spreadsheet solver file:

```text
01_CMP_optimization/002_Excel/CMP_solver.xlsx
```

Use this option for a spreadsheet-based CMP optimization workflow. Open `CMP_solver.xlsx` and follow the instructions inside the workbook to enter inputs and run the solver.

---

## 02_manuscript_reproduction

This folder contains the manuscript-scale simulation and scenario reproduction workflow.

Main script:

```text
02_manuscript_reproduction/Simulation_optimization_script
```

Run from the manuscript reproduction folder:

```bash
cd 02_manuscript_reproduction
Rscript Simulation_optimization_script
```

This workflow generates scenario input data, including additive relationship matrices, true GCA/SCA values, half-sib GCA/SCA estimates, and full-sib GCA/SCA estimates. The same script performs CMP optimization for the study scenarios below.

| Scenario | Decision information used in optimization |
|---|---|
| `sc1_true` | True GCA + true SCA |
| `sc2_hs` | Half-sib GCA only |
| `sc3_30` | Full-sib GCA + SCA estimated from 30 progeny per cross |
| `sc4_100` | Full-sib GCA + SCA estimated from 100 progeny per cross |

This script is computationally intensive. For testing, reduce simulation-scale parameters such as `n_iterations`, `nf`, `snp`, `nr`, `reps_vec`, and `itr_grid` before running the full manuscript-scale workflow.

---

## Software requirements

### R implementation for CMP optimization

Install CRAN packages:

```r
install.packages(c("Matrix", "readxl", "openxlsx"))
```

The script also requires the `gurobi` R package, which is installed with Gurobi Optimizer rather than from CRAN. A valid Gurobi installation and license are required.

### Excel implementation

The Excel workbook requires Microsoft Excel with Solver support enabled.

### Manuscript reproduction workflow

Install CRAN packages:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "Matrix",
  "reshape2",
  "openxlsx",
  "slam",
  "ggplot2"
))
```

Additional dependencies require separate installation or licensing:

- `asreml`: ASReml-R, proprietary software.
- `gurobi`: Gurobi Optimizer and R package, license required.
- `MoBPS`, `miraculix`, and `RandomFieldsUtils`: install according to the MoBPS package guidance.

---

## Contact

Corresponding author: Prof. Milan Lstibůrek  
E-mail: lstiburek@fld.czu.cz

First author: Christi Sagariya  
E-mail: csagariya@gmail.com
