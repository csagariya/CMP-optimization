# A novel convex optimization framework for controlled mass pollination

This repository supports the article **“A novel convex optimization framework for controlled mass pollination: capturing additive and non-additive genetic gain in seed orchards.”**

The repository is organized so that folder names match the workflow. Most users should start with **01_CMP_optimization**, which contains the practical CMP model and example inputs. The manuscript-scale MoBPS simulation and scenario reproduction workflow is separated under **02_manuscript_reproduction**.

---

## Contents

1. [CMP optimization model for users](#01_cmp_optimization)
   - [R implementation](#r-implementation)
   - [Excel implementation](#excel-implementation)
2. [Manuscript reproduction](#02_manuscript_reproduction)
   - [Simulation data generation](#simulation-data-generation)
   - [CMP optimization scenarios](#cmp-optimization-scenarios)
3. [Software requirements](#software-requirements)
4. [Contact](#contact)

---

## Repository structure

```text
repository_root/
├── README.md
├── LICENSE
├── 01_CMP_optimization/
│   ├── 001_R/
│   │   
│   │   ├── run_CMP_optimization.R
│   │   ├── input.xlsx
│   │   
│   └── 002_Excel/
│       ├── CMP_solver.xlsx
│       
└── 02_manuscript_reproduction/
    ├── Simulation_optimization_script.R
    
```

---

## 01_CMP_optimization

This folder contains the main CMP optimization model for new users.

### R implementation

Main file:

```text
01_CMP_optimization/001_R/run_CMP_optimization.R
```

Example input workbook:

```text
01_CMP_optimization/001_R/input_example/input.xlsx
```

Run from the repository root:

```bash
Rscript 01_CMP_optimization/001_R/run_CMP_optimization.R 01_CMP_optimization/001_R/input_example/input.xlsx 5 1000
```

Arguments:

```text
1. Input Excel workbook path
2. Target status number, Ns
3. Total number of operational crosses
```

If no arguments are provided, the script uses:

```text
01_CMP_optimization/001_R/input_example/input.xlsx
Ns = 5
number of crosses = 1000
```

The script writes three output sheets back into the workbook:

| Output sheet | Meaning |
|---|---|
| `Output_Parents_p` | Optimized total parental contribution for each parent |
| `Output_Families_Y` | Optimized family proportions and operational number of crosses |
| `Operational_Plan_Simple` | Editable operational plan for female and male contributions |

### Excel implementation

Spreadsheet solver file:

```text
01_CMP_optimization/002_Excel/CMP_optimization.xlsx
```

Use this option when users prefer a spreadsheet-based CMP optimization workflow. See the `file_description.md` file in the folder for the workbook purpose, required inputs, and expected outputs.

---

## 02_manuscript_reproduction

This folder contains the manuscript-scale workflow.

### Simulation data generation

Main MoBPS simulation/data-generation script:

```text
02_manuscript_reproduction/001_simulation_data/MoBPS_scripts/MoBPS_CMP_simulation_and_data_generation.R
```

The script generates scenario input data such as additive relationship matrices, true GCA/SCA values, half-sib GCA/SCA estimates, and full-sib GCA/SCA estimates.

### CMP optimization scenarios

The uploaded repository contained one integrated manuscript workflow script. The scenario folder therefore contains a runner note pointing to the integrated script and documents the expected inputs and outputs for the four manuscript scenarios:

| Scenario | Decision information used in optimization |
|---|---|
| `sc1_true` | True GCA + true SCA |
| `sc2_hs` | Half-sib GCA only |
| `sc3_30` | Full-sib GCA + SCA estimated from 30 progeny per cross |
| `sc4_100` | Full-sib GCA + SCA estimated from 100 progeny per cross |

---

## Software requirements

### User CMP optimization script

Install R packages:

```r
install.packages(c("Matrix", "readxl", "openxlsx"))
```

The script also requires the `gurobi` R package, which is installed with Gurobi Optimizer rather than from CRAN.

### Manuscript reproduction workflow

Install or configure:

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
