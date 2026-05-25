# A novel convex  optimization framework for controlled mass pollination: capturing additive and non-additive genetic gain in seed orchards
This repository contains an R script for the article A novel convex  optimization framework for controlled mass pollination: capturing additive and non-additive genetic gain in seed orchards. We provide a practical framework for control mass pollination (CMP) optimization in seed orchards to maximise total genetic gain by incorporating both general combining ability (GCA) and specific combining ability (SCA), while constraining relatedness through a target status number \(N_s\).
## Part 1: R script for data generation via simulation and CMP optimization used in this article. 
## R script workflow
```text
R/cmp_sim_opt.R
```
The script simulates a breeding population and evaluates controlled mass pollination (CMP) decisions under four information scenarios:

1. **Scenario 1 (`sc1_true`)**: true GCA and true SCA used for optimization.
2. **Scenario 2 (`sc2_hs`)**: half-sib GCA only used for optimization.
3. **Scenario 3 (`sc3_30`)**: full-sib GCA and SCA estimated from 30 progeny per cross.
4. **Scenario 4 (`sc4_100`)**: full-sib GCA and SCA estimated from 100 progeny per cross.

For each scenario, the workflow optimizes cross contributions subject to a status-number diversity constraint and then evaluates the realized genetic gain using the true simulated GCA/SCA values.

## Requirements

Install R and the required R packages before running the workflow:

```r
install.packages(c(
  "dplyr", "tidyr", "Matrix", "reshape2", "openxlsx", "slam", "ggplot2"
))
```

The following packages require separate installation or licensing steps:

- **ASReml-R** (`asreml`): proprietary software; install according to your license.
- **Gurobi** (`gurobi`): install Gurobi Optimizer and configure a valid license.
- **MoBPS**, **miraculix**, and **RandomFieldsUtils**: install from their documented sources if not available through your standard R repository setup.

## Part 2: R script - CMP optimization for users. 
## Main file

```text
R/cmp_optimization.R
```

The script reads an Excel workbook, solves the optimization problem with Gurobi, and writes the optimized parental contributions and operational crossing plan back into the same workbook.

## Required software

Install R and the following R packages:

```r
install.packages(c("Matrix", "readxl", "openxlsx"))
```

The script also requires the `gurobi` R package, which is installed with the Gurobi Optimizer rather than from CRAN. Install Gurobi first, activate a license, and then install the R package from the Gurobi installation folder.

Example on Windows, after installing Gurobi:

```r
install.packages("C:/gurobi1200/win64/R/gurobi_12.0-0.zip", repos = NULL)
library(gurobi)
gurobi::gurobi_version()
```

Adjust the path and version number to match your local Gurobi installation.

## Input Excel workbook

Prepare one Excel workbook, for example:

```text
input.xlsx
```

The workbook must contain these sheets in this order:

| Sheet | Content | Dimension |
|---|---|---|
| 1 | GCA values | n parents x 1 |
| 2 | SCA matrix | n parents x n parents |
| 3 | Genomic relationship matrix, G | n parents x n parents |
| 4 | Optional cross-limit matrix | n parents x n parents |


## How to run

From the repository folder, run:

```bash
Rscript R/cmp_optimization.R input.xlsx 5 1000
```

Arguments are:

```text
1. input Excel file path
2. target status number Ns
3. total number of operational crosses
```

For example:

```bash
Rscript R/cmp_optimization.R data/input.xlsx 5 1000
```

If no arguments are provided, the defaults are:

```text
input file: input.xlsx
target Ns: 5
number of crosses: 1000
```

## Output sheets

The script writes three output sheets to the same Excel workbook:

| Output sheet | Description |
|---|---|
| `Output_Parents_p` | optimized total parental contribution for each parent |
| `Output_Families_Y` | optimized family proportions and number of crosses |
| `Operational_Plan_Simple` | editable operational plan with female and male contribution formulas |

In `Operational_Plan_Simple`, the column `Target Female (f)` is highlighted. You can edit this column to adjust female contributions. The workbook formulas then update the required male contributions and the split of crosses in `Output_Families_Y`.

## Contact
For any further details/information/inquiries, please contact:

Corresponding author: Prof. Milan Lstibůrek
E-mail: lstiburek@fld.czu.cz

First author: Christi Sagariya
E-mail: csagariya@gmail.com


