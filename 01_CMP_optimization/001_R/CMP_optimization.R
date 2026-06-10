# -----------------------------------------------------------------------------
# Script: Optimization of control mass pollination in seed orchards
# Authors: Christi Sagariya, Yousry A. El-Kassaby, Ye-Ji Kim, Milan Lstibůrek
# Purpose: Maximize genetic gain from GCA and SCA while constraining diversity
#          using a target status number (Ns).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(gurobi)
    library(Matrix)
    library(readxl)
    library(openxlsx)
})

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
# Usage from terminal:
#   Rscript R/cmp_optimization.R input.xlsx 5 1000
# Arguments:
#   1. Excel input file path
#   2. Target status number, Ns
#   3. Total number of operational crosses

args <- commandArgs(trailingOnly = TRUE)

filename <- ifelse(length(args) >= 1, args[1], "input.xlsx")
target_ns <- ifelse(length(args) >= 2, as.numeric(args[2]), 5)
number_crosses <- ifelse(length(args) >= 3, as.numeric(args[3]), 1000)

if (!file.exists(filename)) {
    stop(sprintf("Input file not found: %s", filename))
}

if (is.na(target_ns) || target_ns <= 0) {
    stop("target_ns must be a positive number.")
}

if (is.na(number_crosses) || number_crosses <= 0) {
    stop("number_crosses must be a positive number.")
}

# -----------------------------------------------------------------------------
# 1. Data input
# -----------------------------------------------------------------------------
# Expected Excel workbook structure:
#   Sheet 1: GCA vector, one column, n parents x 1
#   Sheet 2: SCA matrix, n parents x n parents
#   Sheet 3: Genomic relationship matrix G, n parents x n parents
#   Sheet 4: Optional cross-limit matrix, n parents x n parents
#            0 = banned cross, value between 0 and 1 = upper bound,
#            blank = no custom limit.

gca_raw <- read_excel(filename, sheet = 1, col_names = FALSE)
gca <- as.numeric(gca_raw[[1]])
n_parents <- length(gca)

sca_raw <- read_excel(filename, sheet = 2, col_names = FALSE)
sca <- as.matrix(sca_raw)

G_raw <- read_excel(filename, sheet = 3, col_names = FALSE)
G <- as.matrix(G_raw)

y_limits_raw <- tryCatch(
    read_excel(filename, sheet = 4, col_names = FALSE),
    error = function(e) {
        message("Warning: Sheet 4 with cross limits was not found. No custom cross limits will be applied.")
        NULL
    }
)

if (!is.null(y_limits_raw)) {
    y_limits <- as.matrix(y_limits_raw)
    if (nrow(y_limits) != n_parents || ncol(y_limits) != n_parents) {
        stop("The cross-limit matrix dimensions do not match the number of parents.")
    }
} else {
    y_limits <- NULL
}

# -----------------------------------------------------------------------------
# 2. Data validation and preprocessing
# -----------------------------------------------------------------------------
cat("Detected number of parents:", n_parents, "\n")

if (nrow(sca) != n_parents || ncol(sca) != n_parents) {
    stop("SCA matrix dimensions do not match the GCA vector length.")
}

if (nrow(G) != n_parents || ncol(G) != n_parents) {
    stop("G matrix dimensions do not match the GCA vector length.")
}

colnames(sca) <- rownames(sca) <- seq_len(n_parents)
colnames(G) <- rownames(G) <- seq_len(n_parents)
names(gca) <- seq_len(n_parents)

# Convert the genetic relationship matrix to a coancestry matrix.
G <- G / 2

# Force symmetry to protect the optimizer from small numerical errors.
sca <- (sca + t(sca)) / 2
G <- (G + t(G)) / 2

# Force positive semi-definiteness of the coancestry matrix.
eig <- eigen(G)
eig$values[eig$values < 1e-6] <- 1e-6
G <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

# -----------------------------------------------------------------------------
# 3. Build mapping matrices
# -----------------------------------------------------------------------------
# The decision variables are family proportions, Y_ij, for unique crosses i < j.

family_map <- combn(n_parents, 2)
n_families <- ncol(family_map)
cat("Number of decision variables, families:", n_families, "\n")

# Objective coefficient for each family:
# gain_ij = 0.5 * GCA_i + 0.5 * GCA_j + SCA_ij
obj_coeffs <- numeric(n_families)
for (k in seq_len(n_families)) {
    p1 <- family_map[1, k]
    p2 <- family_map[2, k]
    obj_coeffs[k] <- 0.5 * gca[p1] + 0.5 * gca[p2] + sca[p1, p2]
}

# M maps family proportions to parental contributions.
M <- Matrix(0, nrow = n_parents, ncol = n_families, sparse = TRUE)
for (k in seq_len(n_families)) {
    p1 <- family_map[1, k]
    p2 <- family_map[2, k]
    M[p1, k] <- 0.5
    M[p2, k] <- 0.5
}

# Diversity constraint matrix: Qc = M'GM.
Qc <- t(M) %*% G %*% M

# -----------------------------------------------------------------------------
# 4. Construct Gurobi model
# -----------------------------------------------------------------------------
model <- list()
model$obj <- obj_coeffs
model$modelsense <- "max"
model$varnames <- paste0("Cross_", family_map[1, ], "_", family_map[2, ])

# Sum of all family proportions equals 1.
model$A <- matrix(1, nrow = 1, ncol = n_families)
model$rhs <- 1
model$sense <- "="

# Bounds: 0 <= Y_ij <= 1 unless custom upper limits are supplied.
model$lb <- rep(0, n_families)
model$ub <- rep(1, n_families)

if (!is.null(y_limits)) {
    cat("Applying custom cross limits from Sheet 4...\n")
    count_constrained <- 0
    
    for (k in seq_len(n_families)) {
        p1 <- family_map[1, k]
        p2 <- family_map[2, k]
        
        val1 <- y_limits[p1, p2]
        val2 <- y_limits[p2, p1]
        lim1 <- ifelse(is.na(val1), 1, val1)
        lim2 <- ifelse(is.na(val2), 1, val2)
        limit <- min(lim1, lim2)
        
        if (limit < 1) {
            model$ub[k] <- limit
            count_constrained <- count_constrained + 1
        }
    }
    cat(sprintf("Applied custom limits to %d crosses.\n", count_constrained))
}

# Status number constraint. Coancestry <= 1 / (2 Ns).
diversity_limit <- 1 / (2 * target_ns)
model$quadcon <- list(
    list(
        Qc = as.matrix(Qc),
        rhs = diversity_limit,
        sense = "<="
    )
)

# -----------------------------------------------------------------------------
# 5. Solve
# -----------------------------------------------------------------------------
params <- list(OutputFlag = 1)
result <- gurobi(model, params)

if (result$status != "OPTIMAL") {
    stop(sprintf("Optimization failed. Gurobi status: %s", result$status))
}

Y_opt <- result$x
x_opt <- as.vector(M %*% Y_opt)
cat("\nOptimization successful.\n")

# -----------------------------------------------------------------------------
# 6. Post-optimal analysis
# -----------------------------------------------------------------------------
gain_gca_component <- sum(x_opt * gca)
gain_sca_component <- 0
active_idx <- which(Y_opt > 1e-8)

for (k in active_idx) {
    p1 <- family_map[1, k]
    p2 <- family_map[2, k]
    gain_sca_component <- gain_sca_component + Y_opt[k] * sca[p1, p2]
}

total_gain_calculated <- gain_gca_component + gain_sca_component

cat("\n=== RESULTS SUMMARY ===\n")
cat(sprintf("Total genetic gain:       %.5f\n", total_gain_calculated))
cat(sprintf("Contribution from GCA:    %.5f\n", gain_gca_component))
cat(sprintf("Contribution from SCA:    %.5f\n", gain_sca_component))
cat(sprintf("Target status number Ns:  %.2f\n", target_ns))
cat(sprintf("Diversity limit:          %.5f\n", diversity_limit))

# -----------------------------------------------------------------------------
# 7. Generate operational crossing plan
# -----------------------------------------------------------------------------
active_parents <- which(x_opt > 1e-6)
clone_plan <- data.frame(
    Clone_ID = active_parents,
    Optimized_p = round(x_opt[active_parents], 5),
    Target_Female_f = round(x_opt[active_parents], 5)
)

# -----------------------------------------------------------------------------
# 8. Export results to the Excel workbook
# -----------------------------------------------------------------------------
cat("Preparing Excel output...\n")

wb <- loadWorkbook(filename)

sn_p <- "Output_Parents_p"
sn_y <- "Output_Families_Y"
sn_ops <- "Operational_Plan_Simple"

for (sn in c(sn_p, sn_y, sn_ops)) {
    if (sn %in% names(wb)) removeWorksheet(wb, sn)
    addWorksheet(wb, sn)
}

# Operational plan sheet.
headers_ops <- c(
    "Clone ID", "Optimized Total (p)", "Target Female (f)",
    "Required Male (m)", "Status", "Femaleness_Index (Helper)"
)
writeData(wb, sheet = sn_ops, x = t(headers_ops), startCol = 1, startRow = 1, colNames = FALSE)
writeData(wb, sheet = sn_ops, x = clone_plan, startCol = 1, startRow = 2, colNames = FALSE)

n_rows_ops <- nrow(clone_plan)
for (i in seq_len(n_rows_ops)) {
    r <- i + 1
    writeFormula(wb, sn_ops, startCol = 4, startRow = r, x = sprintf("2*B%d - C%d", r, r))
    writeFormula(wb, sn_ops, startCol = 5, startRow = r,
                 x = sprintf("IF(OR(D%d < -1E-6, C%d < -1E-6), \"IMPOSSIBLE\", \"OK\")", r, r))
    writeFormula(wb, sn_ops, startCol = 6, startRow = r,
                 x = sprintf("IFERROR(C%d / (2*B%d), 0.5)", r, r))
}

style_header <- createStyle(textDecoration = "bold", border = "Bottom")
style_input <- createStyle(fgFill = "#FFFF00", border = "TopBottomLeftRight")
addStyle(wb, sn_ops, style_header, rows = 1, cols = 1:6)
if (n_rows_ops > 0) {
    addStyle(wb, sn_ops, style_input, rows = 2:(n_rows_ops + 1), cols = 3, gridExpand = TRUE)
}

# Family plan sheet.
df_Y <- data.frame(
    Parent1 = family_map[1, active_idx],
    Parent2 = family_map[2, active_idx],
    Proportion = round(Y_opt[active_idx], 5)
)
df_Y$Num_Crosses <- round(df_Y$Proportion * number_crosses)
df_Y <- df_Y[order(-df_Y$Num_Crosses), ]

headers_y <- c("Parent1", "Parent2", "Proportion", "Num_Crosses", "As_Female (P1)", "As_Male (P1)")
writeData(wb, sheet = sn_y, x = t(headers_y), startCol = 1, startRow = 1, colNames = FALSE)
writeData(wb, sheet = sn_y, x = df_Y, startCol = 1, startRow = 2, colNames = FALSE)

n_rows_y <- nrow(df_Y)
lookup_range <- sprintf("'%s'!$A$2:$F$%d", sn_ops, n_rows_ops + 1)

for (i in seq_len(n_rows_y)) {
    r <- i + 1
    f_p1_idx <- sprintf("VLOOKUP(A%d, %s, 6, FALSE)", r, lookup_range)
    f_p2_idx <- sprintf("VLOOKUP(B%d, %s, 6, FALSE)", r, lookup_range)
    ratio_calc <- sprintf("(%s) / ((%s) + (%s))", f_p1_idx, f_p1_idx, f_p2_idx)
    formula_E <- sprintf("IFERROR(ROUND(D%d * %s, 0), ROUND(D%d * 0.5, 0))", r, ratio_calc, r)
    formula_F <- sprintf("D%d - E%d", r, r)
    
    writeFormula(wb, sn_y, startCol = 5, startRow = r, x = formula_E)
    writeFormula(wb, sn_y, startCol = 6, startRow = r, x = formula_F)
}

last_row <- n_rows_y + 2
writeData(wb, sheet = sn_y, x = "TOTAL", startCol = 1, startRow = last_row)
writeFormula(wb, sn_y, startCol = 4, startRow = last_row, x = sprintf("SUM(D2:D%d)", last_row - 1))
writeFormula(wb, sn_y, startCol = 5, startRow = last_row, x = sprintf("SUM(E2:E%d)", last_row - 1))
writeFormula(wb, sn_y, startCol = 6, startRow = last_row, x = sprintf("SUM(F2:F%d)", last_row - 1))

style_total <- createStyle(textDecoration = "bold", border = "Top")
addStyle(wb, sn_y, style_header, rows = 1, cols = 1:6)
addStyle(wb, sn_y, style_total, rows = last_row, cols = 1:6)

# Parent contribution sheet.
df_p <- data.frame(Parent_ID = seq_len(n_parents), Contribution = x_opt)
writeData(wb, sheet = sn_p, x = df_p)
addStyle(wb, sn_p, style_header, rows = 1, cols = 1:2)

saveWorkbook(wb, filename, overwrite = TRUE)
cat(sprintf("Excel results saved to: %s\n", filename))
