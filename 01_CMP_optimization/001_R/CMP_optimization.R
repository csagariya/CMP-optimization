# -----------------------------------------------------------------------------
# Script: CMP optimization for seed orchards (user-facing implementation)
# Purpose: Maximize expected full-sib family value from standard GCA + SCA
#          while constraining group coancestry / status number.
#
# IMPORTANT SCALING CONVENTION
# ----------------------------
# GCA values in Sheet 1 are STANDARD diallel GCA effects (g_i), i.e. the
# parental effect appearing in
#
#     family mean = mu + g_i + g_j + s_ij .
#
# Therefore the family objective coefficient is
#
#     g_i + g_j + s_ij
#
# and NOT 0.5 * (g_i + g_j) + s_ij.
# If your input values are additive breeding values (BV_i = 2*g_i), divide
# them by 2 before using this script, or replace the family additive term by
# 0.5*(BV_i + BV_j).
#
# RELATIONSHIP / STATUS-NUMBER CONVENTION
# ---------------------------------------
# Sheet 3 contains an additive relationship matrix A (pedigree or genomic).
# The script forms the coancestry matrix K = A/2 and constrains
#
#     Theta = p' K p <= 1/(2*Ns),
#
# where p is the vector of total parental genetic contributions and
# Ns = 0.5/Theta.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(gurobi)
    library(Matrix)
    library(readxl)
    library(openxlsx)
})

# -----------------------------------------------------------------------------
# 1. User settings / command-line arguments
# -----------------------------------------------------------------------------
# Example:
#   Rscript CMP_optimization.R input.xlsx 5 1000
#
# Arguments:
#   1. Excel input file
#   2. Target status number Ns
#   3. Total number of operational pollinations/crosses (integer)

args <- commandArgs(trailingOnly = TRUE)

filename <- if (length(args) >= 1) args[1] else "input.xlsx"
target_ns <- if (length(args) >= 2) as.numeric(args[2]) else 5
number_crosses <- if (length(args) >= 3) as.integer(round(as.numeric(args[3]))) else 1000L

if (!file.exists(filename)) {
    stop(sprintf("Input file not found: %s", filename), call. = FALSE)
}
if (!is.finite(target_ns) || target_ns <= 0) {
    stop("target_ns must be a positive finite number.", call. = FALSE)
}
if (is.na(number_crosses) || number_crosses <= 0L) {
    stop("number_crosses must be a positive integer.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Read workbook
# -----------------------------------------------------------------------------
# Sheet 1: standard GCA vector g_i, n x 1
# Sheet 2: SCA matrix s_ij, n x n
# Sheet 3: additive relationship matrix A, n x n
# Sheet 4: optional symmetric cross-limit matrix, n x n
#          0 = banned; (0,1] = maximum family proportion; blank = no limit

gca_raw <- read_excel(filename, sheet = 1, col_names = FALSE)
gca <- suppressWarnings(as.numeric(gca_raw[[1]]))

if (length(gca) < 2L || any(!is.finite(gca))) {
    stop("Sheet 1 must contain at least two finite numeric GCA values.", call. = FALSE)
}
n_parents <- length(gca)
names(gca) <- as.character(seq_len(n_parents))

sca <- as.matrix(read_excel(filename, sheet = 2, col_names = FALSE))
A <- as.matrix(read_excel(filename, sheet = 3, col_names = FALSE))

storage.mode(sca) <- "numeric"
storage.mode(A) <- "numeric"

if (!all(dim(sca) == c(n_parents, n_parents))) {
    stop("Sheet 2 (SCA) dimensions must be n parents x n parents.", call. = FALSE)
}
if (!all(dim(A) == c(n_parents, n_parents))) {
    stop("Sheet 3 (additive relationship matrix A) dimensions must be n x n.", call. = FALSE)
}
if (any(!is.finite(sca))) {
    stop("Sheet 2 contains missing or non-finite SCA values.", call. = FALSE)
}
if (any(!is.finite(A))) {
    stop("Sheet 3 contains missing or non-finite relationship values.", call. = FALSE)
}

# Optional cross-specific limits.
y_limits_raw <- tryCatch(
    read_excel(filename, sheet = 4, col_names = FALSE),
    error = function(e) NULL
)

if (!is.null(y_limits_raw)) {
    y_limits <- as.matrix(y_limits_raw)
    storage.mode(y_limits) <- "numeric"
    if (!all(dim(y_limits) == c(n_parents, n_parents))) {
        stop("Sheet 4 cross-limit matrix dimensions must be n x n.", call. = FALSE)
    }
    finite_limits <- y_limits[is.finite(y_limits)]
    if (length(finite_limits) > 0L && any(finite_limits < 0 | finite_limits > 1)) {
        stop("Finite Sheet 4 limits must lie between 0 and 1.", call. = FALSE)
    }
} else {
    y_limits <- NULL
    message("Sheet 4 not found: no custom cross-specific upper limits will be used.")
}

# -----------------------------------------------------------------------------
# 3. Preprocess matrices
# -----------------------------------------------------------------------------
# SCA is treated as reciprocal/symmetric in this formulation.
sca <- (sca + t(sca)) / 2
diag(sca) <- 0

# Symmetrize additive relationship matrix A.
A <- (A + t(A)) / 2

# Convert additive relationship A to coancestry K.
K <- A / 2

# Numerical PSD check. Small negative eigenvalues can occur because of rounding;
# materially negative eigenvalues indicate an invalid relationship matrix.
eig <- eigen(K, symmetric = TRUE)
min_eig <- min(eig$values)
if (min_eig < -1e-7) {
    stop(sprintf(
        "The coancestry matrix K=A/2 is not positive semidefinite (minimum eigenvalue %.6g). Check Sheet 3.",
        min_eig
    ), call. = FALSE)
}
if (min_eig < 0) {
    eig$values[eig$values < 0] <- 0
    K <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
    K <- (K + t(K)) / 2
}

dimnames(sca) <- list(names(gca), names(gca))
dimnames(K) <- list(names(gca), names(gca))

# -----------------------------------------------------------------------------
# 4. Family variables and objective
# -----------------------------------------------------------------------------
# Unique unordered families only: i < j; selfing is excluded by construction.
family_map <- combn(n_parents, 2)
n_families <- ncol(family_map)

# M maps family proportions y to total parental genetic contributions p:
# p_i = 0.5 * sum_{j != i} y_ij.
M <- Matrix(0, nrow = n_parents, ncol = n_families, sparse = TRUE)

cross_gca <- numeric(n_families)
cross_sca <- numeric(n_families)
obj_coeffs <- numeric(n_families)

for (k in seq_len(n_families)) {
    i <- family_map[1, k]
    j <- family_map[2, k]

    M[i, k] <- 0.5
    M[j, k] <- 0.5

    # Standard diallel decomposition: mu + g_i + g_j + s_ij.
    cross_gca[k] <- gca[i] + gca[j]
    cross_sca[k] <- sca[i, j]
    obj_coeffs[k] <- cross_gca[k] + cross_sca[k]
}

# Theta = p' K p = y' (M' K M) y.
Qc <- as.matrix(t(M) %*% K %*% M)
Qc <- (Qc + t(Qc)) / 2

# -----------------------------------------------------------------------------
# 5. Gurobi model
# -----------------------------------------------------------------------------
model <- list()
model$obj <- obj_coeffs
model$modelsense <- "max"
model$varnames <- paste0("Cross_", family_map[1, ], "_", family_map[2, ])

# All family proportions sum to one.
model$A <- matrix(1, nrow = 1, ncol = n_families)
model$rhs <- 1
model$sense <- "="
model$lb <- rep(0, n_families)
model$ub <- rep(1, n_families)

# Optional unordered family upper bounds. If the two reciprocal cells differ,
# the more restrictive bound is used.
if (!is.null(y_limits)) {
    for (k in seq_len(n_families)) {
        i <- family_map[1, k]
        j <- family_map[2, k]
        a <- y_limits[i, j]
        b <- y_limits[j, i]
        lim_a <- if (is.na(a)) 1 else a
        lim_b <- if (is.na(b)) 1 else b
        model$ub[k] <- min(lim_a, lim_b)
    }
}

# Status-number constraint: Theta <= 1/(2 Ns).
diversity_limit <- 1 / (2 * target_ns)
model$quadcon <- list(list(
    Qc = Qc,
    rhs = diversity_limit,
    sense = "<="
))

params <- list(OutputFlag = 1)
result <- gurobi(model, params)

if (is.null(result$status) || result$status != "OPTIMAL") {
    stop(sprintf(
        "Optimization did not return OPTIMAL status. Gurobi status: %s",
        ifelse(is.null(result$status), "NULL", result$status)
    ), call. = FALSE)
}

# -----------------------------------------------------------------------------
# 6. Continuous optimum and consistency checks
# -----------------------------------------------------------------------------
Y_opt <- result$x
Y_opt[abs(Y_opt) < 1e-12] <- 0
p_opt <- as.vector(M %*% Y_opt)

active_idx <- which(Y_opt > 1e-8)

gain_gca <- sum(Y_opt * cross_gca)
gain_sca <- sum(Y_opt * cross_sca)
gain_total <- gain_gca + gain_sca

# Same additive component written through parental contributions:
# sum y_ij(g_i+g_j) = 2 * sum p_i g_i.
gain_gca_check <- 2 * sum(p_opt * gca)
if (abs(gain_gca - gain_gca_check) > 1e-7) {
    stop("Internal GCA scaling check failed.", call. = FALSE)
}

theta <- as.numeric(t(p_opt) %*% K %*% p_opt)
status_number_realized <- 0.5 / theta

cat("\n=== CONTINUOUS OPTIMUM ===\n")
cat(sprintf("Total expected genetic response: %.8f\n", gain_total))
cat(sprintf("GCA component:                  %.8f\n", gain_gca))
cat(sprintf("SCA component:                  %.8f\n", gain_sca))
cat(sprintf("Target status number Ns:        %.4f\n", target_ns))
cat(sprintf("Realized group coancestry:      %.8f\n", theta))
cat(sprintf("Realized status number:         %.4f\n", status_number_realized))
cat(sprintf("Active families:                %d of %d\n", length(active_idx), n_families))

# -----------------------------------------------------------------------------
# 7. Convert continuous proportions to an exact integer crossing plan
# -----------------------------------------------------------------------------
# Largest-remainder allocation guarantees that the integer family counts sum
# exactly to number_crosses.
raw_counts <- Y_opt * number_crosses
int_counts <- floor(raw_counts + 1e-12)
remaining <- number_crosses - sum(int_counts)

if (remaining > 0L) {
    frac <- raw_counts - int_counts
    add_idx <- order(frac, decreasing = TRUE)[seq_len(remaining)]
    int_counts[add_idx] <- int_counts[add_idx] + 1L
} else if (remaining < 0L) {
    frac <- raw_counts - int_counts
    remove_candidates <- which(int_counts > 0L)
    remove_idx <- remove_candidates[order(frac[remove_candidates], decreasing = FALSE)[seq_len(abs(remaining))]]
    int_counts[remove_idx] <- int_counts[remove_idx] - 1L
}

if (sum(int_counts) != number_crosses) {
    stop("Integerization failed to preserve the requested number of crosses.", call. = FALSE)
}

# Check the rounded operational plan against the continuous constraints. The
# QCP itself is continuous, so integer rounding can move the plan slightly.
Y_integer <- int_counts / number_crosses
p_integer <- as.vector(M %*% Y_integer)
theta_integer <- as.numeric(t(p_integer) %*% K %*% p_integer)
status_number_integer <- 0.5 / theta_integer
max_ub_excess <- max(Y_integer - model$ub)

if (theta_integer > diversity_limit + 1e-10) {
    warning(sprintf(
        paste0("Integer rounding gives Theta=%.8f, above the continuous limit %.8f. ",
               "Use more operational crosses or solve an integer-constrained deployment model if the rounded plan must satisfy Ns exactly."),
        theta_integer, diversity_limit
    ))
}
if (max_ub_excess > 1e-10) {
    warning(sprintf(
        paste0("Integer rounding exceeds at least one cross upper bound by %.8g in proportion units. ",
               "Use more operational crosses or an integer-constrained deployment model for hard discrete limits."),
        max_ub_excess
    ))
}

integer_active <- which(int_counts > 0L)

# Reciprocal orientation: split each unordered family approximately 50:50.
# This is an operational translation only; direction-specific biological limits
# are not part of the continuous QCP unless explicitly modeled separately.
p1_as_female <- int_counts[integer_active] %/% 2L
p2_as_female <- int_counts[integer_active] - p1_as_female

family_plan <- data.frame(
    Parent1 = family_map[1, integer_active],
    Parent2 = family_map[2, integer_active],
    Optimized_Proportion = Y_opt[integer_active],
    Num_Crosses = int_counts[integer_active],
    P1_as_Female_P2_as_Male = p1_as_female,
    P2_as_Female_P1_as_Male = p2_as_female,
    GCA_Component = cross_gca[integer_active],
    SCA_Component = cross_sca[integer_active],
    Family_Value = obj_coeffs[integer_active],
    stringsAsFactors = FALSE
)
family_plan <- family_plan[order(-family_plan$Num_Crosses, -family_plan$Optimized_Proportion), ]

# Derive actual integer female/male use by parent.
female_counts <- integer(n_parents)
male_counts <- integer(n_parents)

for (r in seq_len(nrow(family_plan))) {
    p1 <- family_plan$Parent1[r]
    p2 <- family_plan$Parent2[r]
    n12 <- family_plan$P1_as_Female_P2_as_Male[r]
    n21 <- family_plan$P2_as_Female_P1_as_Male[r]

    female_counts[p1] <- female_counts[p1] + n12
    male_counts[p2] <- male_counts[p2] + n12
    female_counts[p2] <- female_counts[p2] + n21
    male_counts[p1] <- male_counts[p1] + n21
}

parent_plan <- data.frame(
    Parent_ID = seq_len(n_parents),
    Optimized_p = p_opt,
    Expected_Total_Appearances = 2 * p_opt * number_crosses,
    Integer_Female_Uses = female_counts,
    Integer_Male_Uses = male_counts,
    Integer_Total_Uses = female_counts + male_counts,
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 8. Export to the same workbook
# -----------------------------------------------------------------------------
wb <- loadWorkbook(filename)

sn_p <- "Output_Parents_p"
sn_y <- "Output_Families_Y"
sn_ops <- "Operational_Plan_Simple"
sn_sum <- "Optimization_Summary"

for (sn in c(sn_p, sn_y, sn_ops, sn_sum)) {
    if (sn %in% names(wb)) removeWorksheet(wb, sn)
    addWorksheet(wb, sn)
}

summary_out <- data.frame(
    Metric = c(
        "Gurobi status",
        "Objective value",
        "GCA component",
        "SCA component",
        "Target Ns",
        "Realized group coancestry Theta",
        "Realized Ns",
        "Requested operational crosses",
        "Integer-plan group coancestry Theta",
        "Integer-plan realized Ns",
        "Integer-plan maximum cross-bound excess",
        "Active continuous families",
        "Active integer families"
    ),
    Value = c(
        result$status,
        gain_total,
        gain_gca,
        gain_sca,
        target_ns,
        theta,
        status_number_realized,
        number_crosses,
        theta_integer,
        status_number_integer,
        max_ub_excess,
        length(active_idx),
        length(integer_active)
    ),
    stringsAsFactors = FALSE
)

writeData(wb, sn_sum, summary_out)
writeData(wb, sn_p, parent_plan)
writeData(wb, sn_y, data.frame(
    Parent1 = family_map[1, active_idx],
    Parent2 = family_map[2, active_idx],
    Proportion = Y_opt[active_idx],
    GCA_Component = cross_gca[active_idx],
    SCA_Component = cross_sca[active_idx],
    Family_Value = obj_coeffs[active_idx]
))
writeData(wb, sn_ops, family_plan)

header_style <- createStyle(textDecoration = "bold", border = "Bottom")
sheet_ncols <- c(
    Optimization_Summary = ncol(summary_out),
    Output_Parents_p = ncol(parent_plan),
    Output_Families_Y = 6L,
    Operational_Plan_Simple = ncol(family_plan)
)
for (sn in names(sheet_ncols)) {
    nc <- sheet_ncols[[sn]]
    addStyle(wb, sn, header_style, rows = 1, cols = seq_len(nc), gridExpand = TRUE)
    freezePane(wb, sn, firstRow = TRUE)
    setColWidths(wb, sn, cols = seq_len(nc), widths = "auto")
}

saveWorkbook(wb, filename, overwrite = TRUE)
cat(sprintf("\nResults written to: %s\n", filename))
