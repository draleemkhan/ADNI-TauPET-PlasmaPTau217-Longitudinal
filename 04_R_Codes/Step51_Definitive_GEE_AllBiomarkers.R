# ============================================================
# STEP 51
# DEFINITIVE GEE AUDIT — ALL BIOMARKER x TIME TERMS
#
# Purpose:
# Reproduce the exact Step29C GEE specifications and extract:
#   1. Tau PET x Time
#   2. Centiloid x Time
#   3. Plasma p-tau217 x Time
#
# Eight models:
#   - Independence / Exchangeable
#   - Unweighted
#   - IPCW untrimmed
#   - IPCW trimmed at 99th percentile
#   - IPCW trimmed at 95th percentile
#
# IMPORTANT:
# No model specification is changed from Step29C.
# ============================================================

rm(list = ls())

cat("\n============================================================\n")
cat(" STEP 51: DEFINITIVE GEE AUDIT — ALL BIOMARKERS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_file <- paste0(
  "D:/output/",
  "Step29C_Weighted_GEE_Input.csv"
)

output_dir <- "D:/output"

if (!requireNamespace("geepack", quietly = TRUE)) {
  stop('Package "geepack" is not installed.')
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

d <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input observations:", nrow(d), "\n")
cat("Unique participants:", length(unique(d$RID)), "\n\n")

stopifnot(nrow(d) == 379L)
stopifnot(length(unique(d$RID)) == 143L)

# geeglm requires contiguous observations within clusters
d <- d[order(d$RID, d$TimeYears), ]

d$RID <- as.integer(d$RID)
d$SexModel <- factor(d$SexModel)

cat("Dataset audit: PASS\n")

# ------------------------------------------------------------
# Exact Step29C model formula
# ------------------------------------------------------------

form <- CDRSB ~
  TimeYears +
  Tau_z +
  Centiloid_z +
  PTau217_z +
  Age +
  SexModel +
  Education +
  Tau_z:TimeYears +
  Centiloid_z:TimeYears +
  PTau217_z:TimeYears

cat("\nModel formula:\n")
print(form)

# ------------------------------------------------------------
# Exact Step29C fitting function
# ------------------------------------------------------------

fit_one <- function(corstr, weights = NULL) {

  args <- list(
    formula = form,
    id = d$RID,
    data = d,
    family = gaussian(link = "identity"),
    corstr = corstr,
    std.err = "san.se"
  )

  if (!is.null(weights)) {
    args$weights <- weights
  }

  do.call(geepack::geeglm, args)
}

# ------------------------------------------------------------
# Fit exact eight Step29C models
# ------------------------------------------------------------

cat("\nFitting eight GEE models...\n")

models <- list(

  Independence_Unweighted =
    fit_one("independence"),

  Independence_Weighted =
    fit_one("independence", d$Weight_Untrimmed),

  Independence_Weighted_Trim99 =
    fit_one("independence", d$Weight_Trim99),

  Independence_Weighted_Trim95 =
    fit_one("independence", d$Weight_Trim95),

  Exchangeable_Unweighted =
    fit_one("exchangeable"),

  Exchangeable_Weighted =
    fit_one("exchangeable", d$Weight_Untrimmed),

  Exchangeable_Weighted_Trim99 =
    fit_one("exchangeable", d$Weight_Trim99),

  Exchangeable_Weighted_Trim95 =
    fit_one("exchangeable", d$Weight_Trim95)
)

cat("Model fitting complete.\n")

# ============================================================
# Helper function to extract one interaction
# ============================================================

extract_term <- function(m, label, biomarker, possible_names) {

  co <- summary(m)$coefficients
  rn <- rownames(co)

  target <- which(rn %in% possible_names)

  if (length(target) != 1L) {
    stop(
      paste(
        "Could not uniquely identify",
        biomarker,
        "interaction in",
        label
      )
    )
  }

  b <- co[target, "Estimate"]

  se_col <- grep(
    "Std.err|Std\\.err",
    colnames(co),
    value = TRUE
  )[1]

  se <- co[target, se_col]

  z <- b / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  data.frame(
    Model = label,
    Biomarker = biomarker,
    Term = rn[target],
    Beta = b,
    RobustSE = se,
    CI_Lower = b - 1.96 * se,
    CI_Upper = b + 1.96 * se,
    Z = z,
    P = p,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Extract all three biomarker x time terms
# ------------------------------------------------------------

results_list <- list()

counter <- 1L

for (model_name in names(models)) {

  m <- models[[model_name]]

  results_list[[counter]] <- extract_term(
    m,
    model_name,
    "Tau PET",
    c(
      "TimeYears:Tau_z",
      "Tau_z:TimeYears"
    )
  )
  counter <- counter + 1L

  results_list[[counter]] <- extract_term(
    m,
    model_name,
    "Centiloid",
    c(
      "TimeYears:Centiloid_z",
      "Centiloid_z:TimeYears"
    )
  )
  counter <- counter + 1L

  results_list[[counter]] <- extract_term(
    m,
    model_name,
    "Plasma p-tau217",
    c(
      "TimeYears:PTau217_z",
      "PTau217_z:TimeYears"
    )
  )
  counter <- counter + 1L
}

AllResults <- do.call(rbind, results_list)
rownames(AllResults) <- NULL

# ------------------------------------------------------------
# Working-correlation estimates
# ------------------------------------------------------------

WorkingCorrelation <- data.frame(

  Model = names(models),

  Alpha = sapply(
    models,
    function(m) {

      a <- m$geese$alpha

      if (length(a) == 0) {
        NA_real_
      } else {
        as.numeric(a[1])
      }
    }
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Tau PET reproduction audit against Step29C
# ------------------------------------------------------------

old_file <- file.path(
  output_dir,
  "Step29C_Weighted_GEE_Retention_Results.csv"
)

if (!file.exists(old_file)) {
  stop("Original Step29C results file not found.")
}

old <- read.csv(
  old_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tau_new <- AllResults[
  AllResults$Biomarker == "Tau PET",
]

tau_new <- tau_new[
  match(old$Model, tau_new$Model),
]

tol <- 1e-8

beta_pass <- all(
  abs(old$TauTimeBeta - tau_new$Beta) < tol
)

se_pass <- all(
  abs(old$RobustSE - tau_new$RobustSE) < tol
)

lower_pass <- all(
  abs(old$CI_Lower - tau_new$CI_Lower) < tol
)

upper_pass <- all(
  abs(old$CI_Upper - tau_new$CI_Upper) < tol
)

p_pass <- all(
  abs(old$P - tau_new$P) < tol
)

reproduction_pass <- all(
  beta_pass,
  se_pass,
  lower_pass,
  upper_pass,
  p_pass
)

# ============================================================
# Output
# ============================================================

cat("\n============================================================\n")
cat(" ALL BIOMARKER x TIME RESULTS\n")
cat("============================================================\n\n")

print(
  AllResults,
  row.names = FALSE,
  digits = 6
)

cat("\n============================================================\n")
cat(" WORKING CORRELATION ESTIMATES\n")
cat("============================================================\n\n")

print(
  WorkingCorrelation,
  row.names = FALSE,
  digits = 6
)

cat("\n============================================================\n")
cat(" STEP29C TAU PET REPRODUCTION AUDIT\n")
cat("============================================================\n\n")

cat(
  "Tau PET beta reproduction: ",
  ifelse(beta_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat(
  "Tau PET SE reproduction:   ",
  ifelse(se_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat(
  "Tau PET lower CI:          ",
  ifelse(lower_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat(
  "Tau PET upper CI:          ",
  ifelse(upper_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat(
  "Tau PET P value:           ",
  ifelse(p_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat(
  "\nOVERALL STEP29C REPRODUCTION: ",
  ifelse(reproduction_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write.csv(
  AllResults,
  file.path(
    output_dir,
    "Step51_GEE_AllBiomarker_TimeInteractions.csv"
  ),
  row.names = FALSE
)

write.csv(
  WorkingCorrelation,
  file.path(
    output_dir,
    "Step51_GEE_WorkingCorrelation.csv"
  ),
  row.names = FALSE
)

saveRDS(
  models,
  file.path(
    output_dir,
    "Step51_GEE_Models.rds"
  )
)

cat("\n============================================================\n")
cat(" STEP 51 COMPLETE\n")
cat("============================================================\n")

cat("Observations: ", nrow(d), "\n", sep = "")
cat(
  "Participants: ",
  length(unique(d$RID)),
  "\n",
  sep = ""
)

cat(
  "Step29C Tau PET reproduction: ",
  ifelse(reproduction_pass, "PASS", "FAIL"),
  "\n",
  sep = ""
)

cat("\nSaved outputs to:\n")
cat(output_dir, "\n")

cat("\nFiles created:\n")
cat("  Step51_GEE_AllBiomarker_TimeInteractions.csv\n")
cat("  Step51_GEE_WorkingCorrelation.csv\n")
cat("  Step51_GEE_Models.rds\n")

cat("============================================================\n")