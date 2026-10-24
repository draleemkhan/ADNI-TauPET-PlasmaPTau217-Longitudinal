# Step29C_Weighted_GEE_Retention_Sensitivity.R
# Weighted GEE with participant-cluster sandwich inference.
#
# Before first use only:
# install.packages("geepack")

input_file <- "D:/output/Step29C_Weighted_GEE_Input.csv"
output_dir <- "D:/output"

if (!requireNamespace("geepack", quietly = TRUE)) {
  stop('Package "geepack" is not installed. Run install.packages("geepack") once, then rerun.')
}

d <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(nrow(d) == 379L)
stopifnot(length(unique(d$RID)) == 143L)

# geeglm requires observations from each cluster to be contiguous.
d <- d[order(d$RID, d$TimeYears), ]

d$RID <- as.integer(d$RID)
d$SexModel <- factor(d$SexModel)

form <- CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z +
  Age + SexModel + Education +
  Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears

fit_one <- function(corstr, weights = NULL) {
  args <- list(
    formula = form,
    id = d$RID,
    data = d,
    family = gaussian(link = "identity"),
    corstr = corstr,
    std.err = "san.se"
  )
  if (!is.null(weights)) args$weights <- weights
  do.call(geepack::geeglm, args)
}

models <- list(
  Independence_Unweighted = fit_one("independence"),
  Independence_Weighted = fit_one("independence", d$Weight_Untrimmed),
  Independence_Weighted_Trim99 = fit_one("independence", d$Weight_Trim99),
  Independence_Weighted_Trim95 = fit_one("independence", d$Weight_Trim95),
  Exchangeable_Unweighted = fit_one("exchangeable"),
  Exchangeable_Weighted = fit_one("exchangeable", d$Weight_Untrimmed),
  Exchangeable_Weighted_Trim99 = fit_one("exchangeable", d$Weight_Trim99),
  Exchangeable_Weighted_Trim95 = fit_one("exchangeable", d$Weight_Trim95)
)

extract_tau <- function(m, label) {
  co <- summary(m)$coefficients
  rn <- rownames(co)
  target <- which(rn %in% c("TimeYears:Tau_z", "Tau_z:TimeYears"))
  if (length(target) != 1L) stop(paste("Tau x time term not found in", label))
  b <- co[target, "Estimate"]
  se_col <- grep("Std.err|Std\\.err", colnames(co), value = TRUE)[1]
  se <- co[target, se_col]
  z <- b / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  data.frame(
    Model = label,
    TauTimeBeta = b,
    RobustSE = se,
    CI_Lower = b - 1.96 * se,
    CI_Upper = b + 1.96 * se,
    P = p,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, Map(extract_tau, models, names(models)))
rownames(res) <- NULL

write.csv(res,
          file.path(output_dir, "Step29C_Weighted_GEE_Retention_Results.csv"),
          row.names = FALSE)

# Working-correlation estimates where applicable
alpha <- data.frame(
  Model = names(models),
  Alpha = sapply(models, function(m) {
    a <- m$geese$alpha
    if (length(a) == 0) NA_real_ else as.numeric(a[1])
  })
)
write.csv(alpha,
          file.path(output_dir, "Step29C_Weighted_GEE_WorkingCorrelation.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat(" STEP 29C: WEIGHTED GEE RETENTION SENSITIVITY\n")
cat("============================================================\n")
cat("Gaussian identity-link GEE; participant clustering; sandwich SEs.\n")
cat("One stabilized participant weight applied to all observations per RID.\n\n")

print(res, row.names = FALSE, digits = 6)

cat("\n============================================================\n")
cat(" WORKING CORRELATION ESTIMATES\n")
cat("============================================================\n")
print(alpha, row.names = FALSE, digits = 6)

cat("\nSaved results to:\n", output_dir, "\n", sep = "")
