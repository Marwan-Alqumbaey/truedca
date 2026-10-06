# Example: incident acute kidney injury (AKI) in the eICU-CRD demo (v2.0.1)
# Reference outcome: KDIGO serum creatinine criteria during the first 7 ICU days.
# Recorded outcome: AKI diagnosis entered during the ICU stay (ICD-9 584 / ICD-10 N17).
# Usage: Rscript analysis/eicu_aki.R <path to eicu demo csv folder>
suppressPackageStartupMessages(library(truedca))
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "../data/eicu-demo"
rd <- function(f) read.csv(gzfile(file.path(dir, paste0(f, ".csv.gz"))), stringsAsFactors = FALSE)
set.seed(2026)
out <- "analysis/results"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
flow <- list()

# ---- cohort -----------------------------------------------------------------
pat <- rd("patient")
pat$age_n <- suppressWarnings(as.numeric(ifelse(pat$age == "> 89", "90", pat$age)))
flow$all_stays <- nrow(pat)
pat <- pat[pat$unitvisitnumber == 1, ]
pat <- pat[order(pat$uniquepid, pat$hospitaldischargeyear, pat$patientunitstayid), ]
pat <- pat[!duplicated(pat$uniquepid), ]
flow$first_icu_stay <- nrow(pat)
pat <- pat[!is.na(pat$age_n) & pat$age_n >= 18, ]
flow$adult <- nrow(pat)
pat <- pat[pat$unitdischargeoffset >= 1440, ]
flow$stay_24h <- nrow(pat)

lab <- rd("lab")
cr <- lab[lab$labname == "creatinine", c("patientunitstayid", "labresultoffset", "labresult")]
cr$labresult <- suppressWarnings(as.numeric(cr$labresult))
cr <- cr[!is.na(cr$labresult) & cr$labresult > 0 & cr$labresultoffset >= -1440 &
         cr$labresultoffset <= 10080 & cr$patientunitstayid %in% pat$patientunitstayid, ]
cr <- cr[order(cr$patientunitstayid, cr$labresultoffset), ]

# KDIGO creatinine criterion: rise >= 0.3 mg/dL within 48 h, or >= 1.5 x the lowest value in the prior 7 days
kdigo <- function(d, from, to) {
  idx <- which(d$labresultoffset > from & d$labresultoffset <= to)
  any(vapply(idx, function(j) {
    s <- d$labresultoffset[j]; v <- d$labresult[j]
    p48 <- d$labresult[d$labresultoffset >= s - 2880 & d$labresultoffset < s]
    p7 <- d$labresult[d$labresultoffset >= s - 10080 & d$labresultoffset < s]
    (length(p48) && v - min(p48) >= 0.3) || (length(p7) && v / min(p7) >= 1.5)
  }, logical(1)))
}
crs <- split(cr, cr$patientunitstayid)
info <- do.call(rbind, lapply(crs, function(d) data.frame(
  patientunitstayid = d$patientunitstayid[1],
  n_cr = sum(d$labresultoffset >= 0),
  n_all = nrow(d),
  truth = kdigo(d, 0, 10080))))
pat <- merge(pat, info, by = "patientunitstayid")
pat <- pat[pat$n_all >= 2 & pat$n_cr >= 1, ]
flow$two_creatinine <- nrow(pat)

dx <- rd("diagnosis")
aki_dx <- grepl("acute renal failure|acute kidney injury", dx$diagnosisstring, ignore.case = TRUE) |
          grepl("584|N17", dx$icd9code)
esrd_dx <- grepl("ESRD|end stage renal|chronic kidney disease.*(stage 5|dialysis)", dx$diagnosisstring, ignore.case = TRUE)
ph <- rd("pastHistory")
dialysis_hx <- unique(ph$patientunitstayid[grepl("dialysis", ph$pasthistorypath, ignore.case = TRUE)])
aps <- rd("apacheApsVar"); pred <- rd("apachePredVar")
pat <- pat[!(pat$patientunitstayid %in% c(dx$patientunitstayid[esrd_dx], dialysis_hx,
                                          aps$patientunitstayid[aps$dialysis == 1])), ]
flow$no_esrd_dialysis <- nrow(pat)
coded <- dx[aki_dx, c("patientunitstayid", "diagnosisoffset")]
coded <- merge(coded, pat[, c("patientunitstayid", "unitdischargeoffset")])
pat$y <- as.integer(pat$patientunitstayid %in%
                    coded$patientunitstayid[coded$diagnosisoffset <= coded$unitdischargeoffset])
pat$truth <- as.integer(pat$truth)

# ---- predictors from day 1 (creatinine, urea and urine output excluded) -----
d <- merge(pat, aps[, c("patientunitstayid", "vent", "meanbp", "heartrate", "respiratoryrate",
                        "temperature", "wbc", "sodium", "glucose", "hematocrit")], all.x = TRUE)
d <- merge(d, pred[, c("patientunitstayid", "diabetes", "cirrhosis", "immunosuppression",
                       "metastaticcancer", "electivesurgery")], all.x = TRUE)
num <- c("meanbp", "heartrate", "respiratoryrate", "temperature", "wbc", "sodium", "glucose", "hematocrit")
for (v in num) { x <- d[[v]]; x[is.na(x) | x < 0] <- NA; d[[v]] <- ifelse(is.na(x), median(x, na.rm = TRUE), x) }
for (v in c("vent", "diabetes", "cirrhosis", "immunosuppression", "metastaticcancer", "electivesurgery"))
  d[[v]] <- ifelse(is.na(d[[v]]) | d[[v]] < 0, 0, d[[v]])
d$male <- as.integer(d$gender == "Male")
d$source <- factor(ifelse(grepl("Emergency", d$unitadmitsource), "ED",
                   ifelse(grepl("Operating|Recovery|PACU", d$unitadmitsource), "Surgery", "Other")))
fml <- y ~ age_n + male + source + vent + meanbp + heartrate + respiratoryrate + temperature +
  wbc + sodium + glucose + hematocrit + diabetes + cirrhosis + immunosuppression + metastaticcancer + electivesurgery

# 10-fold cross-fitted risks; model trained on the recorded outcome, as is common practice
fold <- sample(rep(1:10, length.out = nrow(d)))
d$risk <- NA_real_
for (k in 1:10) {
  fit <- suppressWarnings(glm(fml, binomial, data = d[fold != k, ]))
  d$risk[fold == k] <- predict(fit, d[fold == k, ], type = "response")
}
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }

# ---- descriptive accuracy of the recorded outcome ---------------------------
grid <- seq(0.02, 0.30, by = 0.01)
rate <- function(a, b) if (sum(b)) mean(a[b]) else NA_real_
acc <- do.call(rbind, lapply(c(0.05, 0.10, 0.15, 0.20), function(p) {
  f <- d$risk >= p
  data.frame(threshold = p, flagged = mean(f),
    se1 = rate(d$y == 1, f & d$truth == 1), sp1 = rate(d$y == 0, f & d$truth == 0),
    se0 = rate(d$y == 1, !f & d$truth == 1), sp0 = rate(d$y == 0, !f & d$truth == 0))
}))
summary_tab <- data.frame(
  n = nrow(d), aki_kdigo = sum(d$truth), aki_coded = sum(d$y),
  sensitivity = mean(d$y[d$truth == 1]), specificity = mean(1 - d$y[d$truth == 0]),
  auc_vs_kdigo = auc(d$risk, d$truth), auc_vs_coded = auc(d$risk, d$y))

# ---- decision curves --------------------------------------------------------
tc <- nb_curve(d$risk, d$truth, grid); oc <- nb_curve(d$risk, d$y, grid)
bd <- nb_bounds(d$risk, d$y, grid, se0 = c(0.15, 1), sp0 = c(0.90, 1), se1 = c(0.15, 1), sp1 = c(0.90, 1))
# bias analysis with vague beliefs: sensitivity anywhere from low to high (Beta(2, 2)),
# specificity around 0.95 (Beta(95, 5)); draws incompatible with the data are dropped
ba <- nb_bias_analysis(d$risk, d$y, grid, se0 = c(2, 2), sp0 = c(95, 5), draws = 4000)
curves <- data.frame(threshold = grid, nb_true = tc$nb, nb_all_true = tc$nb_all, delta_true = tc$delta_all,
  nb_obs = oc$nb, nb_all_obs = oc$nb_all, delta_obs = oc$delta_all,
  delta_lower = bd$delta_lower, delta_upper = bd$delta_upper,
  delta_lower_ci = bd$delta_lower_ci, delta_upper_ci = bd$delta_upper_ci, se_tip_upper = bd$se_tip_upper,
  ba_delta = ba$delta_all, ba_delta_lo = ba$delta_lo, ba_delta_hi = ba$delta_hi, ba_p_beats_all = ba$p_beats_all,
  nb_lower = bd$nb_lower, nb_upper = bd$nb_upper, se_tip = bd$se_tip, sp_tip = bd$sp_tip,
  se0_actual = sapply(grid, function(p) rate(d$y == 1, d$risk < p & d$truth == 1)),
  sp1_actual = sapply(grid, function(p) rate(d$y == 0, d$risk >= p & d$truth == 0)))

# ---- chart review emulation: true outcome hidden except for reviewed charts --
m <- 150; pilot_n <- 40; R <- 2000; key <- c(0.05, 0.10, 0.15)
st <- make_strata(d$risk, d$y, breaks = c(0.05, 0.10, 0.15, 0.20))
emul <- do.call(rbind, lapply(seq_len(R), function(i) {
  n <- nrow(d); pilot <- sample.int(n, pilot_n)
  v <- pilot_v(d$risk, d$y, ifelse(seq_len(n) %in% pilot, d$truth, NA))
  des <- list(SRS = rep((m - pilot_n) / n, n),
              Optimal = twophase_design(d$risk, d$y, m - pilot_n, c(0.02, 0.30), "both", v = v, floor = 0.01))
  do.call(rbind, lapply(names(des), function(nm) {
    p2 <- des[[nm]]; in_pilot <- seq_len(n) %in% pilot
    prob <- ifelse(in_pilot, 1, p2)
    seen <- in_pilot | rbinom(n, 1, p2) == 1
    f <- suppressWarnings(twophase_dca(d$risk, d$y, ifelse(seen, d$truth, NA), prob, key, strata = st, B = 2))
    data.frame(rep = i, design = nm, threshold = key, nb = f$nb, nb_lo = f$nb_lo, nb_hi = f$nb_hi,
               delta = f$delta_all, delta_lo = f$delta_lo, delta_hi = f$delta_hi, reviewed = sum(seen))
  }))
}))
tk <- nb_curve(d$risk, d$truth, key)
# largest gain the optimal design could give (Theorem 4(iv)), using v from the full data
bins <- interaction(cut(d$risk, unique(quantile(d$risk, 0:10 / 10)), include.lowest = TRUE), d$y)
vt <- ave(d$truth, bins); vt <- vt * (1 - vt)
gain_pred <- design_gain(d$risk, d$y, m - pilot_n, c(0.02, 0.30), "both", v = vt)
emul_sum <- do.call(rbind, lapply(split(emul, list(emul$design, emul$threshold)), function(s) {
  t <- tk[abs(tk$threshold - s$threshold[1]) < 1e-9, ]
  data.frame(design = s$design[1], threshold = s$threshold[1], reviewed = mean(s$reviewed),
    delta_true = t$delta_all, delta_mean = mean(s$delta), delta_empse = sd(s$delta),
    delta_cover = mean(s$delta_lo <= t$delta_all & t$delta_all <= s$delta_hi),
    nb_true = t$nb, nb_mean = mean(s$nb), nb_empse = sd(s$nb),
    nb_cover = mean(s$nb_lo <= t$nb & t$nb <= s$nb_hi))
}))
emul_sum$gain_predicted <- gain_pred[["gain"]]

write.csv(data.frame(step = names(flow), n = unlist(flow)), file.path(out, "eicu_flow.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out, "eicu_summary.csv"), row.names = FALSE)
write.csv(acc, file.path(out, "eicu_region_accuracy.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "eicu_curves.csv"), row.names = FALSE)
write.csv(emul_sum, file.path(out, "eicu_chart_review.csv"), row.names = FALSE)
write.csv(data.frame(patientunitstayid = d$patientunitstayid, risk = d$risk, y = d$y, truth = d$truth),
          "../data/derived/eicu_aki_analysis.csv", row.names = FALSE)
print(flow); print(summary_tab); print(acc); print(emul_sum)
