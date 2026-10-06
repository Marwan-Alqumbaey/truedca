# Example 2: diabetes (and chronic kidney disease) in NHANES public-use files.
# Recorded outcome: self-reported doctor diagnosis. Reference outcome: the diagnosis,
# or treatment, or a laboratory value (HbA1c >= 6.5%; for CKD eGFR < 60 or ACR >= 30).
# Every recorded case is a reference case, so the recorded outcome has specificity 1
# and only misses undiagnosed disease (the setting of Theorem 2).
# Usage: Rscript analysis/nhanes.R <folder with NHANES .xpt files>
suppressPackageStartupMessages({ library(truedca); library(foreign) })
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "../data/nhanes"
out <- "analysis/results"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(2027)

files <- c("DEMO", "DIQ", "GHB", "BMX", "MCQ", "BPQ", "SMQ", "PAQ", "KIQ_U", "BIOPRO", "ALB_CR")
cycles <- c(G = "2011-12", H = "2013-14", I = "2015-16", J = "2017-18")
have <- names(cycles)[vapply(names(cycles), function(s)
  all(file.exists(file.path(dir, paste0(files, "_", s, ".xpt")))), logical(1))]
if (!length(have)) stop("No complete NHANES cycle found in ", dir)
message("Cycles used: ", paste(cycles[have], collapse = ", "))

read_cycle <- function(s) {
  rd <- function(f) read.xport(file.path(dir, paste0(f, "_", s, ".xpt")))
  keep <- list(DEMO = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "RIDEXPRG", "RIDSTATR"),
               DIQ = c("DIQ010", "DIQ050", "DIQ070"), GHB = "LBXGH", BMX = c("BMXBMI", "BMXWAIST"),
               MCQ = "MCQ300C", BPQ = "BPQ020", SMQ = "SMQ020", PAQ = c("PAQ650", "PAQ665"),
               KIQ_U = c("KIQ022", "KIQ025"), BIOPRO = "LBXSCR", ALB_CR = "URDACT")
  d <- NULL
  for (f in names(keep)) {
    x <- rd(f); x <- x[, c("SEQN", intersect(keep[[f]], names(x))), drop = FALSE]
    d <- if (is.null(d)) x else merge(d, x, by = "SEQN", all.x = TRUE)
  }
  d$cycle <- s
  d
}
raw <- do.call(rbind, lapply(have, read_cycle))

# ---- cohort -----------------------------------------------------------------
flow <- list(all = nrow(raw))
d <- raw[raw$RIDAGEYR >= 20, ]; flow$adults <- nrow(d)
d <- d[d$RIDSTATR == 2 & (is.na(d$RIDEXPRG) | d$RIDEXPRG != 1), ]; flow$examined_not_pregnant <- nrow(d)
d <- d[!is.na(d$LBXGH) & d$DIQ010 %in% 1:3, ]; flow$hba1c_and_answer <- nrow(d)
d <- d[!is.na(d$BMXBMI) & !is.na(d$BMXWAIST), ]; flow$body_measures <- nrow(d)

yes <- function(x) as.integer(!is.na(x) & x == 1)
d$y <- yes(d$DIQ010)                                            # recorded: told by a doctor
d$truth <- as.integer(d$y == 1 | yes(d$DIQ050) == 1 | yes(d$DIQ070) == 1 | d$LBXGH >= 6.5)
d$age <- d$RIDAGEYR; d$female <- as.integer(d$RIAGENDR == 2)
d$race <- factor(d$RIDRETH3, levels = c(3, 1, 2, 4, 6, 7),
                 labels = c("White", "MexAm", "OtherHisp", "Black", "Asian", "Other"))
d$famhx <- factor(ifelse(d$MCQ300C %in% 1, "yes", ifelse(d$MCQ300C %in% 2, "no", "unknown")))
d$htn <- yes(d$BPQ020); d$smoke <- yes(d$SMQ020); d$active <- as.integer(yes(d$PAQ650) | yes(d$PAQ665))
d$pir_miss <- as.integer(is.na(d$INDFMPIR)); d$pir <- ifelse(is.na(d$INDFMPIR), median(d$INDFMPIR, na.rm = TRUE), d$INDFMPIR)

# ---- risk model: no laboratory values, fitted to the reference outcome ------
# Development on the earlier cycles, validation on the later ones; with a single
# cycle, 10-fold cross-fitting instead.
fml <- truth ~ poly(age, 2) + female + race + BMXBMI + BMXWAIST + famhx + htn + smoke + active + pir + pir_miss
if (length(have) >= 2) {
  dev_c <- have[seq_len(floor(length(have) / 2))]
  fit <- glm(fml, binomial, data = d[d$cycle %in% dev_c, ])
  d <- d[!d$cycle %in% dev_c, ]
  d$risk <- predict(fit, d, type = "response")
  split_note <- paste("development", paste(cycles[dev_c], collapse = ", "), "; validation",
                      paste(cycles[setdiff(have, dev_c)], collapse = ", "))
} else {
  fold <- sample(rep(1:10, length.out = nrow(d))); d$risk <- NA_real_
  for (k in 1:10) d$risk[fold == k] <- predict(glm(fml, binomial, data = d[fold != k, ]), d[fold == k, ], type = "response")
  split_note <- "10-fold cross-fitting"
}
flow$validation <- nrow(d)
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }
rate <- function(a, b) if (sum(b)) mean(a[b]) else NA_real_
ci_prop <- function(x, n) { w <- truedca:::wilson(x / n, sqrt(x / n * (1 - x / n) / n), n, 1.96); sprintf("%.2f (%.2f-%.2f)", x / n, w[1], w[2]) }

grid <- seq(0.02, 0.40, by = 0.01)
summary_tab <- data.frame(outcome = "diabetes", split = split_note, n = nrow(d),
  reference_cases = sum(d$truth), recorded_cases = sum(d$y),
  sensitivity = ci_prop(sum(d$y[d$truth == 1]), sum(d$truth)),
  auc_reference = auc(d$risk, d$truth), auc_recorded = auc(d$risk, d$y),
  mean_risk = mean(d$risk))
acc <- do.call(rbind, lapply(c(0.05, 0.10, 0.15, 0.20), function(p) {
  f <- d$risk >= p
  data.frame(outcome = "diabetes", threshold = p, flagged = mean(f),
             se1 = rate(d$y == 1, f & d$truth == 1), se0 = rate(d$y == 1, !f & d$truth == 1),
             n_cases0 = sum(!f & d$truth == 1))
}))

# ---- decision curves, bounds with confidence limits, bias analysis ----------
tc <- nb_curve(d$risk, d$truth, grid); oc <- nb_curve(d$risk, d$y, grid)
bd <- nb_bounds(d$risk, d$y, grid, se0 = c(0.6, 1), sp0 = c(1, 1), se1 = c(0.6, 1), sp1 = c(1, 1))
# external prior: about a quarter of adults with diabetes are undiagnosed -> Se ~ Beta(75, 25);
# self-report is rarely false -> Sp ~ Beta(995, 5)
ba <- nb_bias_analysis(d$risk, d$y, grid, se0 = c(75, 25), sp0 = c(995, 5), draws = 4000)
curves <- data.frame(outcome = "diabetes", threshold = grid,
  nb_true = tc$nb, nb_all_true = tc$nb_all, delta_true = tc$delta_all,
  nb_obs = oc$nb, nb_all_obs = oc$nb_all, delta_obs = oc$delta_all,
  delta_lower = bd$delta_lower, delta_upper = bd$delta_upper,
  delta_lower_ci = bd$delta_lower_ci, delta_upper_ci = bd$delta_upper_ci,
  nb_lower = bd$nb_lower, nb_upper = bd$nb_upper, se_tip = bd$se_tip, se_tip_upper = bd$se_tip_upper,
  se1_tip = bd$se1_tip, se1_tip_lower = bd$se1_tip_lower,
  se0_actual = sapply(grid, function(p) rate(d$y == 1, d$risk < p & d$truth == 1)),
  se1_actual = sapply(grid, function(p) rate(d$y == 1, d$risk >= p & d$truth == 1)),
  ba_delta = ba$delta_all, ba_delta_lo = ba$delta_lo, ba_delta_hi = ba$delta_hi, ba_p_beats_all = ba$p_beats_all,
  ba_nb = ba$nb, ba_nb_lo = ba$nb_lo, ba_nb_hi = ba$nb_hi)

# ---- emulated laboratory sub-study: HbA1c measured in a subsample -----------
# Everyone has the recorded outcome; the reference outcome is revealed only for
# sampled participants. Both designs spend the same expected number of tests.
m <- 400; pilot_n <- 80; R <- 1000; key <- c(0.05, 0.10, 0.20, 0.30)
st <- make_strata(d$risk, d$y, breaks = c(0.05, 0.10, 0.20, 0.30))
n <- nrow(d)
mt <- tapply(d$truth, interaction(cut(d$risk, unique(quantile(d$risk, 0:20 / 20)), include.lowest = TRUE), d$y), mean)
vt <- mt[as.character(interaction(cut(d$risk, unique(quantile(d$risk, 0:20 / 20)), include.lowest = TRUE), d$y))]
vt <- ifelse(is.na(vt), mean(d$truth), vt); vt <- vt * (1 - vt)
gain_pred <- design_gain(d$risk, d$y, m - pilot_n, c(0.02, 0.40), "both", v = vt)
emul <- do.call(rbind, lapply(seq_len(R), function(i) {
  pilot <- sample.int(n, pilot_n); in_pilot <- seq_len(n) %in% pilot
  v <- pilot_v(d$risk, d$y, ifelse(seq_len(n) %in% pilot, d$truth, NA))
  des <- list(SRS = rep((m - pilot_n) / n, n),
              Optimal = twophase_design(d$risk, d$y, m - pilot_n, c(0.02, 0.40), "both", v = v, floor = 0.01),
              # same design with v from the full data: shows how much is lost by guessing v from the pilot
              `Optimal, known v` = twophase_design(d$risk, d$y, m - pilot_n, c(0.02, 0.40), "both", v = vt, floor = 0.01))
  do.call(rbind, lapply(names(des), function(nm) {
    p2 <- des[[nm]]; prob <- ifelse(in_pilot, 1, p2); seen <- in_pilot | rbinom(n, 1, p2) == 1
    f <- suppressWarnings(twophase_dca(d$risk, d$y, ifelse(seen, d$truth, NA), prob, key, strata = st, B = 2))
    data.frame(rep = i, design = nm, threshold = key, nb = f$nb, nb_lo = f$nb_lo, nb_hi = f$nb_hi,
               delta = f$delta_all, delta_lo = f$delta_lo, delta_hi = f$delta_hi, tested = sum(seen))
  }))
}))
tk <- nb_curve(d$risk, d$truth, key)
emul_sum <- do.call(rbind, lapply(split(emul, list(emul$design, emul$threshold)), function(s) {
  t <- tk[abs(tk$threshold - s$threshold[1]) < 1e-9, ]
  data.frame(design = s$design[1], threshold = s$threshold[1], tested = mean(s$tested),
    delta_true = t$delta_all, delta_mean = mean(s$delta), delta_empse = sd(s$delta),
    delta_cover = mean(s$delta_lo <= t$delta_all & t$delta_all <= s$delta_hi),
    nb_true = t$nb, nb_mean = mean(s$nb), nb_empse = sd(s$nb),
    nb_cover = mean(s$nb_lo <= t$nb & t$nb <= s$nb_hi))
}))
emul_sum$gain_predicted <- gain_pred[["gain"]]

# ---- secondary outcome: chronic kidney disease ------------------------------
k <- raw[raw$RIDAGEYR >= 20 & raw$RIDSTATR == 2 & (is.na(raw$RIDEXPRG) | raw$RIDEXPRG != 1) &
         raw$KIQ022 %in% 1:2 & !is.na(raw$LBXSCR) & !is.na(raw$URDACT) & !is.na(raw$BMXBMI) &
         !is.na(raw$LBXGH) & raw$DIQ010 %in% 1:3, ]
fem <- k$RIAGENDR == 2; kap <- ifelse(fem, 0.7, 0.9); alp <- ifelse(fem, -0.241, -0.302)
k$egfr <- 142 * pmin(k$LBXSCR / kap, 1)^alp * pmax(k$LBXSCR / kap, 1)^-1.200 * 0.9938^k$RIDAGEYR * ifelse(fem, 1.012, 1)
k$y <- yes(k$KIQ022)
k$truth <- as.integer(k$y == 1 | yes(k$KIQ025) == 1 | k$egfr < 60 | k$URDACT >= 30)
k$age <- k$RIDAGEYR; k$female <- as.integer(fem)
k$race <- factor(k$RIDRETH3, levels = c(3, 1, 2, 4, 6, 7))
k$dm <- as.integer(yes(k$DIQ010) == 1 | k$LBXGH >= 6.5); k$htn <- yes(k$BPQ020); k$smoke <- yes(k$SMQ020)
kf <- truth ~ poly(age, 2) + female + race + BMXBMI + dm + htn + smoke
if (length(have) >= 2) {
  kfit <- glm(kf, binomial, data = k[k$cycle %in% dev_c, ]); k <- k[!k$cycle %in% dev_c, ]
  k$risk <- predict(kfit, k, type = "response")
} else {
  fold <- sample(rep(1:10, length.out = nrow(k))); k$risk <- NA_real_
  for (j in 1:10) k$risk[fold == j] <- predict(glm(kf, binomial, data = k[fold != j, ]), k[fold == j, ], type = "response")
}
ktc <- nb_curve(k$risk, k$truth, grid); koc <- nb_curve(k$risk, k$y, grid)
kbd <- nb_bounds(k$risk, k$y, grid, se0 = c(0.05, 1), sp0 = c(1, 1), se1 = c(0.05, 1), sp1 = c(1, 1))
summary_tab <- rbind(summary_tab, data.frame(outcome = "ckd", split = split_note, n = nrow(k),
  reference_cases = sum(k$truth), recorded_cases = sum(k$y),
  sensitivity = ci_prop(sum(k$y[k$truth == 1]), sum(k$truth)),
  auc_reference = auc(k$risk, k$truth), auc_recorded = auc(k$risk, k$y), mean_risk = mean(k$risk)))
acc <- rbind(acc, do.call(rbind, lapply(c(0.05, 0.10, 0.15, 0.20), function(p) {
  f <- k$risk >= p
  data.frame(outcome = "ckd", threshold = p, flagged = mean(f),
             se1 = rate(k$y == 1, f & k$truth == 1), se0 = rate(k$y == 1, !f & k$truth == 1),
             n_cases0 = sum(!f & k$truth == 1))
})))
ckd_curves <- data.frame(outcome = "ckd", threshold = grid, nb_true = ktc$nb, nb_all_true = ktc$nb_all,
  delta_true = ktc$delta_all, nb_obs = koc$nb, nb_all_obs = koc$nb_all, delta_obs = koc$delta_all,
  nb_lower = kbd$nb_lower, nb_upper = kbd$nb_upper, nb_lower_ci = kbd$nb_lower_ci, nb_upper_ci = kbd$nb_upper_ci,
  se_tip = kbd$se_tip, se_tip_upper = kbd$se_tip_upper,
  se1_tip = kbd$se1_tip, se1_tip_lower = kbd$se1_tip_lower, se1_tip_upper = kbd$se1_tip_upper,
  se0_actual = sapply(grid, function(p) rate(k$y == 1, k$risk < p & k$truth == 1)),
  se1_actual = sapply(grid, function(p) rate(k$y == 1, k$risk >= p & k$truth == 1)))
flow$ckd_analysis <- nrow(k)

write.csv(data.frame(step = names(flow), n = unlist(flow)), file.path(out, "nhanes_flow.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out, "nhanes_summary.csv"), row.names = FALSE)
write.csv(acc, file.path(out, "nhanes_region_accuracy.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "nhanes_curves.csv"), row.names = FALSE)
write.csv(ckd_curves, file.path(out, "nhanes_ckd_curves.csv"), row.names = FALSE)
write.csv(emul_sum, file.path(out, "nhanes_lab_substudy.csv"), row.names = FALSE)
# derived analysis file (no identifiers beyond the public SEQN) for the repository
write.csv(d[, c("SEQN", "cycle", "risk", "y", "truth")], "../data/derived/nhanes_diabetes_analysis.csv", row.names = FALSE)
print(flow); print(summary_tab); print(acc); print(emul_sum)
print(curves[curves$threshold %in% c(0.05, 0.1, 0.2, 0.3), c("threshold", "delta_true", "delta_obs", "delta_lower_ci", "delta_upper_ci", "se_tip", "se_tip_upper", "se0_actual", "ba_delta", "ba_delta_lo", "ba_delta_hi", "ba_p_beats_all")])
print(ckd_curves[ckd_curves$threshold %in% c(0.05, 0.1, 0.2, 0.3), ])
