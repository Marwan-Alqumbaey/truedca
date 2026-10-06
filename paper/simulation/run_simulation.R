# Simulation study (ADEMP; Morris, White & Crowther 2019)
# Aims: (A) size of error in decision curves built on recorded outcomes;
#       (B) performance of two-phase chart review estimators and designs.
# Run from the paper folder, in four steps (the two scenarios can run in parallel):
#   Rscript simulation/run_simulation.R A        # part A
#   Rscript simulation/run_simulation.R M3       # part B, scenario M3
#   Rscript simulation/run_simulation.R M4       # part B, scenario M4
#   Rscript simulation/run_simulation.R summary  # performance tables
suppressPackageStartupMessages(library(truedca))
step <- commandArgs(trailingOnly = TRUE)[1]
seeds <- c(A = 20261004, M3 = 20261005, M4 = 20261006, summary = 1)
set.seed(seeds[[step]])
dir.create("simulation/results", showWarnings = FALSE)

grid <- seq(0.02, 0.30, by = 0.02)
key <- c(0.05, 0.10, 0.20)

# ---- Data-generating mechanism ----------------------------------------------
# True risk depends on an observed marker z and an unobserved factor u.
# The model uses z only and is calibrated: r = P(T = 1 | z).
cal <- local({
  zz <- seq(-5, 5, by = 0.01)
  uu <- stats::qnorm(seq(0.0005, 0.9995, by = 0.001))
  vals <- sapply(zz, function(z) mean(plogis(-2.7 + 1.0 * z + 0.8 * uu)))
  stats::approxfun(zz, vals, rule = 2)
})
gen <- function(n, mech) {
  z <- rnorm(n); u <- rnorm(n)
  truth <- rbinom(n, 1, plogis(-2.7 + z + 0.8 * u))
  r <- cal(z)
  se <- switch(mech,
    M1 = rep(0.85, n),
    M2 = rep(0.60, n),
    M3 = plogis(0.4 + 1.0 * z),             # sicker patients recorded more often
    M4 = ifelse(r >= 0.10, 0.95, 0.50))     # deployed alert at 10% drives testing
  sp <- switch(mech, M1 = 0.98, M2 = 0.995, M3 = 0.99, M4 = 0.995)
  y <- ifelse(truth == 1, rbinom(n, 1, se), rbinom(n, 1, 1 - sp))
  data.frame(r = r, truth = truth, y = y)
}
mechs <- c(M1 = "Non-differential (Se 0.85, Sp 0.98)",
           M2 = "Low sensitivity (Se 0.60, Sp 0.995)",
           M3 = "Recording rises with severity",
           M4 = "Deployed alert drives testing")

# ---- Part A: population-level distortion (n = 1e6) --------------------------
if (step == "A") {
partA <- do.call(rbind, lapply(names(mechs), function(m) {
  d <- gen(1e6, m)
  tr <- nb_curve(d$r, d$truth, grid); ob <- nb_curve(d$r, d$y, grid)
  rnd <- runif(1e6)                       # uninformative score
  rt <- nb_curve(rnd, d$truth, grid); ro <- nb_curve(rnd, d$y, grid)
  data.frame(mech = m, threshold = grid,
             nb_true = tr$nb, nb_obs = ob$nb,
             delta_true = tr$delta_all, delta_obs = ob$delta_all,
             random_delta_true = rt$delta_all, random_delta_obs = ro$delta_all,
             prevalence = mean(d$truth), recorded_prevalence = mean(d$y))
}))
write.csv(partA, "simulation/results/partA_distortion.csv", row.names = FALSE)
}

# ---- Part B: two-phase estimation ------------------------------------------
n <- 20000; reps <- as.integer(Sys.getenv("SIM_REPS", "500")); B <- 200
truth_of <- function(m) {
  set.seed(seeds[[m]] + 1e6)
  d <- gen(2e6, m); tc <- nb_curve(d$r, d$truth, grid)
  # maximum gain of the optimal design (Theorem 4), using the true v(r, y) on a sample of n
  br <- unique(quantile(d$r, seq(0, 1, 0.01))); mt <- tapply(d$truth, list(cut(d$r, br, include.lowest = TRUE), d$y), mean)
  dg <- d[seq_len(n), ]; vt <- mt[cbind(as.integer(cut(dg$r, br, include.lowest = TRUE)), dg$y + 1)]
  vt <- ifelse(is.na(vt), mean(d$truth), vt); vt <- vt * (1 - vt)
  gains <- do.call(rbind, lapply(c(300, 600, 1200), function(b) do.call(rbind, lapply(c("all", "both"), function(cmp)
    data.frame(mech = m, budget = b, comparator = cmp, t(design_gain(dg$r, dg$y, b - 100, c(0.02, 0.30), cmp, v = vt)))))))
  list(curve = data.frame(threshold = grid, nb = tc$nb, delta_all = tc$delta_all), gains = gains)
}

one_rep <- function(m, budget) {
  d <- gen(n, m)
  st <- make_strata(d$r, d$y, breaks = c(0.05, 0.10, 0.20, 0.30))
  # pilot: simple random sample of 100 charts, used only to plan the design
  pilot <- sample.int(n, 100)
  # guess Var(T | r, Y) from the pilot with the logistic working model
  v <- pilot_v(d$r, d$y, ifelse(seq_len(n) %in% pilot, d$truth, NA))
  main <- budget - 100
  designs <- list(
    SRS = rep(main / n, n),
    Balanced = { w <- 1 / as.numeric(table(st)[as.character(st)]); main * w / sum(w) },
    `Optimal (both)` = twophase_design(d$r, d$y, main, c(0.02, 0.30), "both", v = v, floor = 0.002),
    `Optimal (treat-all)` = twophase_design(d$r, d$y, main, c(0.02, 0.30), "all", v = v, floor = 0.002)
  )
  designs$Balanced <- pmin(designs$Balanced, 1)
  res <- list()
  in_pilot <- seq_len(n) %in% pilot
  for (nm in names(designs)) {
    p2 <- designs[[nm]]
    # condition on the pilot: pilot charts are certain, the rest use their wave-2 probability
    prob <- ifelse(in_pilot, 1, p2)
    seen <- in_pilot | rbinom(n, 1, p2) == 1
    obs_t <- ifelse(seen, d$truth, NA)
    # default analysis (logistic working model, split score interval); for SRS and the
    # optimal design also the Wald interval and the stratum-mean working model, on the same sample
    variants <- list(default = list(interval = "split", working = "logistic", B = B))
    if (nm %in% c("SRS", "Optimal (both)"))
      variants <- c(variants, list(wald = list(interval = "wald", working = "logistic", B = 2),
                                   strata = list(interval = "split", working = "strata", B = 2)))
    s0 <- sample.int(1e9, 1)
    for (vn in names(variants)) {
      a <- variants[[vn]]; set.seed(s0)
      f <- suppressWarnings(twophase_dca(d$r, d$y, obs_t, prob, grid, strata = st, B = a$B,
                                         interval = a$interval, working = a$working))
      lab <- if (vn == "default") nm else paste0(nm, " [", vn, "]")
      res[[lab]] <- data.frame(design = lab, threshold = grid,
        nb = f$nb, nb_se = f$nb_se, nb_lo = f$nb_lo, nb_hi = f$nb_hi,
        nb_slo = if (vn == "default") f$nb_slo else NA, nb_shi = if (vn == "default") f$nb_shi else NA,
        delta = f$delta_all, delta_se = f$delta_se, delta_lo = f$delta_lo,
        delta_hi = f$delta_hi, delta_slo = if (vn == "default") f$delta_slo else NA,
        delta_shi = if (vn == "default") f$delta_shi else NA,
        n_reviewed = sum(seen))
    }
    set.seed(s0 + 1)
  }
  # comparators: recorded outcomes, and correction assuming equal error rates
  ob <- nb_curve(d$r, d$y, grid)
  srs <- seq_len(n) %in% pilot | rbinom(n, 1, main / n) == 1
  se_hat <- mean(d$y[srs & d$truth == 1]); sp_hat <- mean(1 - d$y[srs & d$truth == 0])
  adj <- suppressWarnings(nb_adjust(d$r, d$y, grid, se = se_hat, sp = sp_hat))
  na <- rep(NA_real_, length(grid))
  res$Recorded <- data.frame(design = "Recorded", threshold = grid, nb = ob$nb, nb_se = na,
    nb_lo = na, nb_hi = na, nb_slo = na, nb_shi = na, delta = ob$delta_all, delta_se = na,
    delta_lo = na, delta_hi = na, delta_slo = na, delta_shi = na, n_reviewed = 0)
  res$EqualRates <- data.frame(design = "Equal-rate correction", threshold = grid,
    nb = adj$nb, nb_se = na, nb_lo = na, nb_hi = na, nb_slo = na, nb_shi = na,
    delta = adj$delta_all, delta_se = na, delta_lo = na, delta_hi = na,
    delta_slo = na, delta_shi = na, n_reviewed = sum(srs))
  do.call(rbind, res)
}

if (step %in% c("M3", "M4")) {
  truth_m <- truth_of(step); set.seed(seeds[[step]])
  out <- list()
  for (budget in c(300, 600, 1200)) {
    t0 <- Sys.time()
    ck <- file.path("simulation/results", paste0("ckpt_", step, "_", budget, ".rds"))  # checkpoint per budget
    if (file.exists(ck)) { out[[as.character(budget)]] <- readRDS(ck); next }
    set.seed(seeds[[step]] + budget)
    o <- do.call(rbind, lapply(seq_len(reps), function(i) cbind(rep = i, one_rep(step, budget))))
    o$mech <- step; o$budget <- budget; out[[as.character(budget)]] <- o; saveRDS(o, ck)
    message(step, " budget ", budget, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  }
  saveRDS(list(partB = do.call(rbind, out), truth = truth_m),
          file.path("simulation/results", paste0("partB_", step, ".rds")))
}

# ---- Performance measures with Monte Carlo SEs ------------------------------
if (step == "summary") {
raw <- lapply(c(M3 = "M3", M4 = "M4"), function(m) readRDS(file.path("simulation/results", paste0("partB_", m, ".rds"))))
partB <- do.call(rbind, lapply(raw, `[[`, "partB"))
truth_pop <- lapply(raw, function(x) x$truth$curve)
write.csv(do.call(rbind, lapply(raw, function(x) x$truth$gains)), "simulation/results/design_gain_predicted.csv", row.names = FALSE)
perf <- do.call(rbind, lapply(split(partB, list(partB$mech, partB$budget, partB$design, partB$threshold), drop = TRUE), function(s) {
  tr <- truth_pop[[s$mech[1]]]; tr <- tr[abs(tr$threshold - s$threshold[1]) < 1e-9, ]
  k <- nrow(s)
  b_d <- mean(s$delta) - tr$delta_all; esd_d <- sd(s$delta)
  b_n <- mean(s$nb) - tr$nb; esd_n <- sd(s$nb)
  cov_d <- mean(s$delta_lo <= tr$delta_all & tr$delta_all <= s$delta_hi)
  cov_n <- mean(s$nb_lo <= tr$nb & tr$nb <= s$nb_hi)
  data.frame(mech = s$mech[1], budget = s$budget[1], design = s$design[1], threshold = s$threshold[1],
    delta_bias = b_d, delta_bias_mcse = esd_d / sqrt(k), delta_empse = esd_d,
    delta_modse = sqrt(mean(s$delta_se^2)), delta_cover = cov_d,
    nb_bias = b_n, nb_bias_mcse = esd_n / sqrt(k), nb_empse = esd_n,
    nb_modse = sqrt(mean(s$nb_se^2)), nb_cover = cov_n,
    rmse_delta = sqrt(b_d^2 + esd_d^2), delta_width = mean(s$delta_hi - s$delta_lo),
    nb_width = mean(s$nb_hi - s$nb_lo), mean_reviewed = mean(s$n_reviewed))
}))
# simultaneous coverage over the whole grid
simcov <- do.call(rbind, lapply(split(partB, list(partB$mech, partB$budget, partB$design), drop = TRUE), function(s) {
  if (all(is.na(s$delta_slo))) return(NULL)
  tr <- truth_pop[[s$mech[1]]]
  ok <- sapply(split(s, s$rep), function(r) {
    r <- r[order(r$threshold), ]
    c(all(r$delta_slo <= tr$delta_all & tr$delta_all <= r$delta_shi),
      all(r$nb_slo <= tr$nb & tr$nb <= r$nb_shi))
  })
  data.frame(mech = s$mech[1], budget = s$budget[1], design = s$design[1],
             band_cover_delta = mean(ok[1, ]), band_cover_nb = mean(ok[2, ]))
}))
write.csv(perf, "simulation/results/partB_performance.csv", row.names = FALSE)
write.csv(simcov, "simulation/results/partB_band_coverage.csv", row.names = FALSE)
writeLines(c(paste("R", getRversion()), R.version$platform, paste("truedca", packageVersion("truedca"))),
           "simulation/results/software.txt")
}
