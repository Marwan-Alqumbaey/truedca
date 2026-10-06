# Supplementary simulations: chart review under M1/M2, and sensitivity of the
# optimal design to pilot size, defensive share and pilot model (600 charts, 300 repetitions).
# Usage: Rscript simulation/run_sensitivity.R scenarios | design
suppressPackageStartupMessages(library(truedca))
step <- commandArgs(trailingOnly = TRUE)[1]
src <- readLines("simulation/run_simulation.R")
eval(parse(text = src[grep("^cal <- local", src):(grep("^mechs <- ", src) - 1)]))
grid <- seq(0.02, 0.30, by = 0.02); key <- c(0.10, 0.20); n <- 20000; reps <- 300

run <- function(m, budget, pilot_n, defensive, designs, pilot_type = "logistic") {
  set.seed(20261009); d0 <- gen(2e6, m); tr <- nb_curve(d0$r, d0$truth, key)
  out <- do.call(rbind, lapply(seq_len(reps), function(i) {
    d <- gen(n, m); st <- make_strata(d$r, d$y, breaks = c(0.05, 0.10, 0.20, 0.30))
    pilot <- sample.int(n, pilot_n); in_pilot <- seq_len(n) %in% pilot
    v <- if (pilot_type == "logistic") pilot_v(d$r, d$y, ifelse(in_pilot, d$truth, NA)) else {
      # stratum means shrunk towards the overall pilot mean (2 pseudo-observations)
      m0 <- mean(d$truth[pilot]); k <- tapply(d$truth[pilot], st[pilot], sum); nh <- tapply(d$truth[pilot], st[pilot], length)
      pm <- (ifelse(is.na(k), 0, k) + 2 * m0) / (ifelse(is.na(nh), 0, nh) + 2)
      as.numeric(pm[as.character(st)] * (1 - pm[as.character(st)]))
    }
    main <- budget - pilot_n
    des <- list(SRS = rep(main / n, n),
                Optimal = twophase_design(d$r, d$y, main, c(0.02, 0.30), "both", v = v,
                                          floor = 0.002, defensive = defensive))[designs]
    do.call(rbind, lapply(names(des), function(nm) {
      p2 <- des[[nm]]; prob <- ifelse(in_pilot, 1, p2); seen <- in_pilot | rbinom(n, 1, p2) == 1
      f <- suppressWarnings(twophase_dca(d$r, d$y, ifelse(seen, d$truth, NA), prob, key, strata = st, B = 2))
      data.frame(design = nm, threshold = key, delta = f$delta_all, dlo = f$delta_lo, dhi = f$delta_hi,
                 nb = f$nb, nlo = f$nb_lo, nhi = f$nb_hi)
    }))
  }))
  do.call(rbind, lapply(split(out, list(out$design, out$threshold)), function(s) {
    t <- tr[abs(tr$threshold - s$threshold[1]) < 1e-9, ]
    data.frame(mech = m, pilot = pilot_n, pilot_model = pilot_type, defensive = defensive, design = s$design[1], threshold = s$threshold[1],
      delta_bias = mean(s$delta) - t$delta_all, delta_empse = sd(s$delta),
      delta_cover = mean(s$dlo <= t$delta_all & t$delta_all <= s$dhi),
      nb_bias = mean(s$nb) - t$nb, nb_empse = sd(s$nb), nb_cover = mean(s$nlo <= t$nb & t$nb <= s$nhi))
  }))
}

if (step == "scenarios") {
  res <- rbind(run("M1", 600, 100, 0.3, c("SRS", "Optimal")), run("M2", 600, 100, 0.3, c("SRS", "Optimal")))
  write.csv(res, "simulation/results/sens_scenarios.csv", row.names = FALSE)
} else {
  res <- do.call(rbind, c(
    lapply(c(50, 200), function(pn) run("M4", 600, pn, 0.3, "Optimal")),
    lapply(c(0, 0.6), function(df) run("M4", 600, 100, df, "Optimal")),
    list(run("M4", 600, 100, 0.3, "Optimal", pilot_type = "strata"))))
  write.csv(res, "simulation/results/sens_design.csv", row.names = FALSE)
}
