# Builds the tables of the paper and its supplement (Markdown) from the result files.
# Rows are sorted with method = "radix", so the order is the same in every locale.
# Run from the paper folder: Rscript tables/make_tables.R
md <- function(df, file) {
  out <- c(paste("|", paste(names(df), collapse = " | "), "|"),
           paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
           apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
  writeLines(out, file.path("tables", file))
}
f4 <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 4))
f3 <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 3))
f2 <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 2))
pct <- function(x) paste0(round(100 * x), "%")
mech_lab <- c(M3 = "M3", M4 = "M4")
des_lab <- c(Recorded = "Recorded", `Equal-rate correction` = "Equal-rate", SRS = "SRS", Balanced = "Balanced",
             `Optimal (both)` = "Optimal (both)", `Optimal (treat-all)` = "Optimal (treat-all)")

pb <- read.csv("simulation/results/partB_performance.csv")
re <- function(s, all) {           # variance under SRS / variance under the design, same mech/budget/threshold
  srs <- all[all$design == "SRS" & all$mech == s$mech & all$budget == s$budget & abs(all$threshold - s$threshold) < 1e-9, ]
  srs$delta_empse^2 / s$delta_empse^2
}

# ---- Table 1 -----------------------------------------------------------------
t1 <- pb[pb$budget == 600 & pb$threshold %in% c(0.10, 0.20) & pb$design %in% names(des_lab), ]
t1 <- t1[order(t1$mech, t1$threshold, match(t1$design, names(des_lab))), ]
cr <- grepl("^(SRS|Balanced|Optimal)", t1$design)
md(data.frame(Scenario = mech_lab[t1$mech], Threshold = pct(t1$threshold), Method = des_lab[t1$design],
              Bias = f4(t1$delta_bias), `Empirical SE` = f4(t1$delta_empse), RMSE = f4(t1$rmse_delta),
              Coverage = ifelse(cr, f2(t1$delta_cover), "-"), `Mean width` = ifelse(cr, f3(t1$delta_width), "-"),
              `Relative efficiency` = ifelse(cr, f2(sapply(seq_len(nrow(t1)), function(i) re(t1[i, ], pb))), "-"),
              check.names = FALSE), "table1.md")

# ---- Table S2: full chart review results ------------------------------------
s2 <- pb[grepl("^(SRS|Balanced|Optimal)", pb$design) & !grepl("\\[", pb$design) & pb$threshold %in% c(0.06, 0.10, 0.20, 0.30), ]
s2 <- s2[order(s2$mech, s2$budget, s2$design, s2$threshold, method = "radix"), ]
md(data.frame(Scenario = s2$mech, Charts = s2$budget, Design = s2$design, Threshold = pct(s2$threshold),
              `Bias (gain)` = f4(s2$delta_bias), `SE (gain)` = f4(s2$delta_empse), `Model SE (gain)` = f4(s2$delta_modse),
              `Coverage (gain)` = f2(s2$delta_cover), `Bias (NB)` = f4(s2$nb_bias), `SE (NB)` = f4(s2$nb_empse),
              `Model SE (NB)` = f4(s2$nb_modse), `Coverage (NB)` = f2(s2$nb_cover), check.names = FALSE), "tableS2.md")

# ---- Table S3: bands ----------------------------------------------------------
bc <- read.csv("simulation/results/partB_band_coverage.csv")
bc <- bc[!grepl("\\[", bc$design), ]; bc <- bc[order(bc$mech, bc$budget, bc$design, method = "radix"), ]
md(data.frame(Scenario = bc$mech, Charts = bc$budget, Design = bc$design,
              `Gain over treat-all` = f2(bc$band_cover_delta), `Net benefit` = f2(bc$band_cover_nb), check.names = FALSE),
   "tableS3.md")

# ---- Table S4: interval and working model ------------------------------------
v <- pb[pb$design %in% c("SRS", "SRS [wald]", "SRS [strata]", "Optimal (both)", "Optimal (both) [wald]", "Optimal (both) [strata]") &
        pb$threshold %in% c(0.04, 0.10, 0.20, 0.30), ]
lab <- function(d) ifelse(grepl("wald", d), "Logistic, Wald", ifelse(grepl("strata", d), "Stratum means, split", "Logistic, split (default)"))
v <- v[order(v$mech, v$budget, sub(" \\[.*", "", v$design), lab(v$design), v$threshold, method = "radix"), ]
md(data.frame(Scenario = v$mech, Charts = v$budget, Design = sub(" \\[.*", "", v$design), `Working model, interval` = lab(v$design),
              Threshold = pct(v$threshold), `SE (gain)` = f4(v$delta_empse), `Coverage (gain)` = f2(v$delta_cover),
              `Width (gain)` = f3(v$delta_width), `SE (NB)` = f4(v$nb_empse), `Coverage (NB)` = f2(v$nb_cover),
              `Width (NB)` = f3(v$nb_width), check.names = FALSE), "tableS4.md")

# ---- Table S5: predicted vs realised design gain -----------------------------
gp <- read.csv("simulation/results/design_gain_predicted.csv")
real <- do.call(rbind, lapply(split(pb[pb$design %in% c("SRS", "Optimal (both)", "Optimal (treat-all)"), ],
                                    list(pb$mech[pb$design %in% c("SRS", "Optimal (both)", "Optimal (treat-all)")],
                                         pb$budget[pb$design %in% c("SRS", "Optimal (both)", "Optimal (treat-all)")])), function(s) {
  tot <- function(d, all) sum(s$delta_empse[s$design == d]^2 + if (all) 0 else s$nb_empse[s$design == d]^2)
  data.frame(mech = s$mech[1], budget = s$budget[1],
             both = tot("SRS", FALSE) / tot("Optimal (both)", FALSE), all = tot("SRS", TRUE) / tot("Optimal (treat-all)", TRUE))
}))
g <- merge(reshape(gp[, c("mech", "budget", "comparator", "gain")], idvar = c("mech", "budget"), timevar = "comparator", direction = "wide"), real)
g <- g[order(g$mech, g$budget, method = "radix"), ]
md(data.frame(Scenario = g$mech, Charts = g$budget,
              `Predicted, treat-all design` = f2(g$gain.all), `Realised, treat-all design` = f2(g$all),
              `Predicted, both design` = f2(g$gain.both), `Realised, both design` = f2(g$both), check.names = FALSE), "tableS5.md")

# ---- Table S6: sensitivity ------------------------------------------------------
sc <- read.csv("simulation/results/sens_scenarios.csv"); sd <- read.csv("simulation/results/sens_design.csv")
sc$pilot_model <- "logistic"; ss <- rbind(sc[, names(sd)], sd)
md(data.frame(Scenario = ss$mech, Pilot = ss$pilot, `Pilot model` = ss$pilot_model, `Defensive share` = ss$defensive,
              Design = ifelse(ss$design == "SRS", "Simple random", ss$design), Threshold = pct(ss$threshold),
              `Bias (gain)` = f4(ss$delta_bias), `SE (gain)` = f4(ss$delta_empse), `Coverage (gain)` = f2(ss$delta_cover),
              `Bias (NB)` = f4(ss$nb_bias), `SE (NB)` = f4(ss$nb_empse), `Coverage (NB)` = f2(ss$nb_cover), check.names = FALSE),
   "tableS6.md")

# ---- Table 2: accuracy of the recorded outcome in the examples ---------------
ea <- read.csv("analysis/results/eicu_region_accuracy.csv"); na <- read.csv("analysis/results/nhanes_region_accuracy.csv")
t2 <- rbind(data.frame(Example = "Kidney injury, ICU codes", Threshold = pct(ea$threshold), Flagged = pct(ea$flagged),
                       `Se, flagged` = f2(ea$se1), `Se, unflagged` = f2(ea$se0), `Sp, flagged` = f2(ea$sp1), `Sp, unflagged` = f2(ea$sp0),
                       check.names = FALSE),
            data.frame(Example = ifelse(na$outcome == "diabetes", "Diabetes, self-report", "Kidney disease, self-report"),
                       Threshold = pct(na$threshold), Flagged = pct(na$flagged), `Se, flagged` = f2(na$se1),
                       `Se, unflagged` = f2(na$se0), `Sp, flagged` = "1", `Sp, unflagged` = "1", check.names = FALSE))
md(t2, "table2.md")

# ---- Table 3: what each tool says at a threshold of 10% -----------------------
ec <- read.csv("analysis/results/eicu_curves.csv"); hc <- read.csv("analysis/results/nhanes_curves.csv")
kc <- read.csv("analysis/results/nhanes_ckd_curves.csv")
er <- read.csv("analysis/results/eicu_chart_review.csv"); hr <- read.csv("analysis/results/nhanes_lab_substudy.csv")
at <- function(x) x[abs(x$threshold - 0.10) < 1e-9, ]
e <- at(ec); h <- at(hc); k <- at(kc); e2 <- at(er[er$design == "SRS", ]); h2 <- at(hr[hr$design == "SRS", ])
rng <- function(a, b) paste0(f3(a), " to ", f3(b))
md(data.frame(
  Quantity = c("Comparison", "Recorded-outcome value", "Reference value", "Tipping point (upper or lower 95% limit)",
               "Actual sensitivity", "Bounds, 95% confidence limits", "Bias analysis, median (95% interval)",
               "Review sample: mean estimate, coverage"),
  `Kidney injury (eICU)` = c("Gain over treat-all", f3(e$delta_obs), f3(e$delta_true),
                             paste0("Se0 > ", f2(e$se_tip), " (", f2(e$se_tip_upper), ")"), f2(e$se0_actual),
                             rng(e$delta_lower_ci, e$delta_upper_ci),
                             paste0(f3(e$ba_delta), " (", rng(e$ba_delta_lo, e$ba_delta_hi), ")"),
                             paste0(f3(e2$delta_mean), ", ", f2(e2$delta_cover))),
  `Diabetes (NHANES)` = c("Gain over treat-all", f3(h$delta_obs), f3(h$delta_true),
                          paste0("Se0 > ", f2(h$se_tip), " (", f2(h$se_tip_upper), ")"), f2(h$se0_actual),
                          rng(h$delta_lower_ci, h$delta_upper_ci),
                          paste0(f3(h$ba_delta), " (", rng(h$ba_delta_lo, h$ba_delta_hi), ")"),
                          paste0(f3(h2$delta_mean), ", ", f2(h2$delta_cover))),
  `Kidney disease (NHANES)` = c("Net benefit (vs treat-none)", f3(k$nb_obs), f3(k$nb_true),
                                paste0("Se1 < ", f2(k$se1_tip), " (", f2(k$se1_tip_lower), ")"), f2(k$se1_actual),
                                rng(k$nb_lower_ci, k$nb_upper_ci), "-", "-"),
  check.names = FALSE), "table3.md")

# ---- Tables S7, S8: cohort flows -----------------------------------------------
ef <- read.csv("analysis/results/eicu_flow.csv"); hf <- read.csv("analysis/results/nhanes_flow.csv")
md(data.frame(Step = c("All ICU stays in the demo", "First ICU stay per patient", "Age 18 or older", "ICU stay of at least 24 hours",
                       "At least two creatinine values (from 24 h before to 7 days after admission), one after admission",
                       "No end-stage kidney disease or dialysis (analysis cohort)"),
              Stays = format(ef$n, big.mark = ",")), "tableS7.md")
md(data.frame(Step = c("Participants, NHANES 2011-2018", "Age 20 or older", "Examined, not pregnant",
                       "HbA1c measured and diabetes question answered", "Body mass index and waist measured (diabetes cohort)",
                       "Validation cycles 2015-2018 (diabetes analysis)", 
                       "Kidney disease cohort, validation cycles (steps 2 to 4, with body mass index, kidney question, creatinine and urine albumin; waist not required)"),
              Participants = format(hf$n, big.mark = ",")), "tableS8.md")
cat("tables written\n")
