# Figures 1-4 (base R graphics; vector PDF and 800 dpi TIFF). Run from the paper folder: Rscript analysis/figures.R
dir.create("figures", showWarnings = FALSE)
ink <- "#1f1f1f"; blue <- "#2a78d6"; orange <- "#eb6834"; grey <- "#8c8c8c"
gain <- function(nb, all) nb - pmax(all, 0)          # gain over the better default
# blank out the creation and modification stamps that pdf() writes (same length, so the file stays valid)
strip_stamps <- function(file) {
  b <- readBin(file, "raw", file.info(file)$size)
  txt <- rawToChar(replace(b, b == as.raw(0), as.raw(32)))
  for (key in c("/CreationDate", "/ModDate")) {
    m <- regexpr(paste0(key, " \\(D:[^)]*\\)"), txt, useBytes = TRUE)
    if (m > 0) b[m:(m + attr(m, "match.length") - 1)] <- charToRaw(" ")
  }
  writeBin(b, file)
}
save_fig <- function(name, w, h, expr) {
  pdf(file.path("figures", paste0(name, ".pdf")), w, h); expr(); dev.off()
  strip_stamps(file.path("figures", paste0(name, ".pdf")))
  tiff(file.path("figures", paste0(name, ".tiff")), w, h, units = "in", res = 800, compression = "lzw"); expr(); dev.off()
}
par_std <- function() par(mar = c(4, 4.2, 2.2, 0.8), mgp = c(2.4, 0.7, 0), las = 1, cex.axis = 0.85,
                          tcl = -0.3, bty = "l", col.axis = ink, col.lab = ink, fg = grey)

# ---- Figure 1: distortion of decision curves (simulation part A) -------------
a <- read.csv("simulation/results/partA_distortion.csv")
a$nb_all_true <- a$nb_true - a$delta_true; a$nb_all_obs <- a$nb_obs - a$delta_obs
titles <- c(M1 = "A. Errors unrelated to risk (Se 0.85, Sp 0.98)", M2 = "B. Low sensitivity (Se 0.60, Sp 0.995)",
            M3 = "C. Recording rises with severity", M4 = "D. Deployed alert drives testing")
save_fig("figure1", 7.5, 6.2, function() {
  par(mfrow = c(2, 2)); par_std()
  for (m in names(titles)) {
    s <- a[a$mech == m, ]
    gt <- gain(s$nb_true, s$nb_all_true); go <- gain(s$nb_obs, s$nb_all_obs)
    plot(s$threshold, gt, type = "n", ylim = range(c(gt, go, -0.01, 0.025)), xlab = "Threshold probability",
         ylab = "Gain over better default", main = titles[m], cex.main = 0.85, font.main = 1, adj = 0)
    abline(h = 0, col = grey, lty = 3)
    lines(s$threshold, gt, col = ink, lwd = 2)
    lines(s$threshold, go, col = blue, lwd = 2)
    if (m == "M1") legend("topright", c("True outcome", "Recorded outcome"), col = c(ink, blue), lwd = 2,
                          bty = "n", cex = 0.8, text.col = ink)
  }
})

# ---- Figure 3: eICU example ---------------------------------------------------
e <- read.csv("analysis/results/eicu_curves.csv")
save_fig("figure3", 7.5, 3.6, function() {
  par(mfrow = c(1, 2)); par_std()
  rng <- range(c(e$nb_true, e$nb_all_true, e$nb_obs, e$nb_all_obs, 0), na.rm = TRUE); rng[1] <- max(rng[1], -0.05)
  plot(e$threshold, e$nb_true, type = "n", ylim = rng, xlab = "Threshold probability", ylab = "Net benefit",
       main = "A. Decision curves, KDIGO vs coded AKI", cex.main = 0.85, font.main = 1, adj = 0)
  abline(h = 0, col = grey, lty = 3)
  lines(e$threshold, e$nb_true, col = ink, lwd = 2); lines(e$threshold, e$nb_all_true, col = ink, lty = 2)
  lines(e$threshold, e$nb_obs, col = blue, lwd = 2); lines(e$threshold, e$nb_all_obs, col = blue, lty = 2)
  legend("topright", c("Model, KDIGO", "Treat all, KDIGO", "Model, coded", "Treat all, coded"),
         col = c(ink, ink, blue, blue), lty = c(1, 2, 1, 2), lwd = c(2, 1, 2, 1), bty = "n", cex = 0.72, text.col = ink)
  ok <- e$threshold >= 0.05
  plot(e$threshold[ok], e$se_tip[ok], type = "l", col = blue, lwd = 2, ylim = c(0, 1.05),
       xlab = "Threshold probability", ylab = "Code sensitivity, unflagged patients",
       main = "B. Tipping point vs actual sensitivity", cex.main = 0.85, font.main = 1, adj = 0)
  lines(e$threshold[ok], pmin(e$se_tip_upper[ok], 1.05), col = blue, lty = 2)
  points(e$threshold[ok], e$se0_actual[ok], pch = 16, cex = 0.6, col = ink)
  legend("topright", c("Tipping point, no false codes", "Upper 95% confidence limit", "Actual (vs KDIGO)"),
         col = c(blue, blue, ink), lty = c(1, 2, NA), pch = c(NA, NA, 16), lwd = c(2, 1, 2), bty = "n", cex = 0.72, text.col = ink)
})

# ---- Figure 2: coverage of 95% intervals by review budget (simulation part B) --
pb <- read.csv("simulation/results/partB_performance.csv")
save_fig("figure2", 7.5, 3.6, function() {
  par(mfrow = c(1, 2)); par_std()
  shades <- c(`300` = "#9ec5f4", `600` = "#3987e5", `1200` = "#104281")
  for (est in c("delta", "nb")) {
    plot(NA, xlim = c(0.02, 0.30), ylim = c(0.5, 1), xlab = "Threshold probability",
         ylab = "Coverage of 95% intervals",
         main = if (est == "delta") "A. Gain over treat all" else "B. Net benefit",
         cex.main = 0.85, font.main = 1, adj = 0)
    abline(h = 0.95, col = grey, lty = 2)
    for (b in names(shades)) for (v in c("SRS", "SRS [wald]")) {
      s <- pb[pb$mech == "M4" & pb$design == v & pb$budget == as.numeric(b), ]
      s <- s[order(s$threshold), ]
      lines(s$threshold, s[[paste0(est, "_cover")]], col = shades[[b]], lwd = 2, lty = if (v == "SRS") 1 else 3)
    }
    if (est == "delta") legend("bottomright", c(paste(names(shades), "charts"), "Split score interval", "Wald interval"),
                               col = c(shades, ink, ink), lwd = 2, lty = c(1, 1, 1, 1, 3), bty = "n", cex = 0.72, text.col = ink)
  }
})

# ---- Figure 4: NHANES example ------------------------------------------------
h <- read.csv("analysis/results/nhanes_curves.csv"); hk <- read.csv("analysis/results/nhanes_ckd_curves.csv")
curve_panel <- function(x, title, lab) {
  rng <- range(c(x$nb_true, x$nb_obs, 0)); rng <- c(max(min(rng[1], -0.03), -0.08), rng[2] + 0.01)
  plot(x$threshold, x$nb_true, type = "n", ylim = rng, xlab = "Threshold probability", ylab = "Net benefit",
       main = title, cex.main = 0.85, font.main = 1, adj = 0)
  abline(h = 0, col = grey, lty = 3)
  lines(x$threshold, x$nb_true, col = ink, lwd = 2); lines(x$threshold, x$nb_all_true, col = ink, lty = 2)
  lines(x$threshold, x$nb_obs, col = blue, lwd = 2); lines(x$threshold, x$nb_all_obs, col = blue, lty = 2)
  legend("topright", c("Model, reference", "Treat all, reference", paste("Model,", lab), paste("Treat all,", lab)),
         col = c(ink, ink, blue, blue), lty = c(1, 2, 1, 2), lwd = c(2, 1, 2, 1), bty = "n", cex = 0.8, text.col = ink)
}
save_fig("figure4", 7.5, 3.0, function() {
  par(mfrow = c(1, 3)); par_std(); par(mar = c(4, 4.2, 2.2, 0.5))
  curve_panel(h, "A. Diabetes", "self-report")
  curve_panel(hk, "B. Chronic kidney disease", "self-report")
  ok <- hk$threshold >= 0.05
  plot(hk$threshold[ok], pmin(hk$se1_tip[ok], 1), type = "n", ylim = c(0, 1), xlab = "Threshold probability",
       ylab = "Sensitivity of self-report, flagged", main = "C. Kidney disease: tipping point",
       cex.main = 0.85, font.main = 1, adj = 0)
  polygon(c(hk$threshold[ok], rev(hk$threshold[ok])), pmin(c(hk$se1_tip_lower[ok], rev(hk$se1_tip_upper[ok])), 1),
          col = adjustcolor(blue, 0.2), border = NA)
  lines(hk$threshold[ok], pmin(hk$se1_tip[ok], 1), col = blue, lwd = 2)
  points(hk$threshold[ok], hk$se1_actual[ok], pch = 16, cex = 0.6, col = ink)
  legend("topright", c("Tipping point (95% CI)", "Actual"), col = c(blue, ink), lty = c(1, NA), pch = c(NA, 16),
         lwd = c(2, NA), bty = "n", cex = 0.8, text.col = ink)
})
