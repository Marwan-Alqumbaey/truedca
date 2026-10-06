# Reference outputs from the R package, used by the Python tests.
# Run from python/: Rscript tests/fixtures/make_fixtures.R
library(truedca)
set.seed(42)
n <- 3000
risk <- round(plogis(rnorm(n, -2, 1.1)), 6)
truth <- rbinom(n, 1, risk)
y <- ifelse(truth == 1, rbinom(n, 1, ifelse(risk > 0.15, 0.85, 0.5)), rbinom(n, 1, 0.01))
prob <- ifelse(y == 1, 0.5, ifelse(risk < 0.1, 0.08, 0.2))
seen <- rbinom(n, 1, prob) == 1
tt <- ifelse(seen, truth, NA)
st <- make_strata(risk, y, breaks = c(0.05, 0.1, 0.2))
dir <- "tests/fixtures"
write.csv(data.frame(risk, y, truth = tt, prob, strata = as.integer(st), strata5 = as.integer(make_strata(risk, y))),
          file.path(dir, "data.csv"), row.names = FALSE, na = "")
th <- c(0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
w <- function(x, f) write.csv(as.data.frame(x), file.path(dir, f), row.names = FALSE, na = "")
w(nb_curve(risk, y, th), "nb_curve.csv")
w(suppressWarnings(nb_adjust(risk, y, th, se0 = 0.6, sp0 = 0.99, se1 = 0.85, sp1 = 0.98)), "nb_adjust.csv")
w(nb_bounds(risk, y, th, se0 = c(0.4, 1), sp0 = c(0.97, 1), se1 = c(0.7, 1), sp1 = c(0.98, 1), tip_sp0 = 0.99, tip_sp1 = 0.995), "nb_bounds.csv")
f <- twophase_dca(risk, y, tt, prob, th, strata = st, B = 10, working = "strata")
w(f[, c("threshold", "nb", "nb_se", "nb_lo", "nb_hi", "delta_all", "delta_se", "delta_lo", "delta_hi")], "twophase_strata.csv")
for (iv in c("wald", "wilson")) {
  fw <- twophase_dca(risk, y, tt, prob, th, strata = st, B = 10, working = "strata", interval = iv)
  w(fw[, c("threshold", "nb_lo", "nb_hi", "delta_lo", "delta_hi")], paste0("twophase_", iv, ".csv"))
}
m <- ave(ifelse(is.na(tt), 0, tt) / prob * seen, st, FUN = sum) / ave(seen / prob, st, FUN = sum)
v <- m * (1 - m)
w(data.frame(v = v, all = twophase_design(risk, y, 300, c(0.02, 0.3), "all", v = v, floor = 0.01),
             both = twophase_design(risk, y, 300, c(0.02, 0.3), "both", v = v, floor = 0.01, defensive = 0)), "design.csv")
w(t(design_gain(risk, y, 300, c(0.02, 0.3), "both", v = v)), "design_gain.csv")
lr <- qlogis(pmin(pmax(risk, 1e-6), 1 - 1e-6))
g <- glm(truth ~ lr * y, family = quasibinomial(), weights = 1 / prob, data = data.frame(truth = tt, lr, y, prob)[seen, ])
w(t(coef(g)), "glm_coef.csv")
pil <- ifelse(seq_len(n) %in% 1:150, truth, NA)
w(data.frame(v = pilot_v(risk, y, pil)), "pilot_v.csv")
w(data.frame(truth = pil), "pilot_truth.csv")
# edge case: only one reviewed record in the pilot, so the interaction is not estimable
pil2 <- rep(NA, n); idx <- c(which(y == 0)[1:60], which(y == 1)[1]); pil2[idx] <- truth[idx]
w(data.frame(truth = pil2, v = pilot_v(risk, y, pil2)), "pilot_v_edge.csv")
