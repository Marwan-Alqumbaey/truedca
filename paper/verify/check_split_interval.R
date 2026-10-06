# Coverage of split score and Wald intervals by review fraction (Supplementary S1.7).
# Run from the paper folder: Rscript verify/check_split_interval.R
suppressPackageStartupMessages(library(truedca))
set.seed(1); P <- c(0.1, 0.2)
big <- 2e6; r <- rbeta(big, 2, 8); t <- rbinom(big, 1, r); tr <- nb_curve(r, t, P)
for (f in c(1, 0.5, 0.2, 0.05)) {
  res <- replicate(300, {
    n <- 2000; r <- rbeta(n, 2, 8); t <- rbinom(n, 1, r)
    y <- ifelse(t == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.01))
    prob <- rep(f, n); tt <- if (f == 1) t else ifelse(rbinom(n, 1, prob) == 1, t, NA)
    out <- c()
    for (iv in c("split", "wald")) {
      x <- suppressWarnings(twophase_dca(r, y, tt, prob, P, B = 2, interval = iv))
      out <- c(out, x$nb_lo <= tr$nb & tr$nb <= x$nb_hi, x$delta_lo <= tr$delta_all & tr$delta_all <= x$delta_hi)
    }
    out
  })
  m <- matrix(rowMeans(res), 2); stopifnot(all(m[, 1:2] >= 0.92)); cat(sprintf("f=%.2f  split nb %.2f %.2f d %.2f %.2f | wald nb %.2f %.2f d %.2f %.2f\n", f, m[1,1], m[2,1], m[1,2], m[2,2], m[1,3], m[2,3], m[1,4], m[2,4]))
}
