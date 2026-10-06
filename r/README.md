# truedca

Decision curve analysis when the outcome is recorded with error.

Decision curves show whether a prediction model does more good than harm at a
given risk threshold. They are usually drawn with outcomes taken from diagnosis
codes, chart abstraction or text mining, and these outcomes are not perfect.
`truedca` shows how such errors distort the curve and how to correct it.

## Install

Requires R 4.1 or later.

```r
# install.packages("remotes")
remotes::install_github("Marwan-Alqumbaey/truedca", subdir = "r")
```

## What it does

| Function | Purpose |
| --- | --- |
| `nb_curve()` | Decision curve from the outcomes as recorded |
| `dca_window()` | Thresholds where a recorded-outcome curve cannot show benefit |
| `warp_threshold()` | Threshold at which recorded outcomes actually judge a decision |
| `nb_adjust()` | Corrected curve when error rates are known (may differ by model flag) |
| `nb_bounds()` | Sharp bounds, tipping points and their confidence limits when error rates are only partly known |
| `nb_bias_analysis()` | Probabilistic bias analysis when no chart review is possible |
| `pilot_v()` | Guess outcome uncertainty from a small pilot review |
| `twophase_design()` | Chart review probabilities that make the corrected curve most precise |
| `design_gain()` | How much that design can gain over simple random sampling |
| `twophase_dca()` | Corrected curve from a chart review subsample, with split score intervals and simultaneous bands |
| `make_strata()` | Strata that cross risk groups with the recorded outcome |
| `plot_dca()` | Simple plot |

## Example

```r
library(truedca)
set.seed(1)
n <- 20000
risk  <- plogis(rnorm(n, -2.5, 1.2))
truth <- rbinom(n, 1, risk)
y     <- ifelse(truth == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.01))  # coded outcome

dca_window(se = 0.7, sp = 0.99)        # floor 0.01, ceiling 0.70

# How robust is the gain over treat all if sensitivity is only known to be >= 0.6?
nb_bounds(risk, y, c(0.05, 0.1, 0.2), se0 = c(0.6, 1), sp0 = c(0.99, 1))

# Review 400 charts: a pilot of 100, then 300 placed where they matter most
pilot <- ifelse(seq_len(n) %in% sample.int(n, 100), truth, NA)
v <- pilot_v(risk, y, pilot)
design_gain(risk, y, 300, range = c(0.02, 0.3), comparator = "both", v = v)
prob <- twophase_design(risk, y, 300, range = c(0.02, 0.3), comparator = "both", v = v)
prob[!is.na(pilot)] <- 1                                         # pilot charts are already reviewed
seen <- ifelse(!is.na(pilot) | rbinom(n, 1, prob) == 1, truth, NA)
fit  <- twophase_dca(risk, y, seen, prob, thresholds = seq(0.02, 0.3, 0.02))
plot_dca(fit)
```

## Citation

If you use `truedca`, please cite the accompanying paper (see `citation("truedca")`). This package is part of the repository <https://github.com/Marwan-Alqumbaey/truedca>, which also holds the Python version and all code for the paper.

## License

MIT
