# R packages needed to reproduce the paper. Run from the paper folder: Rscript install.R
install.packages(c("foreign", "digest"), repos = "https://cloud.r-project.org")
install.packages("../r", repos = NULL, type = "source")
