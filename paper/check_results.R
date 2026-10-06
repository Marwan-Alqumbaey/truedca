# Compare result files written by a fresh run with the stored ones.
# Usage, from the paper folder: Rscript check_results.R <stored folder> <new folder> [tolerance]
# Numeric columns must agree to the relative tolerance (default 1e-6); other columns must be identical.
args <- commandArgs(trailingOnly = TRUE)
stored <- args[1]; fresh <- args[2]
tol <- if (length(args) > 2) as.numeric(args[3]) else 1e-6
files <- list.files(stored, pattern = "\\.(csv|md)$")
if (!length(files)) stop("No result files found in ", stored)
bad <- character(0)
for (f in files) {
  new <- file.path(fresh, f)
  if (!file.exists(new)) { bad <- c(bad, paste(f, "is missing")); next }
  if (grepl("\\.md$", f)) {
    if (!identical(readLines(file.path(stored, f)), readLines(new))) bad <- c(bad, paste(f, "differs"))
    next
  }
  a <- read.csv(file.path(stored, f), stringsAsFactors = FALSE)
  b <- read.csv(new, stringsAsFactors = FALSE)
  ok <- identical(dim(a), dim(b)) && identical(names(a), names(b)) &&
    all(vapply(names(a), function(v) {
      if (is.numeric(a[[v]]) && is.numeric(b[[v]])) isTRUE(all.equal(a[[v]], b[[v]], tolerance = tol))
      else identical(a[[v]], b[[v]])
    }, logical(1)))
  if (!ok) bad <- c(bad, paste(f, "differs"))
}
if (length(bad)) stop("Results do not match:\n", paste(bad, collapse = "\n"))
cat("All", length(files), "files match the stored results.\n")
