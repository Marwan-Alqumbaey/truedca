# Download the NHANES public-use files used in the paper (cycles 2011-12 to 2017-18)
# into this folder and check them against SHA256SUMS.txt.
# Run from the repository root: Rscript data/nhanes/download_nhanes.R
dest <- "data/nhanes"
files <- c("DEMO", "DIQ", "GHB", "BMX", "MCQ", "BPQ", "SMQ", "PAQ", "KIQ_U", "BIOPRO", "ALB_CR")
cycles <- c(G = 2011, H = 2013, I = 2015, J = 2017)
options(timeout = 600)
for (s in names(cycles)) for (f in files) {
  name <- paste0(f, "_", s, ".xpt")
  out <- file.path(dest, name)
  if (file.exists(out)) next
  url <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s", cycles[[s]], name)
  message("Downloading ", name)
  download.file(url, out, mode = "wb", quiet = TRUE)
}
sums <- read.table(file.path(dest, "SHA256SUMS.txt"), col.names = c("sha256", "file"))
if (requireNamespace("digest", quietly = TRUE)) {
  got <- vapply(file.path(dest, sums$file), digest::digest, character(1), algo = "sha256", file = TRUE)
  bad <- sums$file[got != sums$sha256]
  if (length(bad)) warning("Checksum differs from the file used in the paper: ", paste(bad, collapse = ", "))
  else message("All ", nrow(sums), " files match the versions used in the paper.")
} else message("Install the 'digest' package to verify checksums.")
