# scripts/update_metrics.R
# install.packages(c("scholar","jsonlite"), repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(scholar)
  library(jsonlite)
})
out_dir <- file.path("data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outfile <- file.path(out_dir, "metrics.json")
SCHOLAR_ID <- "C2VNls8AAAAJ" # Paul Garrett

# Fetch + write ---------------------------------------------------------
ok <- TRUE
prof <- try(get_profile(SCHOLAR_ID), silent = TRUE)
if (inherits(prof, "try-error")) ok <- FALSE

if (ok) {
  metrics <- list(
    source      = "google_scholar",
    updated_utc = format(Sys.time(), "%a %Y-%m-%d", tz = ""),
    h_index     = prof$h_index,
    i10_index   = prof$i10_index,
    citations   = prof$total_cites
  )
  write_json(metrics, outfile, auto_unbox = TRUE, pretty = TRUE)
} else {
  message("GS fetch failed; leaving existing data/metrics.json in place (if any).")
}