#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 05_extract_prime_recommendations.R <prime.stdout> <output.cfg> <output.tsv>")
lines <- readLines(args[[1]], warn = FALSE)
marker <- grep("PLACE THE LINES BELOW IN THE CONFIGURATION FILE", lines, fixed = TRUE)
if (length(marker) != 1L) stop("Prime recommendation marker not found exactly once")
tail_lines <- trimws(lines[seq.int(marker + 1L, length(lines))])
rec <- tail_lines[tail_lines %in% c("moredetail", "moredetailad", "moredetailcvad") |
                    grepl("^(opt|optad|optcvad)\\s*=\\s*[0-9]+$", tail_lines)]
keys <- trimws(sub("=.*$", "", rec))
required <- c("opt", "optad", "optcvad")
allowed <- c(required, "moredetail", "moredetailad", "moredetailcvad")
if (!all(required %in% keys) || any(!keys %in% allowed) || anyDuplicated(keys)) {
  stop("Unexpected or incomplete prime recommendations: ", paste(rec, collapse = " | "))
}
writeLines(rec, args[[2]])
tab <- data.frame(
  parameter = keys,
  recommendation = ifelse(grepl("=", rec, fixed = TRUE), trimws(sub("^[^=]*=", "", rec)), "enabled"),
  stringsAsFactors = FALSE
)
write.table(tab, args[[3]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste(rec, collapse = "\n"), "\n")
