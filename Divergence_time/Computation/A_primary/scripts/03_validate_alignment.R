#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 03_validate_alignment.R <tree> <fasta> <output.tsv>")
}

suppressPackageStartupMessages(library(ape))
tree_path <- normalizePath(args[[1]], mustWork = TRUE)
fasta_path <- normalizePath(args[[2]], mustWork = TRUE)
out_path <- args[[3]]

tr <- read.tree(tree_path)
con <- file(fasta_path, open = "rt")
on.exit(close(con), add = TRUE)

ids <- character()
lengths <- integer()
current <- NA_integer_
repeat {
  lines <- readLines(con, n = 10000L, warn = FALSE)
  if (!length(lines)) break
  for (line in lines) {
    if (startsWith(line, ">")) {
      id <- sub("^>(\\S+).*$", "\\1", line)
      ids <- c(ids, id)
      lengths <- c(lengths, 0L)
      current <- length(ids)
    } else if (nzchar(line)) {
      if (is.na(current)) stop("Sequence text encountered before first FASTA header")
      lengths[[current]] <- lengths[[current]] + nchar(gsub("[[:space:]]", "", line))
    }
  }
}

checks <- data.frame(
  check = c(
    "fasta_records_495", "fasta_unique_ids", "all_sequences_182682",
    "tree_tip_count_495", "tree_unique_tips", "tree_fasta_tip_sets_identical"
  ),
  pass = c(
    length(ids) == 495L,
    anyDuplicated(ids) == 0L,
    length(lengths) == 495L && all(lengths == 182682L),
    Ntip(tr) == 495L,
    anyDuplicated(tr$tip.label) == 0L,
    setequal(ids, tr$tip.label) && length(ids) == length(tr$tip.label)
  ),
  observed = c(
    length(ids), length(unique(ids)),
    if (length(lengths)) paste(range(lengths), collapse = "-") else "none",
    Ntip(tr), length(unique(tr$tip.label)),
    sprintf("tree_only=%d;fasta_only=%d",
            length(setdiff(tr$tip.label, ids)), length(setdiff(ids, tr$tip.label)))
  ),
  expected = c("495", "495", "182682-182682", "495", "495", "tree_only=0;fasta_only=0"),
  stringsAsFactors = FALSE
)

write.table(checks, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(checks$pass)) stop("Alignment validation failed; see ", out_path)
cat("Alignment validation PASS: 495 unique records, all 182682 nt, exact tree tip set.\n")
