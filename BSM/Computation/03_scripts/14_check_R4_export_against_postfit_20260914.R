#!/usr/bin/env Rscript
# Independent direct-object check; no likelihood or optimization is evaluated.
suppressPackageStartupMessages({library(digest); library(jsonlite)})
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)))
job <- dirname(dirname(self))
out <- file.path(job,"08_qa/direct_postfit_export_check_20260914.json")
stopifnot(!file.exists(out))
pins <- c(M0="fe80e8e5791d501536856f1b4c22e6eafcaa66f3fa75b8d9ab232e698bc0e542",
          M1="566b86c007d123be0fccf1064abf076ad26e6c082ceac8133cbed447ef52e08d",
          M2="3e23818be5d98b3150e3be09fbfdaa1e1cf1c178c34415cc8349e80db1ad49e2")
hash <- function(p) digest(file=p,algo="sha256",serialize=FALSE)
results <- list()
for (model in names(pins)) {
  directory <- file.path(job,sprintf("04_runs/%s_R4_postfit_v1",model))
  postfile <- file.path(directory,"postfit_recalculated_ancestral_states.rds")
  matricesfile <- file.path(directory,"ancestral_probability_matrices.rds")
  tsvfile <- file.path(job,sprintf("07_tables/ancestral_R4_20260914_v1/%s_node_posteriors.tsv",model))
  stopifnot(hash(postfile)==pins[[model]])
  post <- readRDS(postfile)
  expected <- post$ML_marginal_prob_each_state_at_branch_top_AT_node
  matrices <- readRDS(matricesfile)
  tsv <- read.delim(tsvfile,check.names=FALSE)
  observed <- as.matrix(tsv[,-1L])
  stopifnot(identical(dim(expected),c(833L,562L)),
            identical(dim(observed),dim(expected)),
            identical(tsv$ape_node,1:833),
            identical(matrices$top,expected),
            all(is.finite(observed)),max(abs(observed-expected))<1e-14)
  results[[model]] <- list(status="PASS",postfit_sha256=hash(postfile),
    separate_matrices_sha256=hash(matricesfile),export_sha256=hash(tsvfile),
    matrix_object_identical=TRUE,max_export_absolute_difference=max(abs(observed-expected)),
    node_count=nrow(expected),state_count=ncol(expected))
  post <- matrices <- expected <- observed <- tsv <- NULL
  invisible(gc())
}
result <- list(status="PASS_DIRECT_POSTFIT_TO_FIGURE_MATRIX_IDENTITY",
  checked_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
  script_sha256=hash(self),new_optimizations=0L,new_likelihood_evaluations=0L,
  models=results,scientific_scope="Technical export identity; original KKT1 remains FALSE")
write_json(result,out,auto_unbox=TRUE,pretty=TRUE,digits=17)
cat(toJSON(result,auto_unbox=TRUE,pretty=TRUE,digits=17),"\n")
