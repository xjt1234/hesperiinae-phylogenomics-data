#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_preflight_mapping.R <input.treefile> <run_dir>")
}

input_tree <- normalizePath(args[[1]], mustWork = TRUE)
run_dir <- normalizePath(args[[2]], mustWork = TRUE)
mapping_dir <- file.path(run_dir, "mapping")
qa_dir <- file.path(run_dir, "qa")
dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(input_tree)
n_tip <- Ntip(tree)
n_node <- Nnode(tree)
if (is.null(tree$edge.length)) stop("Input tree has no branch lengths")

root_nodes <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
if (length(root_nodes) != 1L) stop("Could not identify exactly one root node")
root_node <- root_nodes[[1]]

descendant_tips <- function(node) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0L) {
    if (node <= n_tip) return(tree$tip.label[[node]])
    return(character())
  }
  sort(unique(unlist(lapply(children, descendant_tips), use.names = FALSE)))
}

node_depth <- node.depth.edgelength(tree)
tip_depth <- node_depth[seq_len(n_tip)]

papilionidae_tips <- sort(c(
  "Parnassius_apollo_ncbi",
  "Parnassius_glacialis_ncbi",
  "Graphium_cloanthus_mydata",
  "Graphium_sarpedon_mydata",
  "Papilio_helenus_ncbi",
  "Papilio_machaon_ncbi",
  "Papilio_xuthus_ncbi"
))

pap_root_pair <- c("Parnassius_apollo_ncbi", "Aeromachus_catocyanea_mydata")
pap_family_pair <- c("Parnassius_apollo_ncbi", "Graphium_cloanthus_mydata")
hesp_pair <- c("Aeromachus_catocyanea_mydata", "Acada_biseriata_kawahara2023")

all_required_labels <- unique(c(papilionidae_tips, pap_root_pair, pap_family_pair, hesp_pair))
missing_labels <- setdiff(all_required_labels, tree$tip.label)
if (length(missing_labels)) {
  stop("Required labels missing from input tree: ", paste(missing_labels, collapse = ", "))
}

root_children <- tree$edge[tree$edge[, 1] == root_node, 2]
root_child_sets <- lapply(root_children, descendant_tips)
root_child_sizes <- lengths(root_child_sets)

pap_node <- getMRCA(tree, pap_family_pair)
root_mrca <- getMRCA(tree, pap_root_pair)
hesp_tips <- sort(setdiff(tree$tip.label, papilionidae_tips))
hesp_node <- getMRCA(tree, hesp_pair)

pap_desc <- descendant_tips(pap_node)
root_desc <- descendant_tips(root_mrca)
hesp_desc <- descendant_tips(hesp_node)
hesp_children <- tree$edge[tree$edge[, 1] == hesp_node, 2]
hesp_child_sets <- lapply(hesp_children, descendant_tips)
hesp_child_sizes <- lengths(hesp_child_sets)
hesp_anchor_child_index <- vapply(
  hesp_pair,
  function(anchor) which(vapply(hesp_child_sets, function(x) anchor %in% x, logical(1))),
  integer(1)
)

checks <- data.frame(
  check = c(
    "tip_count_495",
    "unique_tip_labels",
    "rooted",
    "binary",
    "finite_branch_lengths",
    "nonnegative_branch_lengths",
    "input_is_not_ultrametric",
    "root_has_two_children",
    "root_split_7_488",
    "papilionidae_expected_labels_present",
    "papilionidae_anchor_descendants_exact_7",
    "papilionidae_is_root_child",
    "root_anchor_is_root",
    "root_anchor_descendants_exact_495",
    "hesperiinae_anchor_descendants_exact_488",
    "hesperiinae_is_root_child",
    "hesperiinae_has_two_children",
    "hesperiinae_child_split_38_450",
    "hesperiinae_anchors_span_crown_split",
    "baronia_absent"
  ),
  pass = c(
    n_tip == 495L,
    anyDuplicated(tree$tip.label) == 0L,
    is.rooted(tree),
    is.binary.phylo(tree),
    all(is.finite(tree$edge.length)),
    all(tree$edge.length >= 0),
    !is.ultrametric(tree, tol = 1e-8),
    length(root_children) == 2L,
    identical(sort(root_child_sizes), c(7L, 488L)),
    all(papilionidae_tips %in% tree$tip.label),
    identical(pap_desc, papilionidae_tips),
    pap_node %in% root_children,
    identical(root_mrca, root_node),
    identical(root_desc, sort(tree$tip.label)),
    identical(hesp_desc, hesp_tips),
    hesp_node %in% root_children,
    length(hesp_children) == 2L,
    identical(sort(hesp_child_sizes), c(38L, 450L)),
    length(unique(hesp_anchor_child_index)) == 2L,
    !any(grepl("Baronia", tree$tip.label, ignore.case = TRUE))
  ),
  observed = c(
    as.character(n_tip),
    as.character(length(unique(tree$tip.label))),
    as.character(is.rooted(tree)),
    as.character(is.binary.phylo(tree)),
    as.character(all(is.finite(tree$edge.length))),
    sprintf("%.12g", min(tree$edge.length)),
    sprintf("range=%.12g", max(tip_depth) - min(tip_depth)),
    as.character(length(root_children)),
    paste(sort(root_child_sizes), collapse = ","),
    as.character(sum(papilionidae_tips %in% tree$tip.label)),
    as.character(length(pap_desc)),
    as.character(pap_node %in% root_children),
    as.character(root_mrca),
    as.character(length(root_desc)),
    as.character(length(hesp_desc)),
    as.character(hesp_node %in% root_children),
    as.character(length(hesp_children)),
    paste(sort(hesp_child_sizes), collapse = ","),
    paste(hesp_anchor_child_index, collapse = ","),
    as.character(sum(grepl("Baronia", tree$tip.label, ignore.case = TRUE)))
  ),
  expected = c(
    "495", "495", "TRUE", "TRUE", "TRUE", ">=0", "non-ultrametric",
    "2", "7,488", "7", "7 exact labels", "TRUE", as.character(root_node),
    "495 exact labels", "488 exact labels", "TRUE", "2", "38,450",
    "different child clades", "0"
  ),
  stringsAsFactors = FALSE
)

write.table(
  checks,
  file.path(qa_dir, "preflight_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

stats <- data.frame(
  metric = c(
    "input_tree", "tips", "internal_nodes", "edges", "root_node",
    "root_child_sizes", "branch_length_min", "branch_length_max",
    "branch_length_sum", "root_to_tip_min", "root_to_tip_max",
    "root_to_tip_range", "ultrametric_tol_1e-8"
  ),
  value = c(
    input_tree, n_tip, n_node, nrow(tree$edge), root_node,
    paste(sort(root_child_sizes), collapse = ","),
    sprintf("%.15g", min(tree$edge.length)),
    sprintf("%.15g", max(tree$edge.length)),
    sprintf("%.15g", sum(tree$edge.length)),
    sprintf("%.15g", min(tip_depth)),
    sprintf("%.15g", max(tip_depth)),
    sprintf("%.15g", max(tip_depth) - min(tip_depth)),
    as.character(is.ultrametric(tree, tol = 1e-8))
  ),
  stringsAsFactors = FALSE
)
write.table(stats, file.path(qa_dir, "input_tree_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mapping <- data.frame(
  node_id = c("CROWN_PAPILIONOIDEA", "CROWN_PAPILIONIDAE_NONBARONIA", "CROWN_HESPERIINAE"),
  biological_definition = c(
    "crown Papilionoidea (total-tree root)",
    "crown Papilionidae excluding Baronia",
    "crown Hesperiinae"
  ),
  anchor_tip_1 = c(pap_root_pair[[1]], pap_family_pair[[1]], hesp_pair[[1]]),
  anchor_tip_2 = c(pap_root_pair[[2]], pap_family_pair[[2]], hesp_pair[[2]]),
  mrca_node = c(root_mrca, pap_node, hesp_node),
  is_root = c(root_mrca == root_node, pap_node == root_node, hesp_node == root_node),
  descendant_tip_count = c(length(root_desc), length(pap_desc), length(hesp_desc)),
  calibrated = c("yes", "yes", "yes"),
  lower_bound_ma = c(91.5046, 44.1968, 36.210861),
  upper_bound_ma = c(100.8925, 52.9473, 40.662537),
  descendant_set_exact = c(
    identical(root_desc, sort(tree$tip.label)),
    identical(pap_desc, papilionidae_tips),
    identical(hesp_desc, hesp_tips)
  ),
  stringsAsFactors = FALSE
)
write.table(
  mapping,
  file.path(mapping_dir, "calibration_node_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

con <- file(file.path(mapping_dir, "calibration_descendants.txt"), open = "wt")
on.exit(close(con), add = TRUE)
writeLines(c(
  "[CROWN_PAPILIONOIDEA]",
  sprintf("MRCA node: %d; descendant tips: %d", root_mrca, length(root_desc)),
  root_desc,
  "",
  "[CROWN_PAPILIONIDAE_NONBARONIA]",
  sprintf("MRCA node: %d; descendant tips: %d", pap_node, length(pap_desc)),
  pap_desc,
  "",
  "[CROWN_HESPERIINAE]",
  sprintf("MRCA node: %d; descendant tips: %d", hesp_node, length(hesp_desc)),
  hesp_desc
), con)
close(con)
on.exit(NULL, add = FALSE)

if (!all(checks$pass)) {
  failed <- checks$check[!checks$pass]
  stop("Preflight hard gate failed: ", paste(failed, collapse = ", "))
}

cat(sprintf("Preflight PASS: %d/%d checks\n", sum(checks$pass), nrow(checks)))
cat(sprintf("Root node %d: %d descendants\n", root_mrca, length(root_desc)))
cat(sprintf("Papilionidae node %d: %d descendants\n", pap_node, length(pap_desc)))
cat(sprintf(
  "Hesperiinae node %d: %d descendants; calibrated 36.210861-40.662537 Ma\n",
  hesp_node, length(hesp_desc)
))
