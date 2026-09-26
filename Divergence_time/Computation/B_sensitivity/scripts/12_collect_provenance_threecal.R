#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  cat("Usage: 12_collect_provenance_threecal.R <run_root>\n", file = stderr())
  quit(status = 64L)
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

run_root <- normalizePath(args[[1]], mustWork = TRUE)
config_dir <- file.path(run_root, "configs")
log_dir <- file.path(run_root, "logs")
qa_dir <- file.path(run_root, "qa")
input_dir <- file.path(run_root, "input")
software_dir <- file.path(run_root, "software")

required_dirs <- c(config_dir, log_dir, qa_dir, input_dir, software_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs)) {
  fail("Missing required directory/directories: ", paste(missing_dirs, collapse = ", "))
}

artifact_out <- file.path(qa_dir, "artifact_sha256.tsv")
stage_out <- file.path(qa_dir, "stage_summary.tsv")
software_out <- file.path(run_root, "software_versions.txt")
final_outputs <- c(artifact_out, stage_out, software_out)
incomplete_outputs <- paste0(final_outputs, ".incomplete")
existing_outputs <- c(final_outputs, incomplete_outputs)[
  file.exists(c(final_outputs, incomplete_outputs))
]
if (length(existing_outputs)) {
  cat(
    "Refusing to overwrite existing provenance output or incomplete output:\n",
    paste0("  ", existing_outputs, collapse = "\n"),
    "\n",
    file = stderr()
  )
  quit(status = 73L)
}

inside_run <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  startsWith(paste0(resolved, "/"), paste0(run_root, "/"))
}

resolve_run_path <- function(path, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    fail(label, " is empty")
  }
  candidate <- if (startsWith(path, "/")) path else file.path(run_root, path)
  resolved <- normalizePath(candidate, mustWork = TRUE)
  if (!inside_run(resolved)) fail(label, " resolves outside run root: ", resolved)
  resolved
}

relative_path <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(run_root, "/")
  if (!startsWith(resolved, prefix)) fail("Path is outside run root: ", resolved)
  substring(resolved, nchar(prefix) + 1L)
}

sha256_file <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  result <- system2(
    "/usr/bin/sha256sum",
    args = c("--", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || !length(result)) {
    fail("sha256sum failed for ", path, ": ", paste(result, collapse = " | "))
  }
  hash <- sub("[[:space:]].*$", "", result[[1]])
  if (!grepl("^[0-9a-f]{64}$", hash)) fail("Malformed SHA-256 for ", path)
  hash
}

file_bytes <- function(path) {
  value <- file.info(path)$size
  if (length(value) != 1L || is.na(value) || value < 0) {
    fail("Could not determine file size: ", path)
  }
  as.character(value)
}

file_mtime <- function(path) {
  value <- file.info(path)$mtime
  if (length(value) != 1L || is.na(value)) fail("Could not determine mtime: ", path)
  format(value, "%Y-%m-%dT%H:%M:%S%z")
}

parse_config <- function(path) {
  raw <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", raw))
  active <- active[nzchar(active)]
  has_equals <- grepl("=", active, fixed = TRUE)
  keys <- active
  keys[has_equals] <- trimws(sub("=.*$", "", active[has_equals]))
  values <- rep("", length(active))
  values[has_equals] <- trimws(sub("^[^=]*=", "", active[has_equals]))
  list(
    active = active,
    keys = keys,
    values = values,
    get = function(key) values[keys == key]
  )
}

parse_metadata <- function(path) {
  x <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    col.names = c("key", "value"),
    colClasses = "character",
    fill = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(x) != 2L || !nrow(x)) fail("Malformed or empty metadata: ", path)
  required <- c(
    "stage", "pid", "start_time", "binary", "config",
    "end_time", "exit_code", "elapsed_s"
  )
  missing <- setdiff(required, x$key)
  duplicated_required <- required[vapply(required, function(k) sum(x$key == k) != 1L, logical(1))]
  if (length(missing)) fail("Metadata missing keys in ", path, ": ", paste(missing, collapse = ", "))
  if (length(duplicated_required)) {
    fail("Metadata required keys are not unique in ", path, ": ", paste(duplicated_required, collapse = ", "))
  }
  value <- function(key) x$value[x$key == key][[1]]
  list(
    stage = value("stage"),
    pid = value("pid"),
    start_time = value("start_time"),
    binary = value("binary"),
    config = value("config"),
    end_time = value("end_time"),
    exit_code = value("exit_code"),
    elapsed_s = value("elapsed_s")
  )
}

validate_metadata_scalars <- function(meta, path, expected_stage) {
  if (!identical(meta$stage, expected_stage)) {
    fail("Metadata stage mismatch in ", path, ": expected ", expected_stage, ", observed ", meta$stage)
  }
  if (!grepl("^[0-9]+$", meta$pid) || as.numeric(meta$pid) <= 0) {
    fail("Invalid pid in metadata: ", path)
  }
  iso_pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(Z|[+-][0-9]{2}:[0-9]{2})$"
  if (!grepl(iso_pattern, meta$start_time) || !grepl(iso_pattern, meta$end_time)) {
    fail("Invalid or incomplete start/end timestamp in metadata: ", path)
  }
  normalized_start <- sub("Z$", "+0000", meta$start_time)
  normalized_end <- sub("Z$", "+0000", meta$end_time)
  normalized_start <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", normalized_start)
  normalized_end <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", normalized_end)
  start_parsed <- as.POSIXct(normalized_start, format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
  end_parsed <- as.POSIXct(normalized_end, format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
  if (is.na(start_parsed) || is.na(end_parsed)) {
    fail("Could not parse start/end timestamp in metadata: ", path)
  }
  duration_from_times <- as.numeric(difftime(end_parsed, start_parsed, units = "secs"))
  if (duration_from_times < 0) fail("Metadata end_time precedes start_time: ", path)
  if (!grepl("^-?[0-9]+$", meta$exit_code)) fail("Invalid exit_code in metadata: ", path)
  if (as.integer(meta$exit_code) != 0L) {
    fail("Stage did not exit successfully: ", expected_stage, " (exit_code=", meta$exit_code, ")")
  }
  if (!grepl("^[0-9]+$", meta$elapsed_s) || as.numeric(meta$elapsed_s) < 0) {
    fail("Invalid elapsed_s in metadata: ", path)
  }
  if (abs(duration_from_times - as.numeric(meta$elapsed_s)) > 2) {
    fail("Metadata elapsed_s disagrees with start/end timestamps by more than 2 seconds: ", path)
  }
}

stage_name_ok <- function(x) grepl("^[A-Za-z0-9_.-]+$", x)

prime_config <- file.path(config_dir, "prime.cfg")
if (!file.exists(prime_config)) fail("Missing prime config: ", prime_config)

prime_meta_paths <- sort(list.files(
  log_dir,
  pattern = "^prime(_[A-Za-z0-9_.-]+)?\\.metadata\\.tsv$",
  full.names = TRUE
))
if (!file.path(log_dir, "prime.metadata.tsv") %in% prime_meta_paths) {
  fail("Missing required prime metadata: logs/prime.metadata.tsv")
}
prime_stages <- sub("\\.metadata\\.tsv$", "", basename(prime_meta_paths))

cv_configs <- sort(list.files(
  config_dir,
  pattern = "^cv_[A-Za-z0-9_.-]+\\.cfg$",
  full.names = TRUE
))
if (!length(cv_configs)) fail("No CV configs found under configs/cv_*.cfg")
cv_stages <- sub("\\.cfg$", "", basename(cv_configs))
cv_meta_stages <- sub(
  "\\.metadata\\.tsv$", "",
  basename(list.files(
    log_dir,
    pattern = "^cv_[A-Za-z0-9_.-]+\\.metadata\\.tsv$",
    full.names = TRUE
  ))
)
if (!setequal(cv_stages, cv_meta_stages)) {
  fail(
    "CV config/metadata stage sets differ; configs only: ",
    paste(setdiff(cv_stages, cv_meta_stages), collapse = ","),
    "; metadata only: ",
    paste(setdiff(cv_meta_stages, cv_stages), collapse = ",")
  )
}

final_manifest_path <- file.path(qa_dir, "final_variant_manifest.tsv")
if (!file.exists(final_manifest_path)) {
  fail(
    "Missing qa/final_variant_manifest.tsv. Run this collector only after CV selection, ",
    "all final variants, and final-manifest creation are complete."
  )
}
final_manifest <- read.delim(
  final_manifest_path,
  sep = "\t",
  quote = "",
  comment.char = "",
  colClasses = "character",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
manifest_required <- c(
  "variant_id", "variant_role", "smooth", "thorough", "seed",
  "config_path", "tree_path", "stage"
)
missing_manifest_columns <- setdiff(manifest_required, names(final_manifest))
if (length(missing_manifest_columns)) {
  fail("Final manifest is missing columns: ", paste(missing_manifest_columns, collapse = ", "))
}
if (!nrow(final_manifest)) fail("Final manifest has no rows")
if (anyNA(final_manifest[manifest_required]) || any(!nzchar(as.matrix(final_manifest[manifest_required])))) {
  fail("Final manifest contains an empty required field")
}
if (anyDuplicated(final_manifest$stage)) fail("Final manifest stage values are not unique")
if (anyDuplicated(final_manifest$variant_id)) fail("Final manifest variant_id values are not unique")
if (any(!vapply(final_manifest$stage, stage_name_ok, logical(1)))) {
  fail("Final manifest contains an unsafe stage name")
}

stage_specs <- list()
for (i in seq_along(prime_stages)) {
  stage_specs[[length(stage_specs) + 1L]] <- list(
    stage = prime_stages[[i]],
    group = "prime",
    config = normalizePath(prime_config, mustWork = TRUE),
    manifest_tree = NA_character_,
    variant_id = if (identical(prime_stages[[i]], "prime")) "prime" else prime_stages[[i]]
  )
}
for (i in seq_along(cv_stages)) {
  stage_specs[[length(stage_specs) + 1L]] <- list(
    stage = cv_stages[[i]],
    group = "cv",
    config = normalizePath(cv_configs[[i]], mustWork = TRUE),
    manifest_tree = NA_character_,
    variant_id = cv_stages[[i]]
  )
}
for (i in seq_len(nrow(final_manifest))) {
  stage_specs[[length(stage_specs) + 1L]] <- list(
    stage = final_manifest$stage[[i]],
    group = "final",
    config = resolve_run_path(final_manifest$config_path[[i]], "final manifest config_path"),
    manifest_tree = resolve_run_path(final_manifest$tree_path[[i]], "final manifest tree_path"),
    variant_id = final_manifest$variant_id[[i]]
  )
}

all_stages <- vapply(stage_specs, `[[`, character(1), "stage")
if (anyDuplicated(all_stages)) fail("Stage names are duplicated across prime, CV, and final groups")
if (any(!vapply(all_stages, stage_name_ok, logical(1)))) fail("Unsafe stage name discovered")

collected_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
artifact_rows <- list()
stage_rows <- list()

add_artifact <- function(stage, group, category, role, path, timing_note) {
  if (!file.exists(path) || dir.exists(path)) fail("Missing required stage artifact: ", path)
  bytes <- file_bytes(path)
  must_be_nonempty <- category %in% c("config", "binary", "input", "output", "metadata", "control", "software") || role == "resource_usage"
  if (must_be_nonempty && as.numeric(bytes) == 0) {
    fail("Required artifact is empty: ", path)
  }
  artifact_rows[[length(artifact_rows) + 1L]] <<- data.frame(
    stage = stage,
    stage_group = group,
    artifact_category = category,
    artifact_role = role,
    path = relative_path(path),
    sha256 = sha256_file(path),
    bytes = bytes,
    mtime = file_mtime(path),
    hash_collected_at = collected_at,
    hash_timing = "post_run_collection_snapshot",
    launch_time_hash = "NO",
    interpretation = timing_note,
    stringsAsFactors = FALSE
  )
  artifact_rows[[length(artifact_rows)]]
}

for (spec in stage_specs) {
  stage <- spec$stage
  group <- spec$group
  meta_path <- file.path(log_dir, paste0(stage, ".metadata.tsv"))
  stdout_path <- file.path(log_dir, paste0(stage, ".stdout"))
  stderr_path <- file.path(log_dir, paste0(stage, ".stderr"))
  resource_path <- file.path(log_dir, paste0(stage, ".resources.txt"))
  required_stage_files <- c(meta_path, stdout_path, stderr_path, resource_path, spec$config)
  missing_stage_files <- required_stage_files[!file.exists(required_stage_files)]
  if (length(missing_stage_files)) {
    fail("Stage ", stage, " is missing required files: ", paste(missing_stage_files, collapse = ", "))
  }

  meta <- parse_metadata(meta_path)
  validate_metadata_scalars(meta, meta_path, stage)
  meta_binary <- resolve_run_path(meta$binary, paste0(stage, " metadata binary"))
  meta_config <- resolve_run_path(meta$config, paste0(stage, " metadata config"))
  if (!identical(meta_config, spec$config)) {
    fail("Stage ", stage, " metadata config does not match enumerated config")
  }

  cfg <- parse_config(spec$config)
  treefiles <- cfg$get("treefile")
  if (length(treefiles) != 1L) fail("Stage ", stage, " config must contain exactly one treefile")
  input_tree <- resolve_run_path(treefiles[[1]], paste0(stage, " treefile"))
  outfile_values <- cfg$get("outfile")
  cvoutfile_values <- cfg$get("cvoutfile")

  if (identical(group, "prime")) {
    if (length(outfile_values) || length(cvoutfile_values)) {
      fail("Prime config unexpectedly contains outfile or cvoutfile: ", spec$config)
    }
    primary_output <- stdout_path
    configured_outputs <- character()
  } else if (identical(group, "cv")) {
    if (length(outfile_values) != 1L || length(cvoutfile_values) != 1L) {
      fail("CV config must contain exactly one outfile and one cvoutfile: ", spec$config)
    }
    configured_outputs <- c(
      resolve_run_path(outfile_values[[1]], paste0(stage, " outfile")),
      resolve_run_path(cvoutfile_values[[1]], paste0(stage, " cvoutfile"))
    )
    primary_output <- configured_outputs[[1]]
  } else {
    if (length(outfile_values) != 1L || length(cvoutfile_values)) {
      fail("Final config must contain exactly one outfile and no cvoutfile: ", spec$config)
    }
    configured_outputs <- resolve_run_path(outfile_values[[1]], paste0(stage, " outfile"))
    if (!identical(configured_outputs[[1]], spec$manifest_tree)) {
      fail("Final manifest tree_path does not match config outfile for stage ", stage)
    }
    primary_output <- configured_outputs[[1]]
  }

  timing_note <- if (identical(group, "cv")) {
    "CV artifact hash collected after successful completion; not a launch-time SHA-256 attestation"
  } else if (identical(group, "final")) {
    "Final-stage artifact hash collected after successful completion; not a launch-time SHA-256 attestation"
  } else {
    "Prime artifact hash collected after successful completion; not a launch-time SHA-256 attestation"
  }

  config_art <- add_artifact(stage, group, "config", "treePL_config", spec$config, timing_note)
  binary_art <- add_artifact(stage, group, "binary", "treePL_executable", meta_binary, timing_note)
  input_art <- add_artifact(stage, group, "input", "input_tree", input_tree, timing_note)
  if (identical(group, "prime")) {
    stdout_art <- add_artifact(
      stage, group, "output", "prime_stdout_result", stdout_path, timing_note
    )
    primary_art <- stdout_art
  } else {
    primary_art <- add_artifact(
      stage, group, "output", "dated_tree", configured_outputs[[1]], timing_note
    )
    if (identical(group, "cv")) {
      add_artifact(stage, group, "output", "cv_score_table", configured_outputs[[2]], timing_note)
    }
    auxiliary_r8s <- paste0(configured_outputs[[1]], ".r8s")
    if (file.exists(auxiliary_r8s)) {
      add_artifact(stage, group, "output", "treePL_r8s_auxiliary", auxiliary_r8s, timing_note)
    }
    stdout_art <- add_artifact(stage, group, "log", "stdout", stdout_path, timing_note)
  }
  stderr_art <- add_artifact(stage, group, "log", "stderr", stderr_path, timing_note)
  resource_art <- add_artifact(stage, group, "log", "resource_usage", resource_path, timing_note)
  metadata_art <- add_artifact(stage, group, "metadata", "wrapper_metadata", meta_path, timing_note)

  configured_output_rel <- if (length(configured_outputs)) {
    paste(vapply(configured_outputs, relative_path, character(1)), collapse = ";")
  } else {
    "NA (prime result is stdout)"
  }
  configured_output_hashes <- if (length(configured_outputs)) {
    paste(vapply(configured_outputs, sha256_file, character(1)), collapse = ";")
  } else {
    "NA (prime result is stdout)"
  }
  stage_rows[[length(stage_rows) + 1L]] <- data.frame(
    stage = stage,
    stage_group = group,
    variant_id = spec$variant_id,
    status = "PASS",
    pid = meta$pid,
    start_time = meta$start_time,
    end_time = meta$end_time,
    elapsed_s = meta$elapsed_s,
    exit_code = meta$exit_code,
    config_path = config_art$path[[1]],
    config_sha256 = config_art$sha256[[1]],
    binary_path = binary_art$path[[1]],
    binary_sha256 = binary_art$sha256[[1]],
    input_tree_path = input_art$path[[1]],
    input_tree_sha256 = input_art$sha256[[1]],
    configured_output_paths = configured_output_rel,
    configured_output_sha256 = configured_output_hashes,
    primary_result_path = primary_art$path[[1]],
    primary_result_sha256 = primary_art$sha256[[1]],
    stdout_path = stdout_art$path[[1]],
    stdout_sha256 = stdout_art$sha256[[1]],
    stderr_path = stderr_art$path[[1]],
    stderr_sha256 = stderr_art$sha256[[1]],
    resource_log_path = resource_art$path[[1]],
    resource_log_sha256 = resource_art$sha256[[1]],
    metadata_path = metadata_art$path[[1]],
    metadata_sha256 = metadata_art$sha256[[1]],
    hash_collected_at = collected_at,
    hash_timing = "post_run_collection_snapshot",
    launch_time_hashes_recorded = "NO",
    hash_interpretation = timing_note,
    stringsAsFactors = FALSE
  )
}

# Hash the complete analysis input directory once, including the alignment that
# supplies numsites but is not read directly by treePL.
global_input_files <- sort(list.files(input_dir, full.names = TRUE, recursive = TRUE))
global_input_files <- global_input_files[file.exists(global_input_files) & !dir.exists(global_input_files)]
if (!length(global_input_files)) fail("No global input files found under input/")
for (path in global_input_files) {
  add_artifact(
    "GLOBAL", "global", "input", "archived_analysis_input", path,
    "Global input hash collected after all successful stages; not a launch-time SHA-256 attestation"
  )
}

global_control_files <- c(
  final_manifest_path,
  file.path(run_root, "analysis_design.md"),
  file.path(run_root, "scripts", "run_treepl_stage.sh"),
  file.path(run_root, "scripts", "12_collect_provenance_threecal.R")
)
missing_control_files <- global_control_files[!file.exists(global_control_files)]
if (length(missing_control_files)) {
  fail("Missing required analysis-control file(s): ", paste(missing_control_files, collapse = ", "))
}
global_control_files <- unique(global_control_files)
for (path in global_control_files) {
  add_artifact(
    "GLOBAL", "global", "control", "analysis_control", path,
    "Control-file hash collected at provenance collection time"
  )
}

patch_path <- file.path(software_dir, "randomcv_parent_index_fix.patch")
if (!file.exists(patch_path)) fail("Missing archived random-CV patch: ", patch_path)
add_artifact(
  "GLOBAL", "global", "software", "randomcv_source_patch", patch_path,
  "Patch hash collected at provenance collection time"
)

artifact_table <- do.call(rbind, artifact_rows)
stage_table <- do.call(rbind, stage_rows)
artifact_table <- artifact_table[order(
  match(artifact_table$stage_group, c("prime", "cv", "final", "global")),
  artifact_table$stage,
  artifact_table$artifact_category,
  artifact_table$artifact_role,
  artifact_table$path
), , drop = FALSE]
stage_table <- stage_table[order(
  match(stage_table$stage_group, c("prime", "cv", "final")),
  stage_table$stage
), , drop = FALSE]

if (any(stage_table$exit_code != "0") || any(stage_table$status != "PASS")) {
  fail("Internal error: non-passing stage reached output construction")
}
if (anyDuplicated(paste(artifact_table$stage, artifact_table$artifact_role, artifact_table$path, sep = "\t"))) {
  fail("Internal error: duplicate stage/role/path artifact row")
}

unique_binaries <- unique(vapply(stage_specs, function(spec) {
  meta <- parse_metadata(file.path(log_dir, paste0(spec$stage, ".metadata.tsv")))
  resolve_run_path(meta$binary, paste0(spec$stage, " metadata binary"))
}, character(1)))

capture_first <- function(command, arguments = character()) {
  result <- suppressWarnings(system2(command, arguments, stdout = TRUE, stderr = TRUE))
  if (!length(result)) return("not reported")
  gsub("[\t\r\n]+", " ", result[[1]])
}

source_version_path <- file.path(run_root, "software_versions.source_two_cal_run.txt")
treepl_reported_version <- "not independently queried; run-local binary SHA-256 is authoritative"
if (file.exists(source_version_path)) {
  source_lines <- readLines(source_version_path, warn = FALSE)
  version_lines <- trimws(source_lines[grepl("^reported version:", trimws(source_lines))])
  if (length(version_lines)) {
    treepl_reported_version <- sub("^reported version:[[:space:]]*", "", version_lines[[1]])
  }
}

wrapper_path <- file.path(run_root, "scripts", "run_treepl_stage.sh")
collector_path <- file.path(run_root, "scripts", "12_collect_provenance_threecal.R")
software_lines <- c(
  "treePL three-calibration analysis — software and provenance snapshot",
  paste0("run_root\t", run_root),
  paste0("collected_at\t", collected_at),
  "hash_timing\tpost_run_collection_snapshot",
  "launch_time_hashes_recorded\tNO",
  paste0(
    "hash_scope_note\t",
    "Hashes were computed only after successful stage completion. CV and final metadata record ",
    "paths/times/exit status but no launch-time SHA-256 values; therefore these are post-run ",
    "snapshots, not cryptographic proof of the exact bytes present at launch."
  ),
  paste0("treePL_reported_version\t", treepl_reported_version),
  paste0("treePL_binary_count\t", length(unique_binaries)),
  unlist(lapply(seq_along(unique_binaries), function(i) c(
    paste0("treePL_binary_", i, "_path\t", unique_binaries[[i]]),
    paste0("treePL_binary_", i, "_sha256\t", sha256_file(unique_binaries[[i]])),
    paste0("treePL_binary_", i, "_bytes\t", file_bytes(unique_binaries[[i]]))
  ))),
  paste0("randomcv_patch_path\t", relative_path(patch_path)),
  paste0("randomcv_patch_sha256\t", sha256_file(patch_path)),
  paste0("R_version\t", R.version.string),
  paste0("R_platform\t", R.version$platform),
  paste0(
    "ape_version\t",
    if (requireNamespace("ape", quietly = TRUE)) as.character(utils::packageVersion("ape")) else "not installed"
  ),
  paste0("bash_version\t", capture_first("/usr/bin/bash", "--version")),
  paste0("sha256sum_version\t", capture_first("/usr/bin/sha256sum", "--version")),
  paste0("kernel\t", capture_first("/usr/bin/uname", c("-s", "-r", "-m"))),
  paste0("wrapper_path\t", relative_path(wrapper_path)),
  paste0("wrapper_sha256\t", sha256_file(wrapper_path)),
  paste0("collector_path\t", relative_path(collector_path)),
  paste0("collector_sha256\t", sha256_file(collector_path)),
  paste0(
    "archived_source_software_report\t",
    if (file.exists(source_version_path)) relative_path(source_version_path) else "not present"
  ),
  paste0(
    "archived_source_software_report_sha256\t",
    if (file.exists(source_version_path)) sha256_file(source_version_path) else "NA"
  ),
  paste0("stages_verified\t", nrow(stage_table)),
  paste0("prime_stages_verified\t", sum(stage_table$stage_group == "prime")),
  paste0("cv_stages_verified\t", sum(stage_table$stage_group == "cv")),
  paste0("final_stages_verified\t", sum(stage_table$stage_group == "final"))
)

write.table(
  artifact_table,
  file = paste0(artifact_out, ".incomplete"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)
write.table(
  stage_table,
  file = paste0(stage_out, ".incomplete"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)
writeLines(software_lines, paste0(software_out, ".incomplete"), useBytes = TRUE)

for (i in seq_along(final_outputs)) {
  if (!file.rename(incomplete_outputs[[i]], final_outputs[[i]])) {
    fail("Could not finalize provenance output: ", final_outputs[[i]])
  }
}

cat(
  "Provenance collection PASS:", nrow(stage_table), "stages and",
  nrow(artifact_table), "artifact rows\n"
)
