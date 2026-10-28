# Retrieve the 117 Bristol BRIDGE files required by the analysis.
# Run explicitly from the repository root before run_review.R.
source("load_pipeline.R")
cfg <- bristol_config()
bristol_packages("ncdf4")
manifest <- bristol_pipeline_sources(cfg)
manifest <- manifest[order(manifest$stage, manifest$period), ]
manifest$bytes <- NA_real_
manifest$md5 <- NA_character_
manifest$status <- "pending"
dir.create(cfg$raw_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(cfg$raw_dir, "download_manifest.csv")

old_options <- options(timeout = max(600, getOption("timeout", 60)))

for (i in seq_len(nrow(manifest))) {
  target <- manifest$path[i]
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  valid <- file.exists(target) &&
    isTRUE(tryCatch(validate_bristol_file(target), error = function(e) FALSE))
  if (!valid) {
    partial <- paste0(target, ".part")
    if (file.exists(partial)) unlink(partial)
    message(sprintf("Downloading %d/%d: %s", i, nrow(manifest), manifest$filename[i]))
    tryCatch({
      status <- utils::download.file(manifest$url[i], partial,
                                     method = "libcurl", mode = "wb", quiet = TRUE)
      if (status != 0L) stop("Download failed: ", manifest$url[i])
      validate_bristol_file(partial)
      if (!file.copy(partial, target, overwrite = TRUE))
        stop("Could not save: ", target)
      unlink(partial)
    }, error = function(e) {
      if (file.exists(partial)) unlink(partial)
      manifest$status[i] <<- paste("failed:", conditionMessage(e))
    })
  }
  if (file.exists(target) &&
      isTRUE(tryCatch(validate_bristol_file(target), error = function(e) FALSE))) {
    manifest$bytes[i] <- file.info(target)$size
    manifest$md5[i] <- unname(tools::md5sum(target))
    manifest$status[i] <- "complete"
  }
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
}

if (any(manifest$status != "complete"))
  stop("Some climate files are missing or invalid; inspect ", manifest_path,
       " and rerun this script.")
options(old_options)
message("All ", nrow(manifest), " Bristol climate files are ready.")
