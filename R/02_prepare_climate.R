# Fork of prepare_paleoclim_dat.R. Native stages; no temporal interpolation.
bristol_pipeline_sources <- function(cfg = bristol_config()) {
  s <- read_bristol_stages(cfg)
  runs <- rbind(data.frame(stage = cfg$baseline[["stage"]], experiment = cfg$baseline[["experiment"]]), s[, c("stage", "experiment")])
  m <- merge(runs, data.frame(period = c("ann", tolower(month.abb))), sort = FALSE)
  m$filename <- paste0(m$experiment, "a.pdcl", m$period, ".nc")
  m$url <- paste0("https://www.paleo.bristol.ac.uk/ummodel/data/", m$experiment, "/climate/", m$filename)
  m$path <- file.path(cfg$raw_dir, "raw", m$experiment, m$filename)
  m
}

# Explicit optional network step. Analysis never downloads automatically.
download_bristol_baseline <- function(cfg = bristol_config()) {
  bristol_packages("ncdf4")
  m <- bristol_pipeline_sources(cfg)
  m <- m[m$experiment == cfg$baseline[["experiment"]], ]
  old <- options(timeout = max(600, getOption("timeout", 60)))
  on.exit(options(old), add = TRUE)
  m$status <- "pending"
  m$error <- ""
  m$md5 <- NA_character_
  dir.create(cfg$raw_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(m))) {
    partial <- paste0(m$path[i], ".part")
    tryCatch({
      dir.create(dirname(m$path[i]), recursive = TRUE, showWarnings = FALSE)
      if (!file.exists(m$path[i])) {
        status <- utils::download.file(m$url[i], partial, method = "libcurl", mode = "wb")
        if (status != 0L) stop("Download failed")
        validate_bristol_file(partial)
        if (!file.rename(partial, m$path[i])) stop("Cannot finalize download")
      } else validate_bristol_file(m$path[i])
      m$md5[i] <- unname(tools::md5sum(m$path[i]))
      m$status[i] <- "complete"
    }, error = function(e) {
      m$status[i] <<- "failed"
      m$error[i] <<- conditionMessage(e)
    }, finally = { if (file.exists(partial)) unlink(partial) })
    utils::write.csv(m, file.path(cfg$raw_dir, "chattian_download_manifest.csv"), row.names = FALSE)
  }
  if (any(m$status != "complete")) warning("Some baseline files failed; inspect chattian_download_manifest.csv")
  invisible(m)
}

prepare_bristol_climate <- function(cfg = bristol_config()) {
  bristol_packages(c("terra", "ncdf4"))
  s <- read_bristol_stages(cfg)
  output <- bristol_output(cfg, "climate_catalog.rds")
  m <- bristol_pipeline_sources(cfg)
  if (any(!file.exists(m$path))) stop("Missing sources (no automatic downloads): ", paste(m$path[!file.exists(m$path)], collapse = ", "))
  m$md5 <- unname(tools::md5sum(m$path))
  utils::write.csv(m, file.path(cfg$out, "climate_source_manifest.csv"), row.names = FALSE)
  experiments <- c(cfg$baseline[["experiment"]], s$experiment)
  variables <- c(temperature = "temp_mm_1_5m", precipitation = "precip_mm_srf")
  catalog <- list()
  units_audit <- list()
  message("[2/5] Reading Chattian baseline and eight native-stage climates")
  if (is.null(cfg$regions)) stop("Old configuration object: recreate cfg with bristol_config()")
  for (region in names(cfg$regions)) {
    catalog[[region]] <- list()
    for (variable in names(variables)) {
      read_period <- function(experiment, period) {
        path <- m$path[m$experiment == experiment & m$period == period]
        stopifnot(length(path) == 1L)
        decoded <- bristol_validate_decoded(read_bristol_variable(path, variables[[variable]]), path)
        audit_key <- paste(experiment, period, variable, sep = "_")
        units_audit[[audit_key]] <<- data.frame(experiment = experiment, period = period,
          variable = variable, source_file = path, source_units = decoded$source_units,
          output_units = decoded$output_units, conversion_note = decoded$conversion_note)
        r <- decoded$raster
        terra::crop(r, terra::ext(cfg$regions[[region]]), snap = "out")
      }
      previous_ann <- read_period(experiments[1], "ann")
      previous_mon <- do.call(c, lapply(tolower(month.abb), function(p) read_period(experiments[1], p)))
      layers <- bands <- list()
      for (i in 1:8) {
        current_ann <- read_period(experiments[i + 1], "ann")
        current_mon <- do.call(c, lapply(tolower(month.abb), function(p) read_period(experiments[i + 1], p)))
        for (r in list(current_ann, current_mon[[1]], previous_mon[[1]])) {
          if (!isTRUE(terra::compareGeom(previous_ann, r, stopOnError = FALSE))) stop("Stage grids differ; no implicit resampling")
        }
        seasonal <- terra::app(current_mon - previous_mon, function(v) {
          if (any(!is.finite(v))) return(NA_real_)
          sqrt(mean((v - mean(v))^2))
        })
        part <- c(current_ann, current_ann - previous_ann, seasonal)
        types <- c("mean", "signed_change", "seasonal_change")
        names(part) <- paste(s$stage[i], types, sep = "_")
        layers[[i]] <- part
        bands[[i]] <- data.frame(band = (i - 1)*3 + 1:3, stage_id = i, stage = s$stage[i],
                                 baseline_stage = if (i == 1) cfg$baseline[["stage"]] else s$stage[i - 1],
                                 metric = types, variable = variable, units = if (variable == "temperature") "degC" else "mm/day")
        previous_ann <- current_ann
        previous_mon <- current_mon
      }
      r <- do.call(c, layers)
      stats <- terra::global(r, c("min", "max"), na.rm = TRUE)
      if (any(!is.finite(as.matrix(stats)))) stop("Empty/non-finite climate layer")
      path <- file.path(cfg$out, paste0(region, "_", variable, "_stage_predictors.tif"))
      terra::writeRaster(r, path, overwrite = cfg$overwrite, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
      bt <- do.call(rbind, bands)
      bt$minimum <- stats$min
      bt$maximum <- stats$max
      utils::write.csv(bt, sub("[.]tif$", "_bands.csv", path), row.names = FALSE)
      catalog[[region]][[variable]] <- path
    }
  }
  utils::write.csv(do.call(rbind, units_audit), file.path(cfg$out, "climate_units_audit.csv"), row.names = FALSE)
  saveRDS(list(catalog = catalog, stages = s, config = cfg), output)
  invisible(catalog)
}
