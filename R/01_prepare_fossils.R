# Fork of numerical_analyses/prepare_fossil_data.R and the spatial binning
# portion of src/functions_data_prep.R. No source-time analysis.
prepare_bristol_fossils <- function(cfg = bristol_config()) {
  bristol_packages(c("sf", "palaeoverse", "dplyr"))
  stages <- read_bristol_stages(cfg)
  output <- bristol_output(cfg, "stage_fossils.rds")
  message("[1/5] Reading occurrences and assigning midpoint stages")
  x <- utils::read.csv(cfg$fossil_file, stringsAsFactors = FALSE, na.strings = c("", "NA", "\\N"))
  required <- c("GENUS", "ORDER", "MIN_AGE", "MAX_AGE", "LONG", "LAT", "LIDNUM",
    "NOW_DIET_2_GENUS", "diet_status", "trophic_guild")
  stopifnot(all(required %in% names(x)))
  x$occurrence_id <- seq_len(nrow(x))
  for (v in c("MIN_AGE", "MAX_AGE", "LONG", "LAT")) x[[v]] <- suppressWarnings(as.numeric(x[[v]]))
  x$exclusion <- NA_character_
  x$exclusion[is.na(x$GENUS) | trimws(x$GENUS) %in% c("", "indet", "indet.", "gen.", "Gen.")] <- "indeterminate_genus"
  bad_age <- !is.finite(x$MIN_AGE) | !is.finite(x$MAX_AGE) | x$MIN_AGE > x$MAX_AGE | x$MIN_AGE < 0
  x$exclusion[is.na(x$exclusion) & bad_age] <- "invalid_age_range"
  bad_xy <- !is.finite(x$LONG) | !is.finite(x$LAT) | abs(x$LONG) > 180 | abs(x$LAT) > 90
  x$exclusion[is.na(x$exclusion) & bad_xy] <- "invalid_coordinates"
  x$midpoint_ma <- (x$MIN_AGE + x$MAX_AGE)/2
  x$stage_id <- NA_integer_
  # Internal boundaries belong to the younger stage. The oldest outer edge
  # is inclusive; the youngest outer edge is excluded (Quaternary).
  for (i in seq_len(nrow(stages))) {
    hit <- is.finite(x$midpoint_ma) & x$midpoint_ma > stages$younger_ma[i] &
      x$midpoint_ma <= stages$older_ma[i]
    x$stage_id[hit] <- i
  }
  x$exclusion[is.na(x$exclusion) & is.na(x$stage_id)] <- "midpoint_outside_neogene"
  x$uncertain_assignment <- NA
  assigned <- which(!is.na(x$stage_id))
  # Crossing, rather than merely touching, a boundary is uncertain.
  x$uncertain_assignment[assigned] <- x$MAX_AGE[assigned] > stages$older_ma[x$stage_id[assigned]] |
    x$MIN_AGE[assigned] < stages$younger_ma[x$stage_id[assigned]]
  candidate <- x[is.na(x$exclusion), ]
  if (!nrow(candidate)) stop("No eligible occurrences")
  points <- sf::st_as_sf(candidate, coords = c("LONG", "LAT"), remove = FALSE, crs = 4326)
  regions <- sf::st_read(cfg$region_file, quiet = TRUE)
  regions <- sf::st_transform(regions, 4326)
  hits <- sf::st_intersects(points, regions)
  inside <- lengths(hits) > 0
  x$exclusion[match(candidate$occurrence_id[!inside], x$occurrence_id)] <- "outside_original_region_polygons"
  candidate <- candidate[inside, ]
  if (!nrow(candidate)) stop("No occurrences inside original region polygons")
  # Stable version of the original geographic split, independent of row order.
  candidate$continent <- ifelse(candidate$LONG > -50, "Europe", "North America")
  if (any(!is.na(candidate$trophic_guild) &
      !candidate$trophic_guild %in% c("Herbivore", "Carnivore", "Omnivore",
        "Invertivore", "Piscivore")))
    stop("Unexpected NOW-derived trophic guild in fossil master")
  candidate$stage <- stages$stage[candidate$stage_id]
  candidate$duration_ma <- stages$duration_ma[candidate$stage_id]
  message("[1/5] Rebuilding original spatial resolutions (no temporal rounding)")
  grids <- list()
  for (res in names(cfg$spacing)) {
    b <- palaeoverse::bin_space(candidate, lng = "LONG", lat = "LAT",
                               spacing = cfg$spacing[[res]], return = TRUE)
    order_id <- match(candidate$occurrence_id, b$occdf$occurrence_id)
    if (anyNA(order_id) || anyDuplicated(b$occdf$occurrence_id)) stop("Spatial binning changed occurrence identity")
    for (field in c("cell_ID", "cell_centroid_lat", "cell_centroid_lng")) {
      if (!field %in% names(b$occdf)) stop("Unexpected palaeoverse output: ", field)
      candidate[[paste0(field, "_", res)]] <- b$occdf[[field]][order_id]
    }
    grids[[res]] <- list(grid = b$grid, grid_info = b$grid_info)
  }
  summary_rows <- lapply(c(NA_integer_, stages$stage_id), function(i) {
    z <- if (is.na(i)) candidate else candidate[candidate$stage_id == i, ]
    data.frame(stage = if (is.na(i)) "Overall" else stages$stage[i],
               n_occurrences = nrow(z), n_uncertain = sum(z$uncertain_assignment),
               percent_uncertain = if (nrow(z)) 100 * mean(z$uncertain_assignment) else NA_real_)
  })
  utils::write.csv(do.call(rbind, summary_rows), file.path(cfg$out, "assignment_uncertainty.csv"), row.names = FALSE)
  x$exclusion[is.na(x$exclusion)] <- "retained"
  utils::write.csv(as.data.frame(table(x$exclusion)), file.path(cfg$out, "occurrence_filter_counts.csv"), row.names = FALSE)
  utils::write.csv(x[, c("occurrence_id", "LIDNUM", "MIN_AGE", "MAX_AGE", "midpoint_ma", "stage_id", "uncertain_assignment", "exclusion")],
                   file.path(cfg$out, "occurrence_assignment_audit.csv"), row.names = FALSE)
  result <- list(occurrences = candidate, stages = stages, grids = grids, config = cfg)
  saveRDS(result, output)
  invisible(result)
}
