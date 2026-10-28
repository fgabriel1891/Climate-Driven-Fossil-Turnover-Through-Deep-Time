# Fork of src/gdm_functions.R: original cell x time design, now stage bins.
extract_bristol_hex <- function(r, cell_ids) {
  polygons <- sf::st_sf(cell_id = cell_ids, geometry = h3jsr::cell_to_polygon(cell_ids, simple = TRUE))
  sf::st_crs(polygons) <- 4326
  v <- terra::vect(polygons)
  if (!isTRUE(terra::same.crs(r, v))) v <- terra::project(v, terra::crs(r))
  z <- as.data.frame(terra::extract(r, v, fun = mean, exact = TRUE, na.rm = TRUE)[, -1, drop = FALSE])
  names(z) <- names(r)
  z
}

prepare_bristol_assemblages <- function(cfg = bristol_config()) {
  bristol_packages(c("dplyr", "terra", "sf", "h3jsr"))
  bristol_require_outputs(cfg, c("stage_fossils.rds", "climate_catalog.rds"), "prepare_bristol_fossils(cfg) and prepare_bristol_climate(cfg)")
  output <- bristol_output(cfg, "assemblage_catalog.rds")
  if (cfg$coordinate_mode != "present_day") stop("Only original fixed-coordinate convention is implemented")
  fossils <- readRDS(file.path(cfg$out, "stage_fossils.rds"))
  climate <- readRDS(file.path(cfg$out, "climate_catalog.rds"))
  if (!identical(fossils$stages, climate$stages)) stop("Fossil/climate stages differ")
  s <- fossils$stages
  catalog <- coverage <- list()
  message("[3/5] Building stage assemblages and extracting six climate predictors")
  for (resolution in names(cfg$spacing)) {
    for (continent in c("Europe", "North America")) {
      region <- if (continent == "Europe") "europe" else "north_america"
      for (guild in c("all", "herb", "carn")) {
        key <- paste(resolution, region, guild, sep = "__")
        x <- fossils$occurrences[fossils$occurrences$continent == continent, ]
        if (guild != "all") {
          wanted <- if (guild == "herb") "Herbivore" else "Carnivore"
          x <- x[!is.na(x$trophic_guild) & x$trophic_guild == wanted, ]
        }
        x$cell_id <- as.character(x[[paste0("cell_ID_", resolution)]])
        x <- x[!is.na(x$cell_id), ]
        coverage[[key]] <- data.frame(model_group = key, assemblages = 0, eligible = 0, climate_complete = 0)
        if (!nrow(x)) next
        sites <- x |>
          dplyr::group_by(cell_id, stage_id) |>
          dplyr::summarise(genera = list(sort(unique(GENUS))), richness = dplyr::n_distinct(GENUS),
                           n_records = dplyr::n(), n_locations = dplyr::n_distinct(LIDNUM, na.rm = TRUE),
                           lng = dplyr::first(.data[[paste0("cell_centroid_lng_", resolution)]]),
                           lat = dplyr::first(.data[[paste0("cell_centroid_lat_", resolution)]]), .groups = "drop")
        coverage[[key]]$assemblages <- nrow(sites)
        sites <- as.data.frame(sites[sites$richness >= cfg$min_genera, ])
        coverage[[key]]$eligible <- nrow(sites)
        if (!nrow(sites)) next
        sites$site <- paste0(sites$cell_id, "_stage", sites$stage_id)
        sites$stage <- s$stage[sites$stage_id]
        sites$duration_ma <- s$duration_ma[sites$stage_id]
        sites$time_ma <- s$stage_midpoint_ma[sites$stage_id]
        sites$sampling <- log1p(sites$n_locations)
        for (variable in c("temperature", "precipitation")) {
          path <- climate$catalog[[region]][[variable]]
          r <- terra::rast(path)
          bt <- utils::read.csv(sub("[.]tif$", "_bands.csv", path))
          for (metric in c("mean", "signed_change", "seasonal_change")) sites[[paste(variable, metric, sep = "_")]] <- NA_real_
          for (stage in sort(unique(sites$stage_id))) {
            rows <- which(sites$stage_id == stage)
            idx <- match(paste(stage, c("mean", "signed_change", "seasonal_change")), paste(bt$stage_id, bt$metric))
            if (anyNA(idx)) stop("Missing stage predictors")
            part <- r[[bt$band[idx]]]
            names(part) <- paste(variable, bt$metric[idx], sep = "_")
            z <- extract_bristol_hex(part, sites$cell_id[rows])
            for (field in names(z)) sites[rows, field] <- z[[field]]
          }
        }
        climate_cols <- grep("^(temperature|precipitation)_", names(sites), value = TRUE)
        sites$climate_complete <- apply(as.matrix(sites[, climate_cols]), 1, function(z) all(is.finite(z)))
        # Euclidean geography in km, using continent-specific azimuthal equidistant projections.
        crs <- if (continent == "Europe") "+proj=aeqd +lat_0=52 +lon_0=15 +datum=WGS84 +units=m" else "+proj=aeqd +lat_0=45 +lon_0=-100 +datum=WGS84 +units=m"
        points <- sf::st_as_sf(sites[, c("lng", "lat")], coords = c("lng", "lat"), crs = 4326)
        xy <- sf::st_coordinates(sf::st_transform(points, crs))/1000
        sites$x_km <- xy[, 1]
        sites$y_km <- xy[, 2]
        path <- file.path(cfg$out, paste0("assemblages_", key, ".rds"))
        saveRDS(list(sites = sites, key = key, continent = continent, guild = guild, resolution = resolution,
                     geographic_crs = crs, stages = s, config = cfg), path)
        utils::write.csv(sites[, setdiff(names(sites), "genera")], sub("[.]rds$", ".csv", path), row.names = FALSE)
        catalog[[key]] <- path
        coverage[[key]]$climate_complete <- sum(sites$climate_complete)
      }
    }
  }
  utils::write.csv(do.call(rbind, coverage), file.path(cfg$out, "assemblage_coverage.csv"), row.names = FALSE)
  saveRDS(catalog, output)
  invisible(catalog)
}
