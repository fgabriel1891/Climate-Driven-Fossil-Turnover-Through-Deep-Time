# Figure 2 renderer. Run from the repository root after climate
# preparation. Sourcing this file defines functions only; see README.md.
source("figure_style/publication_style.R")

fig2_weighted_quantile <- function(x, w, probabilities = c(.1, .9)) {
  good <- is.finite(x) & is.finite(w) & w > 0
  if (!any(good)) return(rep(NA_real_, length(probabilities)))
  x <- x[good]; w <- w[good]
  order_id <- order(x)
  x <- x[order_id]; w <- w[order_id]
  x[vapply(probabilities, function(p) which(cumsum(w)/sum(w) >= p)[1], integer(1))]
}

fig2_prepare_data <- function(cfg) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install/use terra in the current R environment")
  if (!file.exists(file.path(cfg$out, "climate_catalog.rds"))) stop("Run prepare_bristol_climate(cfg) first")
  climate <- readRDS(file.path(cfg$out, "climate_catalog.rds"))
  stages <- climate$stages
  stopifnot(nrow(stages) == 8L, identical(stages$stage_id, 1:8),
            all(c("europe", "north_america") %in% names(climate$catalog)))
  metrics <- c("mean", "signed_change", "seasonal_change")
  all_rows <- list()
  maps <- list()
  for (region in c("europe", "north_america")) for (variable in c("temperature", "precipitation")) {
    path <- climate$catalog[[region]][[variable]]
    if (is.null(path) || !file.exists(path)) stop("Missing Bristol ", variable, " stack for ", region)
    bands_path <- sub("[.]tif$", "_bands.csv", path)
    if (!file.exists(bands_path)) stop("Missing band map: ", bands_path)
    stack <- terra::rast(path)
    bands <- utils::read.csv(bands_path)
    expected_units <- if (variable == "temperature") "degC" else "mm/day"
    if (nrow(bands) != terra::nlyr(stack) || any(bands$units != expected_units))
      stop("Unexpected band count or units: ", path)
    area <- terra::values(terra::cellSize(stack[[1]], unit = "km"), mat = FALSE)
    for (metric in metrics) {
      idx <- match(paste(1:8, metric), paste(bands$stage_id, bands$metric))
      if (anyNA(idx)) stop(variable, " metric missing for at least one stage: ", metric)
      layers <- stack[[bands$band[idx]]]
      values <- terra::values(layers, mat = TRUE)
      if (ncol(values) != 8L) stop("Expected eight stage layers")
      for (stage_id in 1:8) {
        x <- values[, stage_id]
        good <- is.finite(x) & is.finite(area) & area > 0
        if (!any(good)) stop("No valid raster cells for ", region, " / ", variable, " / ", metric, " / ", stages$stage[stage_id])
        q <- fig2_weighted_quantile(x, area)
        all_rows[[length(all_rows) + 1L]] <- data.frame(
          region = region, variable = variable, metric = metric, stage_id = stage_id, stage = stages$stage[stage_id],
          midpoint_ma = stages$stage_midpoint_ma[stage_id], duration_ma = stages$duration_ma[stage_id],
          mean = stats::weighted.mean(x[good], area[good]),
          p10 = q[1], p90 = q[2], n_cells = sum(good),
          units = expected_units)
      }
      # Match the existing Bristol publication-map definitions.
      if (metric == "mean") {
        duration <- stages$duration_ma
        summary_raster <- terra::app(layers, function(v) {
          if (any(!is.finite(v))) return(NA_real_)
          sum(v * duration)/sum(duration)
        })
      } else if (metric == "signed_change") {
        mean_idx <- match(paste(c(1, 8), "mean"), paste(bands$stage_id, bands$metric))
        summary_raster <- stack[[bands$band[mean_idx[2]]]] - stack[[bands$band[mean_idx[1]]]]
      } else {
        summary_raster <- terra::app(layers, function(v) {
          if (any(!is.finite(v))) return(NA_real_)
          mean(v)
        })
      }
      names(summary_raster) <- "value"
      maps[[paste(region, variable, metric, sep = "__")]] <- summary_raster
    }
  }
  list(trends = do.call(rbind, all_rows), maps = maps, stages = stages)
}

fig2_panel_title <- function(variable, metric, variable_names, metric_names) {
  if (metric == "mean") return(paste(metric_names[["mean"]], tolower(variable_names[[variable]])))
  if (metric == "signed_change")
    return(paste(variable_names[[variable]], tolower(metric_names[["change"]])))
  paste(variable_names[[variable]], tolower(metric_names[["seasonality"]]))
}

fig2_plot_trend <- function(data, variable, metric, stages, variable_names, metric_names) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install/use ggplot2")
  x <- data[data$variable == variable & data$metric == metric, ]
  x$region <- factor(x$region, levels = c("europe", "north_america"), labels = c("Europe", "North America"))
  units <- if (variable == "temperature") "deg C" else "mm/day"
  title <- fig2_panel_title(variable, metric, variable_names, metric_names)
  stage_labels <- c("Aq", "Bur", "Lan", "Ser", "Tor", "Mes", "Zan", "Pia")
  p <- ggplot2::ggplot(x, ggplot2::aes(midpoint_ma, mean, colour = region, group = region)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = p10, ymax = p90, fill = region),
      alpha = .16, colour = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = .8) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = figpub_colours$region, name = NULL) +
    ggplot2::scale_fill_manual(values = figpub_colours$region) +
    ggplot2::scale_x_reverse(breaks = stages$stage_midpoint_ma,
      labels = stage_labels, expand = ggplot2::expansion(mult = c(.03, .03)),
      guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::labs(title = title,
      x = NULL, y = units) + figpub_theme()
  if (metric == "signed_change") p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = .3, linetype = "dashed")
  p
}

fig2_plot_map <- function(maps, variable, metric, coastlines, land, land_polar,
                          land_only, variable_names, metric_names) {
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE)) stop("Use sf and terra")
  # Merge only the modelled European/North American windows; unmodelled cells
  # remain NA. The polar projection is for display, not climate inference.
  r <- terra::merge(maps[[paste("europe", variable, metric, sep = "__")]],
    maps[[paste("north_america", variable, metric, sep = "__")]])
  polar_crs <- "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m"
  if (land_only) r <- terra::mask(r, land, touches = FALSE)
  polar <- terra::project(r, polar_crs, method = "bilinear", mask = TRUE)
  if (land_only) polar <- terra::mask(polar, land_polar, touches = FALSE)
  xy <- as.data.frame(polar, xy = TRUE, na.rm = TRUE)
  if (!nrow(xy)) stop("Polar map has no finite data: ", variable, " / ", metric)
  names(xy)[3] <- "value"
  map_bbox <- sf::st_bbox(c(xmin = terra::xmin(polar), ymin = terra::ymin(polar),
    xmax = terra::xmax(polar), ymax = terra::ymax(polar)), crs = sf::st_crs(polar_crs))
  coast <- suppressWarnings(sf::st_crop(coastlines, map_bbox))
  metric_label <- fig2_panel_title(variable, metric, variable_names, metric_names)
  p <- ggplot2::ggplot(xy, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile(width = terra::res(polar)[1], height = terra::res(polar)[2]) +
    ggplot2::geom_sf(data = coast, inherit.aes = FALSE, colour = "grey20",
      linewidth = .25, fill = NA, show.legend = FALSE) +
    ggplot2::coord_sf(crs = sf::st_crs(polar_crs), default_crs = sf::st_crs(polar_crs),
      datum = NA, expand = FALSE) +
    ggplot2::labs(title = metric_label,
      x = NULL, y = NULL, fill = if (variable == "temperature") "deg C" else "mm/day") +
    figpub_map_theme() +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      barwidth = grid::unit(2.5, "cm"), barheight = grid::unit(.2, "cm"),
      title.position = "top"))
  if (metric == "signed_change") {
    limit <- max(abs(xy$value))
    if (limit == 0) limit <- 1
    p <- p + ggplot2::scale_fill_gradient2(low = figpub_colours$diverging[["low"]],
      mid = figpub_colours$diverging[["mid"]], high = figpub_colours$diverging[["high"]],
      midpoint = 0, limits = c(-limit, limit), na.value = "transparent")
  } else if (metric == "mean") {
    p <- p + ggplot2::scale_fill_viridis_c(option = "C", na.value = "transparent")
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(option = "C", na.value = "transparent")
  }
  p
}

render_fig2_climate <- function(cfg = bristol_config(),
                               out_dir = "fig2_temperature/output/climate_six_panels", overwrite = FALSE,
                               trends_title = "Neogene climate across eight stages",
                               maps_title = "Spatial climate summaries across the Neogene",
                               variable_names = c(temperature = "Temperature", precipitation = "Precipitation"),
                               metric_names = c(mean = "Mean", change = "Change", seasonality = "Seasonality"),
                               annotations = c(trends = TRUE, maps = TRUE), land_only = TRUE) {
  for (package in c("terra", "sf", "ggplot2", "patchwork", "rnaturalearth")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package: ", package)
  }
  if (!all(c("temperature", "precipitation") %in% names(variable_names)) ||
      !all(c("mean", "change", "seasonality") %in% names(metric_names)))
    stop("variable_names and metric_names require their documented named entries")
  if (!is.logical(annotations) || length(annotations) < 1L || length(annotations) > 2L || anyNA(annotations))
    stop("annotations must be TRUE/FALSE or c(trends = TRUE/FALSE, maps = TRUE/FALSE)")
  if (length(annotations) == 1L) annotations <- rep(annotations, 2L)
  if (is.null(names(annotations))) names(annotations) <- c("trends", "maps")
  if (!all(c("trends", "maps") %in% names(annotations)))
    stop("annotations must name trends and maps")
  if (!is.logical(land_only) || length(land_only) != 1L || is.na(land_only))
    stop("land_only must be TRUE or FALSE")
  prepared <- fig2_prepare_data(cfg)
  metrics <- c("mean", "signed_change", "seasonal_change")
  panel_specs <- expand.grid(metric = metrics, variable = c("temperature", "precipitation"),
    stringsAsFactors = FALSE)
  upper <- lapply(seq_len(nrow(panel_specs)), function(i)
    fig2_plot_trend(prepared$trends, panel_specs$variable[i], panel_specs$metric[i],
      prepared$stages, variable_names, metric_names))
  coastlines <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
  coastlines <- sf::st_transform(coastlines, "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m")
  land <- land_polar <- NULL
  if (land_only) {
    land <- terra::vect(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"))
    land <- terra::crop(land, terra::ext(-180, 180, 0, 90))
    land_polar <- terra::project(land, "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m")
  }
  lower <- lapply(seq_len(nrow(panel_specs)), function(i)
    fig2_plot_map(prepared$maps, panel_specs$variable[i], panel_specs$metric[i],
      coastlines, land, land_polar, land_only, variable_names, metric_names))
  names(upper) <- names(lower) <- paste(panel_specs$variable, panel_specs$metric, sep = "_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  planned <- file.path(out_dir, c("fig02_climate_trends.png", "fig02_climate_trends.pdf",
    "fig02_climate_maps.png", "fig02_climate_maps.pdf", "fig02_trend_data.csv"))
  if (!overwrite && any(file.exists(planned))) stop("Figure outputs exist; choose another out_dir or overwrite=TRUE")
  trends <- patchwork::wrap_plots(upper, ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  maps <- patchwork::wrap_plots(lower, ncol = 3)
  if (annotations[["trends"]]) trends <- trends +
    figpub_annotation(trends_title,
      "Change: incoming stage difference. Seasonality: SD of 12 monthly changes. Ribbons: spatial 10th-90th percentiles")
  if (annotations[["maps"]]) maps <- maps +
    figpub_annotation(maps_title,
      "Change: Piacenzian minus Aquitanian. Seasonality: mean incoming monthly-pattern change. Present-day display mask")
  figpub_save(out_dir, "fig02_climate_trends", trends, height = 6.3, overwrite = overwrite)
  figpub_save(out_dir, "fig02_climate_maps", maps, height = 7.3, overwrite = overwrite)
  utils::write.csv(prepared$trends, file.path(out_dir, "fig02_trend_data.csv"), row.names = FALSE)
  invisible(list(trends = trends, maps = maps, trend_data = prepared$trends))
}

# Keep the prior entry point working for users who already sourced this script.
render_fig2_temperature <- render_fig2_climate
