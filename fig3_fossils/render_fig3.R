# Manuscript Fig. 3 and temporal-turnover supplement from saved Bristol outputs.
# Sourcing this file defines functions only; render_fig3() writes outputs.
source("figure_style/publication_style.R")

fig3_theme <- function() {
  figpub_theme()
}

fig3_required <- function(path) {
  if (!file.exists(path)) stop("Missing saved Bristol output: ", path, call. = FALSE)
  path
}

fig3_hex_data <- function(cfg, region) {
  key <- paste("h3", region, "all", sep = "__")
  obj <- readRDS(fig3_required(file.path(cfg$out, paste0("assemblages_", key, ".rds"))))
  sites <- obj$sites
  genera <- split(sites$genera, sites$cell_id)
  ids <- names(genera)
  richness <- vapply(genera, function(x) length(unique(unlist(x, use.names = FALSE))), integer(1))
  polygons <- sf::st_sf(cell_id = ids, genus_richness = richness,
    geometry = h3jsr::cell_to_polygon(ids, simple = TRUE))
  sf::st_crs(polygons) <- 4326
  polygons
}

fig3_map <- function(hexes, coastlines, region, richness_limit, palette = "C") {
  centre <- if (region == "Europe") c(lon = 15, lat = 52) else c(lon = -100, lat = 45)
  projection <- sf::st_crs(sprintf(
    "+proj=laea +lat_0=%s +lon_0=%s +datum=WGS84 +units=m +no_defs",
    centre[["lat"]], centre[["lon"]]))
  hexes <- sf::st_transform(hexes, projection)
  box <- sf::st_bbox(hexes)
  half_span <- max(box[["xmax"]] - box[["xmin"]],
    box[["ymax"]] - box[["ymin"]]) * 0.56
  xmid <- (box[["xmin"]] + box[["xmax"]]) / 2
  ymid <- (box[["ymin"]] + box[["ymax"]]) / 2
  xlim <- xmid + c(-half_span, half_span)
  ylim <- ymid + c(-half_span, half_span)
  coast <- suppressWarnings(sf::st_crop(sf::st_transform(coastlines, projection),
    sf::st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
      crs = projection)))
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = coast, colour = "grey70", linewidth = .3) +
    ggplot2::geom_sf(data = hexes, ggplot2::aes(fill = genus_richness),
      colour = "grey35", linewidth = .06) +
    ggplot2::coord_sf(crs = projection, xlim = xlim, ylim = ylim,
      expand = FALSE, datum = NA) +
    ggplot2::scale_fill_viridis_c(option = palette, limits = c(0, richness_limit),
      name = "Genus richness") +
    ggplot2::labs(title = region, x = NULL, y = NULL) + figpub_map_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 16))
}

fig3_matrix_data <- function(cfg) {
  x <- utils::read.csv(fig3_required(file.path(cfg$tables, "turnover_stage_matrix.csv")))
  x <- x[x$resolution == "h3" & x$guild == "all", ]
  if (!all(c("Europe", "North America") %in% x$continent) || nrow(x) != 128L)
    stop("Expected complete 8 by 8 H3 all-mammal turnover matrices for both continents")
  x
}

fig3_matrix_plot <- function(data, stages, region,
    colours = c("#ff9f1c", "#ffbf69", "#ffff74", "#cbf3f0", "#2ec4b6")) {
  x <- data[data$continent == region, ]
  stage_labels <- c("Aq", "Bur", "Lan", "Ser", "Tor", "Mes", "Zan", "Pia")
  x$stage1 <- factor(stage_labels[x$stage1], levels = stage_labels)
  x$stage2 <- factor(stage_labels[x$stage2], levels = rev(stage_labels))
  ggplot2::ggplot(x, ggplot2::aes(stage1, stage2, fill = mean_turnover)) +
    ggplot2::geom_tile(colour = "white", linewidth = .35) +
    ggplot2::scale_fill_gradientn(colours = colours, limits = c(0, 1), na.value = "grey92",
      name = "Mean Simpson turnover") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = region, x = NULL, y = NULL) + fig3_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

fig3_turnover_data <- function(cfg, region, guild) {
  key <- paste("h3", region, guild, sep = "__")
  obj <- readRDS(fig3_required(file.path(cfg$out, paste0("turnover_", key, ".rds"))))
  d <- as.matrix(obj$distance)
  sites <- obj$sites
  if (nrow(d) != nrow(sites) || ncol(d) != nrow(sites))
    stop("Turnover matrix and site metadata differ: ", key)
  pair <- which(upper.tri(d), arr.ind = TRUE)
  gap <- round(abs(sites$time_ma[pair[, 1]] - sites$time_ma[pair[, 2]]), 6)
  turnover <- d[pair]
  keep <- is.finite(gap) & gap > 0 & is.finite(turnover)
  if (!any(keep)) stop("No cross-stage turnover pairs: ", key)
  pairs <- data.frame(temporal_separation_ma = gap[keep], turnover = turnover[keep])
  out <- dplyr::summarise(dplyr::group_by(pairs, temporal_separation_ma),
    mean_turnover = mean(turnover), n_pairs = dplyr::n(), .groups = "drop")
  out$continent <- obj$continent
  out$guild <- guild
  as.data.frame(out)
}

fig3_turnover_plot <- function(data, region) {
  x <- data[data$continent == region, ]
  x$guild <- factor(x$guild, levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  ggplot2::ggplot(x, ggplot2::aes(temporal_separation_ma, mean_turnover,
    colour = guild, group = guild)) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = figpub_colours$guild, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = region, x = "Stage-midpoint separation (Myr)",
      y = "Mean Simpson turnover") + fig3_theme()
}

render_fig3 <- function(cfg = bristol_config(), out_dir = "fig3_fossils/output",
                        overwrite = FALSE,
                        main_title = "Neogene mammal richness and compositional turnover",
                        supplement_title = "Temporal separation and mammal turnover",
                        map_palette = "C",
                        matrix_colours = rev(c("#F7FBFF", "#C6DBEF", "#6BAED6",
                          "#2171B5", "#08306B"))) {
  for (package in c("ggplot2", "patchwork", "dplyr", "sf", "h3jsr", "rnaturalearth"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package: ", package)
  regions <- c("europe", "north_america")
  hexes <- lapply(regions, function(region) fig3_hex_data(cfg, region))
  names(hexes) <- regions
  coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
  max_richness <- max(vapply(hexes, function(x) max(x$genus_richness), numeric(1)))
  maps <- list(fig3_map(hexes$europe, coast, "Europe", max_richness,
      palette = map_palette),
    fig3_map(hexes$north_america, coast, "North America", max_richness,
      palette = map_palette))
  matrix_data <- fig3_matrix_data(cfg)
  stages <- read_bristol_stages(cfg)
  matrices <- list(fig3_matrix_plot(matrix_data, stages, "Europe",
      colours = matrix_colours),
    fig3_matrix_plot(matrix_data, stages, "North America",
      colours = matrix_colours))
  turnover <- do.call(rbind, lapply(regions, function(region)
    do.call(rbind, lapply(c("all", "herb", "carn"), function(guild)
      fig3_turnover_data(cfg, region, guild)))))
  decay <- list(fig3_turnover_plot(turnover, "Europe"),
    fig3_turnover_plot(turnover, "North America"))
  figure <- patchwork::wrap_plots(c(maps, matrices), ncol = 2,
      widths = c(1, 1), heights = c(1, 1), guides = "collect") +
    figpub_annotation(main_title,
      "H3 genus richness (A-B); mean all-mammal Simpson turnover between stage pairs (C-D)") &
    ggplot2::theme(legend.position = "bottom")
  supplement <- patchwork::wrap_plots(decay, ncol = 2, guides = "collect") +
    figpub_annotation(supplement_title,
      "Points summarize distinct stage-midpoint distances; lines connect observed means") &
    ggplot2::theme(legend.position = "bottom")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  planned <- file.path(out_dir, c("fig03_fossil_richness_turnover_matrix.png",
    "fig03_fossil_richness_turnover_matrix.pdf",
    "figS_temporal_turnover.png", "figS_temporal_turnover.pdf",
    "figS_turnover_by_midpoint_distance.csv"))
  if (!overwrite && any(file.exists(planned)))
    stop("Fig. 3 output already exists; set overwrite=TRUE or choose another out_dir")
  figpub_save(out_dir, "fig03_fossil_richness_turnover_matrix", figure,
    height = 7.8, overwrite = overwrite)
  figpub_save(out_dir, "figS_temporal_turnover", supplement,
    height = 4.1, overwrite = overwrite)
  utils::write.csv(turnover, file.path(out_dir, "figS_turnover_by_midpoint_distance.csv"),
    row.names = FALSE)
  invisible(list(figure = figure, supplement = supplement,
    turnover_matrix = matrix_data, temporal_turnover = turnover))
}
