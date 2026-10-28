# Publication counterparts of the parent project's illustration scripts.
# All functions use saved Bristol products. Source-time execution is disabled.
pub_theme <- function() {
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey95"),
    strip.text = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"))
}

save_publication <- function(plot, stem, cfg, width = 12, height = 8) {
  dir.create(cfg$publication_figs, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(cfg$publication_figs, paste0(stem, ".", cfg$figure_formats))
  if (any(file.exists(paths)) && !cfg$overwrite) stop("Figure already exists: ", stem)
  for (path in paths) ggplot2::ggsave(path, plot = plot, width = width, height = height,
                                     units = "in", dpi = cfg$figure_dpi, bg = "white", limitsize = FALSE)
  invisible(paths)
}

publication_fossils <- function(cfg, tables) {
  fossils <- readRDS(file.path(cfg$out, "stage_fossils.rds"))
  x <- fossils$occurrences
  stages <- fossils$stages
  x$stage <- factor(stages$stage[x$stage_id], levels = stages$stage)
  locations <- unique(x[, c("continent", "LONG", "LAT", "stage")])
  aggregate <- locations |> dplyr::group_by(continent, LONG, LAT) |>
    dplyr::summarise(n_stages = dplyr::n_distinct(stage), .groups = "drop")
  map <- ggplot2::ggplot(aggregate, ggplot2::aes(LONG, LAT)) +
    ggplot2::geom_point(ggplot2::aes(colour = n_stages), alpha = .7, size = 1.1) +
    ggplot2::facet_wrap(~continent, scales = "free") + ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(x = "Longitude", y = "Latitude", colour = "Stages recorded", title = "Neogene fossil distribution") + pub_theme()
  save_publication(map, "fossil_distribution_aggregated", cfg, 12, 5)
  by_stage <- ggplot2::ggplot(locations, ggplot2::aes(LONG, LAT)) +
    ggplot2::geom_point(size = .65, alpha = .65, colour = "#355C7D") +
    ggplot2::facet_wrap(~continent + stage, scales = "free", ncol = 4, drop = FALSE) +
    ggplot2::labs(x = "Longitude", y = "Latitude", title = "Fossil distribution by stage") + pub_theme()
  save_publication(by_stage, "fossil_distribution_by_time", cfg, 14, 11)
  rich <- tables$genus_species_richness_by_stage
  rich$stage <- factor(rich$stage, levels = stages$stage)
  long <- rbind(data.frame(rich, rank = "Genera", richness = rich$n_genera),
                data.frame(rich, rank = "Identified species", richness = rich$n_species))
  p <- ggplot2::ggplot(long, ggplot2::aes(stage, richness, colour = guild, group = guild)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_grid(rank ~ continent, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(all = "#444444", herb = "#009E73", carn = "#D55E00")) +
    ggplot2::labs(x = NULL, y = "Observed richness", title = "Taxonomic richness across Neogene stages") +
    pub_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_publication(p, "genus_species_richness_by_time", cfg)
  p <- ggplot2::ggplot(rich, ggplot2::aes(n_localities, n_genera, colour = guild)) +
    ggplot2::geom_point(ggplot2::aes(shape = stage), size = 2) + ggplot2::facet_wrap(~continent) +
    ggplot2::scale_shape_manual(values = c(0:6, 8)) +
    ggplot2::labs(x = "Recorded localities", y = "Observed genera", title = "Richness and sampling by stage") + pub_theme()
  save_publication(p, "genus_species_richness_sampling", cfg)
  catalog <- readRDS(file.path(cfg$out, "assemblage_catalog.rds"))
  all_sites <- dplyr::bind_rows(lapply(catalog, function(path) {
    o <- readRDS(path)
    if (o$guild != "all") return(NULL)
    data.frame(continent = o$continent, resolution = o$resolution, stage = o$sites$stage,
               richness = o$sites$richness)
  }))
  p <- ggplot2::ggplot(all_sites, ggplot2::aes(resolution, richness)) + ggplot2::geom_boxplot(outlier.size = .4) +
    ggplot2::facet_wrap(~continent, scales = "free_y") +
    ggplot2::labs(x = "Spatial resolution", y = "Genera per eligible cell-stage assemblage") + pub_theme()
  save_publication(p, "genus_richness_h0_h3", cfg)
  invisible(map)
}

publication_climate <- function(cfg) {
  object <- readRDS(file.path(cfg$out, "climate_catalog.rds"))
  summaries <- maps <- list()
  s <- object$stages
  summary_rows <- list()
  for (variable in c("temperature", "precipitation")) for (metric in c("mean", "signed_change", "seasonal_change")) {
    region_maps <- list()
    for (region in names(object$catalog)) {
      r <- terra::rast(object$catalog[[region]][[variable]])
      bt <- utils::read.csv(sub("[.]tif$", "_bands.csv", object$catalog[[region]][[variable]]))
      idx <- bt$band[match(paste(1:8, metric), paste(bt$stage_id, bt$metric))]
      if (anyNA(idx)) stop("Climate band map incomplete")
      layers <- r[[idx]]
      w <- terra::values(terra::cellSize(layers[[1]], unit = "km"), mat = FALSE)
      values <- terra::values(layers, mat = TRUE)
      for (i in 1:8) {
        ok <- is.finite(values[, i]) & is.finite(w)
        if (!any(ok)) stop("No finite climate values")
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(region, variable, metric, stage_id = i,
          stage = s$stage[i], age_ma = s$stage_midpoint_ma[i], units = if (variable == "temperature") "degC" else "mm/day",
          spatial_mean = stats::weighted.mean(values[ok, i], w[ok]), minimum = min(values[ok, i]), maximum = max(values[ok, i]),
          n_cells = sum(ok))
      }
      if (metric == "mean") {
        summary_map <- terra::app(layers, function(v) if (any(!is.finite(v))) NA_real_ else sum(v * s$duration_ma)/sum(s$duration_ma))
      } else if (metric == "signed_change") {
        means <- r[[bt$band[bt$metric == "mean"]]]
        summary_map <- means[[8]] - means[[1]]
      } else {
        summary_map <- terra::app(layers, function(v) if (any(!is.finite(v))) NA_real_ else mean(v))
      }
      names(summary_map) <- "value"
      dat <- as.data.frame(summary_map, xy = TRUE, na.rm = TRUE)
      p <- ggplot2::ggplot(dat, ggplot2::aes(x, y, fill = value)) + ggplot2::geom_raster() +
        ggplot2::coord_equal() + ggplot2::labs(title = gsub("_", " ", region), x = "Longitude", y = "Latitude",
          fill = if (variable == "temperature") "degC" else "mm/day") + pub_theme()
      if (metric == "signed_change") {
        limit <- max(abs(dat$value))
        if (limit == 0) limit <- 1
        p <- p + ggplot2::scale_fill_gradient2(low = "#272757", mid = "grey95", high = "#D55E00", limits = c(-limit, limit))
      } else p <- p + ggplot2::scale_fill_viridis_c()
      region_maps[[region]] <- p
      names(layers) <- s$stage
      dd <- as.data.frame(layers, xy = TRUE, na.rm = FALSE)
      stage_map <- dplyr::bind_rows(lapply(1:8, function(i) data.frame(x = dd$x, y = dd$y, value = dd[[i + 2]], stage = factor(s$stage[i], levels = s$stage))))
      p_stage <- ggplot2::ggplot(stage_map, ggplot2::aes(x, y, fill = value)) + ggplot2::geom_raster() +
        ggplot2::facet_wrap(~stage, ncol = 4) + ggplot2::coord_equal() + ggplot2::scale_fill_viridis_c(na.value = "transparent") +
        ggplot2::labs(title = paste(region, variable, metric), fill = if (variable == "temperature") "degC" else "mm/day") + pub_theme()
      save_publication(p_stage, paste(region, variable, metric, "by_stage", sep = "_"), cfg, 14, 8)
    }
    title <- switch(metric, mean = "Duration-weighted mean of eight stage climatologies",
                    signed_change = "Piacenzian minus Aquitanian mean climate",
                    seasonal_change = "Mean incoming seasonal-pattern change (includes Chattian baseline)")
    panel <- patchwork::wrap_plots(region_maps, nrow = 1) + patchwork::plot_annotation(title = paste(variable, title, sep = ": "))
    stem <- if (variable == "temperature") switch(metric, mean = "neogene_average_summary", signed_change = "neogene_cooling_summary", seasonal_change = "neogene_seasonal_change_summary") else paste0("precipitation_", metric, "_summary")
    save_publication(panel, stem, cfg, 12, 6)
    maps[[paste(variable, metric)]] <- panel
  }
  summary <- do.call(rbind, summary_rows)
  manuscript_write_table(summary, "regional_climate_by_stage", cfg)
  trend <- ggplot2::ggplot(summary, ggplot2::aes(age_ma, spatial_mean, colour = region)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::scale_x_reverse() +
    ggplot2::facet_wrap(~variable + metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(x = "Stage midpoint (Ma; bin descriptor, not simulation age)", y = "Area-weighted regional mean") + pub_theme()
  save_publication(trend, "neogene_climate_trends", cfg, 14, 7)
  composite <- patchwork::wrap_plots(maps, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
  save_publication(composite, "climate_multipanel_neogene", cfg, 16, 18)
  invisible(summary)
}

publication_turnover <- function(cfg, tables, fossil_map = NULL) {
  catalog <- readRDS(file.path(cfg$out, "turnover_catalog.rds"))
  stages <- read_bristol_stages(cfg)
  ordinations <- panels <- contour_panels <- list()
  set.seed(cfg$random_seed)
  for (region in c("europe", "north_america")) {
    key <- paste(cfg$figure_resolution, region, "all", sep = "__")
    if (is.null(catalog[[key]])) next
    obj <- readRDS(catalog[[key]])
    if (nrow(obj$sites) < 4L) next
    ord <- vegan::metaMDS(obj$distance, k = 2, trymax = 40, trace = FALSE)
    points <- as.data.frame(vegan::scores(ord, display = "sites"))
    sites <- obj$sites[match(rownames(points), obj$sites$site), ]
    points$stage <- factor(sites$stage, levels = stages$stage)
    points$temperature <- sites$temperature_mean
    points$site <- sites$site
    ordinations[[region]] <- list(ordination = ord, points = points)
    p <- ggplot2::ggplot(points, ggplot2::aes(NMDS1, NMDS2, colour = stage)) + ggplot2::geom_point(size = 1.6, alpha = .65) +
      ggplot2::scale_colour_viridis_d(drop = FALSE) + ggplot2::coord_equal() + pub_theme() +
      ggplot2::labs(title = obj$continent, subtitle = sprintf("NMDS stress: %.3f", ord$stress), colour = "Stage")
    panels[[region]] <- p
    ok <- is.finite(points$temperature)
    if (sum(ok) > 5 && length(unique(points$temperature[ok])) > 3) {
      surface <- vegan::ordisurf(as.matrix(points[ok, c("NMDS1", "NMDS2")]), points$temperature[ok], knots = 4, plot = FALSE)
      grid <- expand.grid(NMDS1 = surface$grid$x, NMDS2 = surface$grid$y)
      grid$temperature <- as.vector(surface$grid$z)
      contour_panels[[region]] <- p + ggplot2::geom_contour(data = grid,
        ggplot2::aes(NMDS1, NMDS2, z = temperature), inherit.aes = FALSE, colour = "grey20", bins = 7, linewidth = .35)
      ordinations[[region]]$surface <- surface
    }
  }
  saveRDS(ordinations, file.path(cfg$out, "publication_nmds.rds"))
  if (length(panels)) {
    p <- patchwork::wrap_plots(panels, nrow = 1, guides = "collect")
    save_publication(p, "turnover_mds_all_only", cfg, 13, 6)
    if (!is.null(fossil_map)) save_publication(patchwork::wrap_plots(fossil_map, p, ncol = 1), "fossil_distribution_with_nmds", cfg, 13, 11)
  }
  if (length(contour_panels)) save_publication(patchwork::wrap_plots(contour_panels, nrow = 1, guides = "collect"), "turnover_mds_temperature_contours", cfg, 13, 6)
  d <- tables$turnover_same_stage_spatial
  d <- d[d$resolution == cfg$figure_resolution, ]
  d$stage <- factor(d$stage, levels = stages$stage)
  p <- ggplot2::ggplot(d, ggplot2::aes(stage, mean_turnover, colour = guild, group = guild)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75), width = .15) +
    ggplot2::facet_wrap(~continent) + ggplot2::coord_cartesian(ylim = c(0, 1)) + pub_theme() +
    ggplot2::labs(x = NULL, y = "Same-stage spatial beta.sim (mean; bars = IQR)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_publication(p, paste0("fossil_turnover_summary_", cfg$figure_resolution), cfg)
  d <- tables$turnover_same_cell_by_stage_lag
  d <- d[d$resolution == cfg$figure_resolution, ]
  p <- ggplot2::ggplot(d, ggplot2::aes(stage_lag, mean_turnover, colour = guild)) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~continent) + ggplot2::labs(x = "Number of stage boundaries (not Myr)", y = "Mean same-cell beta.sim") + pub_theme()
  save_publication(p, paste0("fossil_turnover_temporal_lag_", cfg$figure_resolution), cfg)
  d <- tables$turnover_stage_matrix
  d <- d[d$resolution == cfg$figure_resolution & d$guild == "all", ]
  d$stage1 <- factor(stages$stage[d$stage1], levels = stages$stage)
  d$stage2 <- factor(stages$stage[d$stage2], levels = rev(stages$stage))
  p <- ggplot2::ggplot(d, ggplot2::aes(stage1, stage2, fill = mean_turnover)) + ggplot2::geom_tile() +
    ggplot2::facet_wrap(~continent) + ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey95") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Mean beta.sim", title = "Stage-aggregated turnover matrix") + pub_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_publication(p, paste0("simpson_dissimilarity_matrix_", cfg$figure_resolution, "_all"), cfg)
}

publication_gdm <- function(cfg, tables) {
  catalog <- readRDS(file.path(cfg$out, "gdm_catalog.rds"))
  rows <- list()
  for (id in names(catalog)) {
    obj <- readRDS(catalog[[id]])
    sp <- gdm::isplineExtract(obj$model)
    labels <- colnames(sp$x)
    if (is.null(labels)) labels <- obj$model$predictors
    key <- strsplit(obj$group, "__", fixed = TRUE)[[1]]
    for (j in seq_len(ncol(sp$x))) rows[[length(rows) + 1L]] <- data.frame(
      model = id, resolution = key[1], continent = if (key[2] == "europe") "Europe" else "North America",
      guild = key[3], variant = obj$variant, geo = obj$geo,
      predictor = labels[j], x = sp$x[, j], y = sp$y[, j])
  }
  if (!length(rows)) stop("No fitted models for publication splines")
  curves <- do.call(rbind, rows)
  manuscript_write_table(curves, "all_gdm_spline_coordinates", cfg)
  predictor_labels <- c(
    temperature_mean = "Mean temperature (°C)",
    precipitation_mean = "Mean precipitation (mm/day)",
    temperature_signed_change = "Temperature change (°C)",
    precipitation_signed_change = "Precipitation change (mm/day)",
    temperature_seasonal_change = "Temperature seasonal-pattern change (°C)",
    precipitation_seasonal_change = "Precipitation seasonal-pattern change (mm/day)"
  )
  for (guild in unique(curves$guild)) for (geo in cfg$geo_models) {
    d <- curves[curves$guild == guild & curves$geo == geo & curves$variant == "combined" & curves$predictor %in% bristol_climate_names(), ]
    if (!nrow(d)) next
    # Draw coarser resolutions first so the H3 reference curve remains visible on top.
    d$resolution <- factor(d$resolution, levels = c("h0", "h1", "h2", "h3"))
    d <- d[order(d$resolution), ]
    predictor_order <- names(predictor_labels)
    d$predictor_label <- factor(unname(predictor_labels[d$predictor]),
                                levels = unname(predictor_labels[predictor_order]))
    if (anyNA(d$predictor_label)) stop("Missing display label for a climate predictor")
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y, alpha = resolution, linewidth = resolution, group = model)) +
      ggplot2::geom_line(colour = "black") +
      ggplot2::scale_alpha_manual(name = "Spatial resolution", values = c(h0 = .16, h1 = .28, h2 = .45, h3 = 1),
                                  breaks = c("h3", "h2", "h1", "h0"), labels = c("H3", "H2", "H1", "H0")) +
      ggplot2::scale_linewidth_manual(values = c(h0 = .35, h1 = .45, h2 = .6, h3 = 1.05), guide = "none") +
      ggplot2::facet_wrap(~continent + predictor_label, scales = "free", ncol = 3) +
      ggplot2::labs(x = "Climate predictor", y = "Cumulative GDM contribution",
                    title = paste("Bristol GDM", guild, "geography", geo)) + pub_theme()
    save_publication(p, paste0("all_splines_", guild, "_geo", if (geo) "T" else "F"), cfg, 14, 13)
  }
  status <- tables$Table_S2_GDM_summary
  if (is.null(status)) return(invisible(curves))
  status <- status[status$ok & is.finite(status$explained), ]
  p <- ggplot2::ggplot(status, ggplot2::aes(resolution, explained, colour = guild, group = interaction(guild, variant))) +
    ggplot2::geom_line(ggplot2::aes(linetype = variant)) + ggplot2::geom_point(ggplot2::aes(shape = variant), size = 2) +
    ggplot2::facet_grid(geo ~ continent) + ggplot2::labs(x = "Spatial resolution", y = "In-sample deviance explained (%)") + pub_theme()
  save_publication(p, "dev_explained_geo_on_off", cfg)
  p <- ggplot2::ggplot(status, ggplot2::aes(sites, explained, colour = guild, shape = variant)) +
    ggplot2::geom_point(size = 2) + ggplot2::facet_grid(geo ~ continent) +
    ggplot2::labs(x = "Fitted cell-stage assemblages", y = "In-sample deviance explained (%)") + pub_theme()
  save_publication(p, "deviance_relationship", cfg)
  invisible(curves)
}

publication_figures <- function(cfg = bristol_config()) {
  bristol_packages(c("ggplot2", "patchwork", "dplyr", "terra", "vegan", "gdm"))
  bristol_require_outputs(cfg, c("manuscript_tables.rds", "climate_catalog.rds", "gdm_catalog.rds"), "analysis and manuscript_tables(cfg)")
  marker <- bristol_output(cfg, "publication_complete.rds")
  tables <- readRDS(file.path(cfg$out, "manuscript_tables.rds"))
  message("[Publication] Fossils and richness")
  map <- publication_fossils(cfg, tables)
  message("[Publication] Climate maps and stage trends")
  publication_climate(cfg)
  message("[Publication] Turnover, NMDS and temperature contours")
  publication_turnover(cfg, tables, map)
  message("[Publication] GDM splines and fit comparisons")
  publication_gdm(cfg, tables)
  files <- list.files(cfg$publication_figs, recursive = TRUE, full.names = TRUE)
  manuscript_write_table(data.frame(file = files, bytes = file.info(files)$size), "publication_inventory", cfg)
  saveRDS(list(config = cfg, session = utils::sessionInfo(), time = Sys.time()), marker)
  invisible(files)
}
