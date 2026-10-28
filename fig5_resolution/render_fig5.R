# Fig. 5: temperature and precipitation GDM contributions across resolutions.
# Source from the repository root; rendering is explicit.
source("figure_style/publication_style.R")

fig5_data <- function(cfg, guilds = c("all", "herb", "carn")) {
  path <- file.path(cfg$tables, "all_gdm_spline_coordinates.csv")
  if (!file.exists(path)) {
    stop("Missing saved GDM spline coordinates: ", path, call. = FALSE)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("model", "resolution", "continent", "guild", "variant",
                "geo", "predictor", "y")
  if (!all(required %in% names(x))) {
    stop("Unexpected GDM spline columns: ", path, call. = FALSE)
  }

  metrics <- list(
    Mean = c("temperature_mean", "precipitation_mean"),
    Change = c("temperature_signed_change", "precipitation_signed_change"),
    Seasonality = c("temperature_seasonal_change", "precipitation_seasonal_change")
  )
  predictors <- unlist(metrics, use.names = FALSE)
  keep <- x$resolution %in% c("h3", "h2", "h1", "h0") &
    x$variant == "combined" & as.logical(x$geo) & x$guild %in% guilds &
    x$continent %in% c("Europe", "North America") &
    x$predictor %in% predictors
  x <- x[!is.na(keep) & keep, ]
  if (!nrow(x) || any(!is.finite(x$y))) {
    stop("No complete combined-model climate spline data with geography included.",
         call. = FALSE)
  }

  # Maximum cumulative I-spline contribution over each predictor's fitted range.
  summary <- stats::aggregate(
    y ~ model + resolution + continent + guild + predictor,
    data = x, FUN = max
  )
  summary$hypothesis <- vapply(summary$predictor, function(p) {
    names(metrics)[vapply(metrics, function(fields) p %in% fields, logical(1))][1]
  }, character(1))
  summary$variable <- ifelse(grepl("^temperature_", summary$predictor),
    "Temperature", "Precipitation")
  summary <- summary[c("resolution", "continent", "guild", "hypothesis", "variable", "y")]
  names(summary)[names(summary) == "y"] <- "contribution"

  keys <- c("resolution", "continent", "guild", "hypothesis", "variable")
  expected <- expand.grid(resolution = c("h3", "h2", "h1", "h0"),
    continent = c("Europe", "North America"), guild = guilds,
    hypothesis = c("Mean", "Change", "Seasonality"),
    variable = c("Temperature", "Precipitation"), stringsAsFactors = FALSE)
  key <- function(d) do.call(paste, c(lapply(d[keys], as.character), sep = "__"))
  if (anyDuplicated(summary[keys]) || !setequal(key(summary), key(expected)) ||
      anyNA(summary$hypothesis) || any(!is.finite(summary$contribution))) {
    stop("Expected 48 complete climate contributions per requested guild across H3-H0.",
         call. = FALSE)
  }
  summary$resolution <- factor(summary$resolution,
    levels = c("h3", "h2", "h1", "h0"))
  summary$continent <- factor(summary$continent,
    levels = c("Europe", "North America"))
  summary$hypothesis <- factor(summary$hypothesis,
    levels = c("Mean", "Change", "Seasonality"))
  summary$variable <- factor(summary$variable,
    levels = c("Temperature", "Precipitation"))
  summary$guild_label <- factor(summary$guild,
    levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  summary$spacing_km <- unname(cfg$spacing[as.character(summary$resolution)])
  if (any(!is.finite(summary$spacing_km))) {
    stop("Missing H3 spacing in cfg$spacing", call. = FALSE)
  }
  summary
}

fig5_panel <- function(data, region, metric, upper, tag = NULL,
                       secondary_alpha = .45) {
  panel_data <- data[as.character(data$continent) == region &
    as.character(data$hypothesis) == metric, ]
  if (nrow(panel_data) != 8L * length(unique(data$guild))) {
    stop("Expected eight contributions per guild for ", region, " / ", metric,
      call. = FALSE)
  }
  resolution_levels <- c("h3", "h2", "h1", "h0")
  panel_data$resolution <- factor(as.character(panel_data$resolution),
    levels = resolution_levels)
  panel_data <- panel_data[order(panel_data$guild_label, panel_data$variable,
    panel_data$resolution), ]
  # Draw the all-mammal groups last, above the lighter guild curves.
  panel_data$curve_group <- interaction(panel_data$variable,
    factor(panel_data$guild, levels = c("herb", "carn", "all")), drop = TRUE)
  panel_data <- panel_data[order(panel_data$curve_group, panel_data$resolution), ]
  spacing <- panel_data$spacing_km[
    match(resolution_levels, as.character(panel_data$resolution))]
  resolution_labels <- stats::setNames(
    paste0(toupper(resolution_levels), "\n", round(spacing), " km"),
    resolution_levels)
  dodge <- ggplot2::position_dodge(width = .24)
  guild_breaks <- levels(panel_data$guild_label)[
    levels(panel_data$guild_label) %in% as.character(panel_data$guild_label)]

  ggplot2::ggplot(panel_data, ggplot2::aes(
    x = resolution, y = contribution, colour = guild_label,
    linetype = variable, alpha = guild_label, group = curve_group
  )) +
    ggplot2::geom_line(ggplot2::aes(linewidth = guild_label), position = dodge) +
    ggplot2::geom_point(ggplot2::aes(size = guild_label), position = dodge) +
    ggplot2::scale_alpha_manual(values = c("All mammals" = 1,
      Herbivores = secondary_alpha, Carnivores = secondary_alpha), guide = "none") +
    ggplot2::scale_linewidth_manual(values = c("All mammals" = 1,
      Herbivores = .55, Carnivores = .55), guide = "none") +
    ggplot2::scale_size_manual(values = c("All mammals" = 2.2,
      Herbivores = 1.5, Carnivores = 1.5), guide = "none") +
    ggplot2::scale_colour_manual(
      values = figpub_colours$guild, breaks = guild_breaks,
      labels = guild_breaks, name = "Guild", drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = c(Temperature = "solid", Precipitation = "dashed"),
      breaks = c("Temperature", "Precipitation"),
      labels = c("Temperature", "Precipitation"), name = "Climate", drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, nrow = 1,
        override.aes = list(linetype = "solid", alpha = 1, linewidth = .8, size = 2)),
      linetype = ggplot2::guide_legend(order = 2, nrow = 1,
        override.aes = list(colour = "#4D4D4D", alpha = 1, linewidth = .8, size = 2))
    ) +
    ggplot2::scale_x_discrete(
      limits = resolution_levels, breaks = resolution_levels,
      labels = resolution_labels, drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      name = paste0("GDM contribution (", metric, ")"),
      breaks = function(limits) {
        ticks <- pretty(c(0, upper), n = 4)
        ticks[ticks >= 0 & ticks <= upper]
      },
      expand = ggplot2::expansion(mult = c(.06, .08))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::labs(
      title = region, tag = tag,
      x = "Spatial resolution"
    ) +
    figpub_theme() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(size = 13, face = "bold"),
      plot.tag.position = c(.01, .99),
      legend.position = "bottom", legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(1.4, "lines"),
      legend.key.height = grid::unit(.8, "lines"),
      legend.spacing.x = grid::unit(2, "pt"),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(0, 2, 0, 2),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(4, 4, 3, 4)
    )
}

render_fig5 <- function(cfg = bristol_config(), out_dir = "fig5_resolution/output",
                        overwrite = FALSE, guilds = c("all", "herb", "carn"),
                        layout = c("combined", "separate"), secondary_alpha = .45) {
  for (package in c("ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required package: ", package, call. = FALSE)
    }
  }
  if (!is.character(guilds) || !length(guilds) || anyNA(guilds) ||
      anyDuplicated(guilds) || !all(guilds %in% c("all", "herb", "carn")))
    stop("guilds must contain unique values from all, herb, and carn.")
  layout <- match.arg(layout)
  if (!is.numeric(secondary_alpha) || length(secondary_alpha) != 1L ||
      !is.finite(secondary_alpha) || secondary_alpha < 0 || secondary_alpha > 1)
    stop("secondary_alpha must be a number between 0 and 1.")
  guild_labels <- c(all = "All mammals", herb = "Herbivores", carn = "Carnivores")
  stems <- c(all = "fig05_spatial_resolution_all_mammals",
    herb = "fig05_spatial_resolution_herbivores",
    carn = "fig05_spatial_resolution_carnivores")
  export_stems <- if (layout == "combined") "fig05_spatial_resolution" else stems[guilds]
  planned <- c(unlist(lapply(export_stems, function(stem)
    file.path(out_dir, paste0(stem, c(".png", ".pdf"))))),
    file.path(out_dir, "fig05_spatial_resolution_data.csv"))
  if (!overwrite && any(file.exists(planned))) {
    stop("Fig. 5 output already exists; set overwrite=TRUE or choose another out_dir",
         call. = FALSE)
  }

  data <- fig5_data(cfg, guilds = guilds)
  upper <- max(data$contribution)
  if (upper <= 0) upper <- 1
  panel_keys <- expand.grid(
    metric = c("Mean", "Change", "Seasonality"),
    region = c("Europe", "North America"), stringsAsFactors = FALSE)
  groups <- if (layout == "combined") "combined" else guilds
  figures <- stats::setNames(vector("list", length(groups)), groups)
  for (group in groups) {
    plot_data <- if (group == "combined") data else data[data$guild == group, ]
    panels <- lapply(seq_len(nrow(panel_keys)), function(i) {
      fig5_panel(plot_data, panel_keys$region[i], panel_keys$metric[i], upper,
        tag = LETTERS[i], secondary_alpha = secondary_alpha)
    })
    title <- "Climate contributions across spatial scales"
    if (group != "combined") title <- paste(title, guild_labels[[group]], sep = "\n")
    figure <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
      figpub_annotation(title) &
      ggplot2::theme(legend.position = "bottom", legend.box = "horizontal",
        legend.direction = "horizontal")
    stem <- if (group == "combined") "fig05_spatial_resolution" else stems[[group]]
    figpub_save(out_dir, stem, figure,
      height = 6.5, overwrite = overwrite)
    figures[[group]] <- figure
  }
  utils::write.csv(data, file.path(out_dir, "fig05_spatial_resolution_data.csv"),
    row.names = FALSE)
  reference <- if (layout == "combined") "combined" else
    if ("all" %in% guilds) "all" else guilds[1L]
  invisible(list(figure = figures[[reference]], figures = figures, data = data))
}
