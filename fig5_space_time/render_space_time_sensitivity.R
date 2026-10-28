# Space/time counterpart of fig5_resolution/render_fig5.R.
# Source from the repository root; rendering is explicit.
source("figure_style/publication_style.R")

space_time_sensitivity_data <- function(cfg, guilds = c("all", "herb", "carn")) {
  path <- file.path(cfg$tables, "all_gdm_spline_coordinates.csv")
  if (!file.exists(path)) stop("Missing saved spline coordinates: ", path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("model", "resolution", "continent", "guild", "variant",
    "geo", "predictor", "y")
  if (!all(required %in% names(x))) stop("Unexpected spline columns: ", path)
  resolutions <- c("h3", "h2", "h1", "h0")
  regions <- c("Europe", "North America")
  keep <- x$resolution %in% resolutions & x$continent %in% regions &
    x$guild %in% guilds & x$variant == "combined" & as.logical(x$geo) &
    x$predictor %in% c("Geographic", "time_ma")
  x <- x[!is.na(keep) & keep, ]
  if (!nrow(x) || any(!is.finite(x$y)))
    stop("Missing or non-finite space/time spline contributions.")

  # Same statistic as the climate sensitivity: maximum cumulative spline height.
  data <- stats::aggregate(y ~ model + resolution + continent + guild + predictor,
    data = x, FUN = max)
  expected <- expand.grid(resolution = resolutions, continent = regions,
    guild = guilds, predictor = c("Geographic", "time_ma"), stringsAsFactors = FALSE)
  keys <- c("resolution", "continent", "guild", "predictor")
  key <- function(d) do.call(paste, c(lapply(d[keys], as.character), sep = "__"))
  if (anyDuplicated(data[keys]) || !setequal(key(data), key(expected)))
    stop("Expected 16 complete space/time contributions per requested guild.")
  names(data)[names(data) == "y"] <- "contribution"
  data$component <- factor(ifelse(data$predictor == "Geographic", "Space", "Time"),
    levels = c("Space", "Time"))
  data$resolution <- factor(data$resolution, levels = resolutions)
  data$continent <- factor(data$continent, levels = regions)
  data$guild_label <- factor(data$guild, levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  data$spacing_km <- unname(cfg$spacing[as.character(data$resolution)])
  if (any(!is.finite(data$spacing_km))) stop("Missing spacing in cfg$spacing.")
  data[order(data$continent, data$guild_label, data$component, data$resolution), ]
}

space_time_sensitivity_panel <- function(data, region, upper, tag,
                                         secondary_alpha = .45) {
  d <- data[as.character(data$continent) == region, ]
  if (nrow(d) != 8L * length(unique(data$guild)))
    stop("Expected eight space/time contributions per guild for ", region)
  d$curve_group <- interaction(d$component,
    factor(d$guild, levels = c("herb", "carn", "all")), drop = TRUE)
  d <- d[order(d$curve_group, d$resolution), ]
  resolutions <- c("h3", "h2", "h1", "h0")
  spacing <- d$spacing_km[match(resolutions, as.character(d$resolution))]
  labels <- stats::setNames(paste0(toupper(resolutions), "\n",
    round(spacing), " km"), resolutions)
  dodge <- ggplot2::position_dodge(width = .24)
  guild_breaks <- levels(d$guild_label)[
    levels(d$guild_label) %in% as.character(d$guild_label)]
  ggplot2::ggplot(d, ggplot2::aes(resolution, contribution,
      colour = guild_label, linetype = component, alpha = guild_label,
      group = curve_group)) +
    ggplot2::geom_line(ggplot2::aes(linewidth = guild_label), position = dodge) +
    ggplot2::geom_point(ggplot2::aes(size = guild_label), position = dodge) +
    ggplot2::scale_alpha_manual(values = c("All mammals" = 1,
      Herbivores = secondary_alpha, Carnivores = secondary_alpha), guide = "none") +
    ggplot2::scale_linewidth_manual(values = c("All mammals" = 1,
      Herbivores = .55, Carnivores = .55), guide = "none") +
    ggplot2::scale_size_manual(values = c("All mammals" = 2.2,
      Herbivores = 1.5, Carnivores = 1.5), guide = "none") +
    ggplot2::scale_colour_manual(values = figpub_colours$guild,
      breaks = guild_breaks, labels = guild_breaks, name = "Guild", drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(Space = "solid", Time = "dashed"),
      breaks = c("Space", "Time"), labels = c("Space", "Time"),
      name = "Distance", drop = FALSE) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, nrow = 1,
        override.aes = list(linetype = "solid", alpha = 1, linewidth = .8, size = 2)),
      linetype = ggplot2::guide_legend(order = 2, nrow = 1,
        override.aes = list(colour = "#4D4D4D", alpha = 1, linewidth = .8, size = 2))) +
    ggplot2::scale_x_discrete(limits = resolutions, breaks = resolutions,
      labels = labels, drop = FALSE) +
    ggplot2::scale_y_continuous(name = "GDM contribution",
      breaks = function(limits) {
        ticks <- pretty(c(0, upper), n = 4)
        ticks[ticks >= 0 & ticks <= upper]
      }, expand = ggplot2::expansion(mult = c(.06, .08))) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::labs(title = region, x = "Spatial resolution", tag = tag) +
    figpub_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
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
      plot.margin = ggplot2::margin(4, 4, 3, 4))
}

render_space_time_sensitivity <- function(cfg = bristol_config(),
    out_dir = "fig5_space_time/output", overwrite = FALSE,
    guilds = c("all", "herb", "carn"), secondary_alpha = .45) {
  for (package in c("ggplot2", "patchwork"))
    if (!requireNamespace(package, quietly = TRUE))
      stop("Required package: ", package)
  if (!is.character(guilds) || !length(guilds) || anyNA(guilds) ||
      anyDuplicated(guilds) || !all(guilds %in% c("all", "herb", "carn")))
    stop("guilds must contain unique values from all, herb, and carn.")
  if (!is.numeric(secondary_alpha) || length(secondary_alpha) != 1L ||
      !is.finite(secondary_alpha) || secondary_alpha < 0 || secondary_alpha > 1)
    stop("secondary_alpha must be a number between 0 and 1.")
  stem <- "space_time_resolution_sensitivity"
  planned <- file.path(out_dir, paste0(stem, c(".png", ".pdf", "_data.csv")))
  if (!overwrite && any(file.exists(planned)))
    stop("Space/time sensitivity output exists; set overwrite=TRUE or choose another out_dir.")
  data <- space_time_sensitivity_data(cfg, guilds = guilds)
  upper <- max(data$contribution)
  if (upper <= 0) upper <- 1
  panels <- list(
    space_time_sensitivity_panel(data, "Europe", upper, "A", secondary_alpha),
    space_time_sensitivity_panel(data, "North America", upper, "B", secondary_alpha))
  figure <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
    figpub_annotation("Space and time contributions across spatial scales") &
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal",
      legend.direction = "horizontal")
  figpub_save(out_dir, stem, figure, height = 4.2, overwrite = overwrite)
  utils::write.csv(data, planned[3], row.names = FALSE)
  invisible(list(figure = figure, data = data))
}
