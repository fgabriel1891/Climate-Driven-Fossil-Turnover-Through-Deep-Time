# Six-panel Figure 04: h3 combined-model temperature and precipitation splines.
# Source from the repository root; rendering is explicit.
source("fig4_predicted/render_fig4.R")

fig04_h3_spline_panel <- function(data, variable, metric, tag, upper) {
  d <- data[data$variable == variable, , drop = FALSE]
  if (!nrow(d)) stop("No spline coordinates found for ", variable, call. = FALSE)
  d <- d[d$metric == metric, , drop = FALSE]
  if (!nrow(d)) stop("No ", metric, " spline coordinates for ", variable, call. = FALSE)
  expected <- expand.grid(
    continent = c("Europe", "North America"),
    guild = c("all", "herb", "carn"),
    stringsAsFactors = FALSE)
  keys <- c("continent", "guild")
  observed <- unique(d[keys])
  if (nrow(observed) != nrow(expected) ||
      !setequal(do.call(paste, observed), do.call(paste, expected)))
    stop("Incomplete h3 spline curves for ", variable, " / ", metric, call. = FALSE)

  d$continent <- factor(d$continent, levels = c("Europe", "North America"))
  d$guild <- factor(d$guild, levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  metric_label <- c(mean = "Mean", signed_change = "Change",
    seasonal_change = "Seasonality")[[metric]]
  x_title <- if (variable == "temperature") "Temperature (deg C)" else "Precipitation (mm/day)"
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = guild,
      linetype = continent, group = interaction(guild, continent))) +
    ggplot2::geom_line(linewidth = .75) +
    ggplot2::scale_colour_manual(values = figpub_colours$guild, name = "Guild") +
    ggplot2::scale_linetype_manual(values = c(Europe = "solid", `North America` = "dashed"),
      name = "Continent") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.02, .05))) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::labs(x = x_title, y = "Cumulative GDM contribution",
      title = metric_label, tag = tag) +
    figpub_theme() +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 13, face = "bold"),
      plot.tag.position = c(.01, .99),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom")
  if (metric == "signed_change") p <- p + ggplot2::geom_vline(
    xintercept = 0, linewidth = .3, linetype = "dashed", colour = "grey45")
  p
}

fig04_h3_spline_figure <- function(data) {
  required <- c("resolution", "variant", "geo", "variable", "continent", "guild", "metric", "x", "y")
  if (!all(required %in% names(data))) stop("Unexpected spline coordinate table", call. = FALSE)
  d <- data[data$resolution == "h3" & data$variant == "combined" &
      data$geo %in% c(TRUE, "TRUE") & data$variable %in% c("temperature", "precipitation") &
      data$guild %in% c("all", "herb", "carn"), , drop = FALSE]
  if (any(!is.finite(d$x)) || any(!is.finite(d$y)))
    stop("Non-finite spline coordinates", call. = FALSE)

  panel_keys <- data.frame(
    variable = c(rep("temperature", 3), rep("precipitation", 3)),
    metric = rep(c("mean", "signed_change", "seasonal_change"), 2),
    stringsAsFactors = FALSE)
  panels <- lapply(seq_len(nrow(panel_keys)), function(i) {
    variable <- panel_keys$variable[i]
    upper <- max(d$y[d$variable == variable], na.rm = TRUE)
    if (upper <= 0) upper <- 1
    fig04_h3_spline_panel(d, variable, panel_keys$metric[i], LETTERS[i], upper)
  })
  figure <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    figpub_annotation(
      "Climate contributions to fossil turnover",
      caption = "h3 combined GDMs with geographic distance. The top row shows temperature and the bottom row precipitation. Lines show cumulative I-spline contributions, not predicted turnover values; colours distinguish guilds and line types distinguish continents.")
  figure & ggplot2::theme(
    legend.position = "bottom", legend.box = "horizontal",
    legend.direction = "horizontal")
}

render_fig04_h3_splines <- function(cfg = bristol_config(),
    out_dir = "fig04_h3_splines/output", overwrite = FALSE) {
  for (package in c("ggplot2", "patchwork"))
    if (!requireNamespace(package, quietly = TRUE))
      stop("Required package: ", package, call. = FALSE)
  data <- fig4_read_splines(cfg)
  figure <- fig04_h3_spline_figure(data)
  figpub_save(out_dir, "fig04_h3_climate_splines", figure,
    height = 6.5, width = 7.5, dpi = 400, overwrite = overwrite)
  invisible(figure)
}
