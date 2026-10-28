# Fig. 4: conditional GDM prediction contrasts and supplementary I-splines.
# Source from the repository root; rendering is explicit.
source("figure_style/publication_style.R")

fig4_required <- function(path) {
  if (!file.exists(path)) stop("Missing saved Bristol output: ", path, call. = FALSE)
  path
}

fig4_weighted_quantile <- function(x, w, probs = c(.25, .75)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]
  order_id <- order(x)
  x <- x[order_id]; w <- w[order_id]
  cumulative <- cumsum(w)/sum(w)
  x[vapply(probs, function(p) which(cumulative >= p)[1], integer(1))]
}

fig4_model_contrasts <- function(path, region, guild, resolution) {
  object <- readRDS(fig4_required(path))
  if (!identical(object$variant, "combined") || !isTRUE(object$geo))
    stop("Fig. 4 requires the combined model with geographic distance: ", path)
  table <- object$gdm_table
  if (!all(c("weights", "s1.xCoord", "s2.xCoord") %in% names(table)))
    stop("Unexpected saved GDM site-pair table: ", path)
  weights <- table$weights
  if (any(!is.finite(weights)) || any(weights <= 0)) stop("Invalid site-pair weights: ", path)
  full <- as.numeric(stats::predict(object$model, data = table))
  if (length(full) != nrow(table) || any(!is.finite(full)))
    stop("Invalid full-model predictions: ", path)
  rows <- list()
  for (metric in c("mean", "signed_change", "seasonal_change")) {
    for (scenario in c("Temperature", "Precipitation", "Both")) {
      variables <- switch(scenario, Temperature = "temperature",
        Precipitation = "precipitation", Both = c("temperature", "precipitation"))
      fields <- paste(variables, metric, sep = "_")
      left <- paste0("s1.", fields)
      right <- paste0("s2.", fields)
      if (!all(c(left, right) %in% names(table)))
        stop("Focal predictor absent from saved model: ", paste(fields, collapse = ", "), " in ", path)
      counterfactual <- table
      for (j in seq_along(fields)) counterfactual[[right[j]]] <- table[[left[j]]]
      equalized <- as.numeric(stats::predict(object$model, data = counterfactual))
      if (length(equalized) != nrow(table) || any(!is.finite(equalized)))
        stop("Invalid equalized predictions: ", path)
      delta <- full - equalized
      q <- fig4_weighted_quantile(delta, weights)
      rows[[length(rows) + 1L]] <- data.frame(
        continent = if (region == "europe") "Europe" else "North America",
        guild = guild, resolution = resolution, hypothesis = metric, scenario = scenario,
        mean_delta = stats::weighted.mean(delta, weights),
        q25 = q[1], q75 = q[2],
        mean_full = stats::weighted.mean(full, weights),
        mean_equalized = stats::weighted.mean(equalized, weights),
        n_pairs = length(delta), sum_weights = sum(weights),
        model_file = path)
    }
  }
  do.call(rbind, rows)
}

fig4_main_plot <- function(contrasts) {
  required <- expand.grid(continent = c("Europe", "North America"),
    guild = c("all", "herb", "carn"), resolution = c("h3", "h2", "h1", "h0"),
    hypothesis = c("mean", "signed_change", "seasonal_change"),
    scenario = c("Temperature", "Precipitation", "Both"),
    stringsAsFactors = FALSE)
  keys <- c("continent", "guild", "resolution", "hypothesis", "scenario")
  if (!all(keys %in% names(contrasts)) || nrow(contrasts) != nrow(required) ||
      anyDuplicated(contrasts[keys]) ||
      !setequal(do.call(paste, contrasts[keys]), do.call(paste, required[keys])))
    stop("Fig. 4 requires all four H3 resolutions for every continent, guild, metric, and scenario")
  contrasts$continent <- factor(contrasts$continent, levels = c("Europe", "North America"))
  contrasts$guild <- factor(contrasts$guild, levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  contrasts$hypothesis <- factor(contrasts$hypothesis,
    levels = c("mean", "signed_change", "seasonal_change"),
    labels = c("Mean", "Change", "Seasonality"))
  contrasts$scenario <- factor(contrasts$scenario,
    levels = c("Temperature", "Precipitation", "Both"))
  contrasts$resolution <- factor(contrasts$resolution, levels = c("h3", "h2", "h1", "h0"))
  reference <- contrasts[contrasts$resolution == "h3", ]
  climate_colours <- c(Temperature = "#7B3294", Precipitation = "#008C7A",
    Both = "#595959")
  ggplot2::ggplot(reference,
    ggplot2::aes(hypothesis, mean_delta, fill = scenario)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = .35) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = .8),
      width = .72, colour = "white", linewidth = .2) +
    ggplot2::facet_grid(continent ~ guild) +
    ggplot2::scale_fill_manual(values = climate_colours, name = "Climate equalized") +
    ggplot2::labs(x = "Climate metric", y = "Predicted turnover contrast",
      title = "Climate contrasts in predicted turnover",
      caption = "Bars: fitting-weighted means from the h3 reference models.") +
    figpub_theme()
}

fig4_spline_figure <- function(data, variable) {
  subset <- data[data$variable == variable, ]
  if (!nrow(subset)) stop("No spline data for ", variable)
  y_max <- max(subset$y, na.rm = TRUE)
  subset$continent <- factor(subset$continent,
    levels = c("Europe", "North America"))
  subset$metric <- factor(subset$metric,
    levels = c("mean", "signed_change", "seasonal_change"),
    labels = c("Mean", "Change", "Seasonality"))
  subset$guild <- factor(subset$guild, levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores"))
  tags <- expand.grid(metric = levels(subset$metric),
    continent = levels(subset$continent), stringsAsFactors = FALSE)
  tags$continent <- factor(tags$continent, levels = levels(subset$continent))
  tags$metric <- factor(tags$metric, levels = levels(subset$metric))
  tags$label <- LETTERS[seq_len(nrow(tags))]
  zero <- unique(subset[subset$metric == "Change", c("continent", "metric")])
  plot <- ggplot2::ggplot(subset,
      ggplot2::aes(x, y, colour = guild, group = guild)) +
    ggplot2::geom_line(linewidth = .7) +
    ggplot2::geom_text(data = tags,
      ggplot2::aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.35, vjust = 1.35, fontface = "bold", size = 4) +
    ggplot2::scale_colour_manual(values = figpub_colours$guild, name = NULL) +
    ggplot2::facet_grid(continent ~ metric, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = c(0, y_max)) +
    ggplot2::labs(title = paste("Bristol GDM I-splines:", variable),
      x = if (variable == "temperature") "Temperature (deg C)" else "Precipitation (mm/day)",
      y = "GDM contribution",
      caption = "Cumulative predictor contribution from combined temperature-precipitation models") +
    figpub_theme() +
    ggplot2::theme(legend.position = "bottom",
      strip.text.y = ggplot2::element_text(size = 10, face = "bold"))
  if (nrow(zero)) plot <- plot +
    ggplot2::geom_vline(data = zero, ggplot2::aes(xintercept = 0),
      inherit.aes = FALSE, linewidth = .3, linetype = "dashed", colour = "grey45")
  plot
}

fig4_read_splines <- function(cfg) {
  x <- utils::read.csv(fig4_required(file.path(cfg$tables, "all_gdm_spline_coordinates.csv")))
  x <- x[x$resolution == "h3" & x$variant == "combined" & x$geo %in% c(TRUE, "TRUE") &
    x$guild %in% c("all", "herb", "carn"), ]
  x <- x[grepl("^(temperature|precipitation)_(mean|signed_change|seasonal_change)$", x$predictor), ]
  x$variable <- sub("_.*$", "", x$predictor)
  x$metric <- sub("^(temperature|precipitation)_", "", x$predictor)
  if (length(unique(x$model)) != 6L || any(!is.finite(x$x)) || any(!is.finite(x$y)))
    stop("Incomplete H3 combined geographic-model spline coordinates")
  x
}

render_fig4 <- function(cfg = bristol_config(), out_dir = "fig4_predicted/output",
                        overwrite = FALSE) {
  for (package in c("gdm", "ggplot2", "patchwork"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package: ", package)
  planned <- file.path(out_dir, c("fig04_predicted_turnover.png", "fig04_predicted_turnover.pdf",
    "figS_temperature_gdm_splines.png", "figS_temperature_gdm_splines.pdf",
    "figS_precipitation_gdm_splines.png", "figS_precipitation_gdm_splines.pdf",
    "fig04_predicted_turnover_contrasts.csv"))
  if (!overwrite && any(file.exists(planned)))
    stop("Fig. 4 output already exists; set overwrite=TRUE or choose another out_dir")
  rows <- list()
  for (resolution in c("h3", "h2", "h1", "h0"))
    for (region in c("europe", "north_america")) for (guild in c("all", "herb", "carn")) {
    name <- paste0("gdm_", resolution, "__", region, "__", guild,
      "__combined__geoTRUE.rds")
    rows[[length(rows) + 1L]] <- fig4_model_contrasts(
      file.path(cfg$out, name), region, guild, resolution)
  }
  contrasts <- do.call(rbind, rows)
  main <- fig4_main_plot(contrasts)
  splines <- fig4_read_splines(cfg)
  supplements <- list(temperature = fig4_spline_figure(splines, "temperature"),
    precipitation = fig4_spline_figure(splines, "precipitation"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  figpub_save(out_dir, "fig04_predicted_turnover", main,
    height = 5.8, overwrite = overwrite)
  for (variable in names(supplements))
    figpub_save(out_dir, paste0("figS_", variable, "_gdm_splines"),
      supplements[[variable]], height = 7.2, overwrite = overwrite)
  utils::write.csv(contrasts, file.path(out_dir, "fig04_predicted_turnover_contrasts.csv"),
    row.names = FALSE)
  invisible(list(main = main, supplements = supplements, contrasts = contrasts))
}

# Redraw only the main figure after a visual edit, using the saved summary.
# This does not recompute GDM predictions or redraw supplementary splines.
render_fig4_main_from_csv <- function(out_dir = "fig4_predicted/output", overwrite = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required package: ggplot2")
  data <- utils::read.csv(fig4_required(file.path(out_dir,
    "fig04_predicted_turnover_contrasts.csv")))
  plot <- fig4_main_plot(data)
  figpub_save(out_dir, "fig04_predicted_turnover", plot,
    height = 5.8, overwrite = overwrite)
  invisible(plot)
}

# Restyle the main figure and both spline supplements from saved CSVs only.
render_fig4_from_saved <- function(cfg = bristol_config(),
                                   out_dir = "fig4_predicted/output", overwrite = FALSE) {
  for (package in c("ggplot2", "patchwork"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package: ", package)
  planned <- file.path(out_dir, paste0(c("fig04_predicted_turnover",
    "figS_temperature_gdm_splines", "figS_precipitation_gdm_splines"),
    rep(c(".png", ".pdf"), each = 3L)))
  if (!overwrite && any(file.exists(planned)))
    stop("Fig. 4 output already exists; set overwrite=TRUE or choose another out_dir")
  data <- utils::read.csv(fig4_required(file.path(out_dir,
    "fig04_predicted_turnover_contrasts.csv")))
  splines <- fig4_read_splines(cfg)
  main <- fig4_main_plot(data)
  supplements <- list(temperature = fig4_spline_figure(splines, "temperature"),
    precipitation = fig4_spline_figure(splines, "precipitation"))
  figpub_save(out_dir, "fig04_predicted_turnover", main,
    height = 5.8, overwrite = overwrite)
  for (variable in names(supplements))
    figpub_save(out_dir, paste0("figS_", variable, "_gdm_splines"),
      supplements[[variable]], height = 7.2, overwrite = overwrite)
  invisible(list(main = main, supplements = supplements))
}
