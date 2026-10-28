# Response-scale test of the three conceptual climate-turnover hypotheses.
# Run from the repository root after the H3 combined models have been fitted.
source("figure_style/publication_style.R")
hypothesis_guild_colours <- c(
  "All mammals" = unname(figpub_colours$guild[["All mammals"]]),
  "Herbivores" = unname(figpub_colours$guild[["Herbivores"]]),
  "Carnivores" = unname(figpub_colours$guild[["Carnivores"]]))
hypothesis_metric_labels <- c(mean = "Mean climate",
  signed_change = "Climate change", seasonal_change = "Seasonal variability")

hypothesis_prediction_background <- function(reference_object) {
  tab <- reference_object$gdm_table
  geographic <- sqrt((tab$s1.xCoord - tab$s2.xCoord)^2 +
    (tab$s1.yCoord - tab$s2.yCoord)^2)
  temporal <- abs(tab$s1.time_ma - tab$s2.time_ma)
  list(geographic_km = stats::median(geographic, na.rm = TRUE),
       temporal_ma = stats::median(temporal, na.rm = TRUE))
}

hypothesis_predictor_values <- function(object, field) {
  if (field %in% names(object$site_data)) return(object$site_data[[field]])
  tab <- object$gdm_table
  endpoints <- paste0(c("s1.", "s2."), field)
  if (!all(endpoints %in% names(tab)))
    stop("Predictor values are unavailable for ", field)
  c(tab[[endpoints[1]]], tab[[endpoints[2]]])
}

hypothesis_prediction_curve <- function(object, background, field, metric, n = 151L,
    orientation = 1) {
  tab <- object$gdm_table
  sites <- object$site_data
  predictor_values <- hypothesis_predictor_values(object, field)
  bounds <- range(predictor_values, na.rm = TRUE)
  if (length(bounds) != 2L || any(!is.finite(bounds))) stop("Invalid predictor range: ", field)
  width <- diff(bounds)
  center <- mean(bounds)
  signed <- identical(metric, "signed_change")
  contrast <- if (signed) {
    seq(-width, width, length.out = n)
  } else {
    seq(0, width, length.out = n)
  }

  # Start with the exact saved GDM schema, then construct a synthetic pair whose
  # non-focal climate predictors are equal and whose other backgrounds are typical.
  base <- tab[1L, , drop = FALSE]
  base$distance <- 0
  base$weights <- 1
  fields <- setdiff(sub("^s1[.]", "", grep("^s1[.]", names(base), value = TRUE)),
    c("xCoord", "yCoord"))
  for (p in fields) {
    value <- stats::median(hypothesis_predictor_values(object, p), na.rm = TRUE)
    base[[paste0("s1.", p)]] <- value
    base[[paste0("s2.", p)]] <- value
  }

  distances <- sqrt((tab$s1.xCoord - tab$s2.xCoord)^2 +
    (tab$s1.yCoord - tab$s2.yCoord)^2)
  geo <- max(min(distances, na.rm = TRUE),
    min(max(distances, na.rm = TRUE), background$geographic_km))
  base$s1.xCoord <- 0
  base$s1.yCoord <- 0
  base$s2.xCoord <- geo
  base$s2.yCoord <- 0

  temporal <- NA_real_
  if ("time_ma" %in% fields) {
    age_bounds <- range(sites$time_ma, na.rm = TRUE)
    temporal <- min(background$temporal_ma, diff(age_bounds))
    midpoint <- mean(age_bounds)
    base$s1.time_ma <- midpoint + temporal / 2
    base$s2.time_ma <- midpoint - temporal / 2
  }

  included <- all(paste0(c("s1.", "s2."), field) %in% names(base))
  raw_contrast <- contrast * orientation
  left <- center - raw_contrast / 2
  right <- center + raw_contrast / 2
  newdata <- base[rep(1L, length(contrast)), , drop = FALSE]
  if (included) {
    newdata[[paste0("s1.", field)]] <- left
    newdata[[paste0("s2.", field)]] <- right
  }
  if (!identical(names(newdata), names(tab))) stop("Prediction table columns do not match saved GDM")
  prediction <- as.numeric(stats::predict(object$model, data = newdata))
  if (length(prediction) != length(contrast) || any(!is.finite(prediction)))
    stop("Invalid GDM prediction for ", field)
  if (any(prediction < 0 | prediction > 1))
    stop("Predictions fell outside the Simpson turnover scale for ", field)

  data.frame(x = contrast, predicted_turnover = prediction,
    site1_predictor = left, site2_predictor = right,
    predictor_min = bounds[1], predictor_max = bounds[2],
    predictor_in_model = included, geographic_km = geo, temporal_ma = temporal)
}

hypothesis_test_panel <- function(data, metric, continent, variable) {
  d <- data[data$metric == metric & data$continent == continent, ]
  if (!nrow(d)) stop("Missing prediction curves for ", continent, " / ", metric)
  # Draw the guild curves first so the heavier All mammals line stays on top.
  d$guild_label <- factor(d$guild, levels = c("herb", "carn", "all"),
    labels = c("Herbivores", "Carnivores", "All mammals"))
  metric_labels <- if (variable == "temperature") {
    c(mean = "Mean temperature", signed_change = "Signed temperature change",
      seasonal_change = "Temperature seasonal variability")
  } else {
    c(mean = "Mean precipitation", signed_change = "Signed precipitation change",
      seasonal_change = "Precipitation seasonal variability")
  }
  metric_label <- metric_labels[[metric]]
  units <- if (variable == "temperature") "°C" else "mm/day"
  ggplot2::ggplot(d, ggplot2::aes(x, predicted_turnover,
      colour = guild_label, linewidth = guild_label, alpha = guild_label,
      group = guild_label)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = hypothesis_guild_colours,
      breaks = c("All mammals", "Herbivores", "Carnivores"), name = NULL) +
    ggplot2::scale_linewidth_manual(values = c("All mammals" = 1.15,
      Herbivores = .65, Carnivores = .65), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("All mammals" = 1,
      Herbivores = .68, Carnivores = .68), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(.7, 1), breaks = seq(.7, 1, .05),
      expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(title = continent, x = paste0(metric_label, " (", units, ")"),
      y = "Predicted Simpson turnover") + figpub_theme() 
}

hypothesis_pca_extremes <- function(cfg, region, guild, axis, metric) {
  path <- file.path(cfg$out, paste0("pca_loadings_h3__", region, "__", guild, ".csv"))
  if (!file.exists(path)) stop("Missing saved PCA loadings: ", path)
  loadings <- utils::read.csv(path, stringsAsFactors = FALSE)
  loadings <- loadings[loadings$axis == axis & loadings$retained, ]
  temp <- loadings[grepl("^temperature_", loadings$variable), ]
  precip <- loadings[grepl("^precipitation_", loadings$variable), ]
  if (nrow(temp) != 1L || nrow(precip) != 1L)
    stop("Expected one temperature and one precipitation loading for ", axis)

  directions <- switch(metric,
    mean = list(temperature = c("warmer", "cooler"),
                precipitation = c("wetter", "drier")),
    signed_change = list(
      temperature = c("higher T change (warming direction)",
                      "lower T change (cooling direction)"),
      precipitation = c("higher P change (wetting direction)",
                        "lower P change (drying direction)")),
    seasonal_change = list(
      temperature = c("greater T seasonal-pattern change",
                      "smaller T seasonal-pattern change"),
      precipitation = c("greater P seasonal-pattern change",
                        "smaller P seasonal-pattern change")))
  # PC signs are arbitrary. Reorient each PC1 so increasing scores mean a
  # positive temperature loading; GDM pair distances are unchanged by reflection.
  orientation <- if (temp$loading >= 0) 1 else -1
  precip_loading <- precip$loading * orientation
  toward_high <- function(loading, labels) if (loading >= 0) labels[1] else labels[2]
  toward_low <- function(loading, labels) if (loading >= 0) labels[2] else labels[1]
  label <- paste0("Low PC1: ", toward_low(1, directions$temperature), "; ",
    toward_low(precip_loading, directions$precipitation),
    "\nHigh PC1: ", toward_high(1, directions$temperature), "; ",
    toward_high(precip_loading, directions$precipitation))
  list(orientation = orientation, label = label)
}

hypothesis_pca_panel <- function(data, metric, continent) {
  d <- data[data$axis == "PC1", ]
  d <- d[d$metric == metric & d$continent == continent, ]
  if (!nrow(d)) stop("Missing PC1 prediction curves for ", continent, " / ", metric)
  d$guild_label <- factor(d$guild, levels = c("herb", "carn", "all"),
    labels = c("Herbivores", "Carnivores", "All mammals"))
 
  metric_label <- hypothesis_metric_labels[[metric]]
   
   
  ggplot2::ggplot(d, ggplot2::aes(x, predicted_turnover, colour = guild_label,
      linewidth = guild_label, alpha = guild_label, group = guild_label)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = hypothesis_guild_colours,
      breaks = c("All mammals", "Herbivores", "Carnivores"), name = NULL) +
    ggplot2::scale_linewidth_manual(values = c("All mammals" = 1.15,
      Herbivores = .65, Carnivores = .65), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("All mammals" = 1,
      Herbivores = .68, Carnivores = .68), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(.7, 1), breaks = seq(.7, 1, .05),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = continent,
      x = paste0(metric_label, "-PC1"),
      y = "Predicted Simpson turnover") + figpub_theme()
}

dual_axis_hypothesis_plot <- function(curves) {
  temperature_max <- max(abs(curves$x[curves$variable == "temperature"]), na.rm = TRUE)
  precipitation_max <- max(abs(curves$x[curves$variable == "precipitation"]), na.rm = TRUE)
  if (!is.finite(temperature_max) || temperature_max <= 0 ||
      !is.finite(precipitation_max) || precipitation_max <= 0)
    stop("Cannot build dual x-axes from zero or missing contrast ranges")

  d <- curves
  d$metric <- factor(d$metric, levels = names(hypothesis_metric_labels),
    labels = unname(hypothesis_metric_labels))
  d$guild_label <- factor(d$guild, levels = c("herb", "carn", "all"),
    labels = c("Herbivores", "Carnivores", "All mammals"))
  d$x_norm <- ifelse(d$variable == "temperature",
    d$x / temperature_max, d$x / precipitation_max)
  temperature <- d[d$variable == "temperature", ]
  precipitation <- d[d$variable == "precipitation", ]

  ggplot2::ggplot(temperature, ggplot2::aes(x_norm, predicted_turnover,
      colour = guild_label, linewidth = guild_label, alpha = guild_label,
      group = guild_label)) +
    ggplot2::geom_line(ggplot2::aes(linetype = "Temperature")) +
    # This layer is drawn after temperature and uses the secondary top x-axis.
    ggplot2::geom_line(data = precipitation,
      ggplot2::aes(linetype = "Precipitation")) +
    ggplot2::scale_colour_manual(values = hypothesis_guild_colours,
      breaks = c("All mammals", "Herbivores", "Carnivores"), name = "Guild") +
    ggplot2::scale_linewidth_manual(values = c("All mammals" = 1.15,
      Herbivores = .65, Carnivores = .65), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("All mammals" = 1,
      Herbivores = .68, Carnivores = .68), guide = "none") +
    ggplot2::scale_linetype_manual(values = c(Temperature = "solid",
      Precipitation = "dotted"), name = "Climate variable") +
    ggplot2::scale_x_continuous(
      breaks = seq(-1, 1, by = .25),
      labels = function(x) format(signif(x * temperature_max, 3), trim = TRUE),
      sec.axis = ggplot2::sec_axis(~ . * precipitation_max,
        name = "Precipitation (mm/day)",
        breaks = seq(-1, 1, by = .25) * precipitation_max)) +
    ggplot2::scale_y_continuous(limits = c(.7, 1), breaks = seq(.7, 1, .05),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::facet_grid(continent ~ metric, scales = "free_x") +
    ggplot2::labs(x = "Temperature (deg C)",
      y = "Predicted Simpson turnover") + figpub_theme()
}

hypothesis_shape_summary <- function(curves) {
  keys <- c("continent", "guild", "variable", "metric", "predictor")
  group_id <- do.call(interaction, c(unname(curves[keys]),
    list(drop = TRUE, lex.order = TRUE)))
  groups <- split(curves, group_id)
  do.call(rbind, lapply(groups, function(d) {
    magnitude <- abs(d$x)
    low <- which.min(magnitude)
    high <- which(magnitude == max(magnitude))
    rho <- if (length(unique(magnitude)) > 1L &&
      length(unique(d$predicted_turnover)) > 1L) {
      stats::cor(magnitude, d$predicted_turnover, method = "spearman")
    } else NA_real_
    symmetry_error <- NA_real_
    if (unique(d$metric) == "signed_change") {
      negative <- d[d$x < 0, c("x", "predicted_turnover")]
      positive <- d[d$x > 0, c("x", "predicted_turnover")]
      negative <- negative[order(abs(negative$x)), ]
      positive <- positive[order(positive$x), ]
      if (nrow(negative) && nrow(positive)) {
        n <- min(nrow(negative), nrow(positive))
        symmetry_error <- max(abs(negative$predicted_turnover[seq_len(n)] -
          positive$predicted_turnover[seq_len(n)]))
      }
    }
    data.frame(continent = d$continent[1], guild = d$guild[1],
      variable = d$variable[1], metric = d$metric[1], predictor = d$predictor[1],
      predicted_at_zero_contrast = d$predicted_turnover[low],
      mean_predicted_at_max_contrast = mean(d$predicted_turnover[high]),
      max_minus_zero_prediction = mean(d$predicted_turnover[high]) -
        d$predicted_turnover[low],
      spearman_rho_vs_contrast_magnitude = rho,
      signed_change_symmetry_max_abs_difference = symmetry_error)
  }))
}

render_hypothesis_test <- function(cfg = bristol_config(),
    out_dir = "fig_hypothesis_test/output", overwrite = FALSE, n = 151L) {
  for (package in c("gdm", "ggplot2", "patchwork"))
    if (!requireNamespace(package, quietly = TRUE))
      stop("Required package: ", package, call. = FALSE)
  if (length(n) != 1L || !is.finite(n) || n < 3L)
    stop("n must be at least 3", call. = FALSE)
  bristol_require_outputs(cfg, "gdm_catalog.rds", "fit_bristol_gdms(cfg)")
  models <- readRDS(file.path(cfg$out, "gdm_catalog.rds"))
  variables <- c(temperature = "Temperature", precipitation = "Precipitation")
  metrics <- c("mean", "signed_change", "seasonal_change")
  planned_figures <- unlist(lapply(names(variables), function(v) {
    file.path(out_dir, paste0("h3_", v, "_contrast_hypothesis_test",
      c(".png", ".pdf")))
  }), use.names = FALSE)
  planned <- c(planned_figures,
    file.path(out_dir, c("h3_temperature_precipitation_contrast_hypothesis_test.png",
      "h3_temperature_precipitation_contrast_hypothesis_test.pdf",
      "h3_climate_pca_contrast_hypothesis_test.png",
      "h3_climate_pca_contrast_hypothesis_test.pdf",
      "h3_contrast_hypothesis_curves.csv",
      "h3_contrast_hypothesis_background.csv",
      "h3_contrast_hypothesis_shape_summary.csv",
      "h3_pca_contrast_hypothesis_curves.csv",
      "h3_pca_contrast_hypothesis_shape_summary.csv")))
  if (!overwrite && any(file.exists(planned)))
    stop("Output exists; set overwrite=TRUE or choose another out_dir.", call. = FALSE)

  rows <- pca_rows <- references <- list()
  for (region in c("europe", "north_america")) {
    continent <- if (region == "europe") "Europe" else "North America"
    reference_id <- paste("h3", region, "all", "combined", "geoTRUE", sep = "__")
    if (is.null(models[[reference_id]]) || !file.exists(models[[reference_id]]))
      stop("Missing H3 all-mammal geography-on background model: ", reference_id,
        call. = FALSE)
    reference <- readRDS(models[[reference_id]])
    background <- hypothesis_prediction_background(reference)
    for (guild in c("all", "herb", "carn")) {
      id <- paste("h3", region, guild, "combined", "geoTRUE", sep = "__")
      if (is.null(models[[id]]) || !file.exists(models[[id]]))
        stop("Missing H3 combined geographic model: ", id, call. = FALSE)
      object <- readRDS(models[[id]])
      if (!identical(object$variant, "combined") || !isTRUE(object$geo))
        stop("Expected a combined H3 geography-on model: ", id, call. = FALSE)
      for (variable in names(variables)) for (metric in metrics) {
        field <- paste(variable, metric, sep = "_")
        curve <- hypothesis_prediction_curve(object, background, field, metric, as.integer(n))
        rows[[length(rows) + 1L]] <- cbind(data.frame(continent, guild, variable,
          metric, predictor = field, model_id = id), curve)
        references[[length(references) + 1L]] <- data.frame(continent, guild,
          predictor = field, predictor_in_model = unique(curve$predictor_in_model),
          predictor_min = unique(curve$predictor_min), predictor_max = unique(curve$predictor_max),
          geographic_km = unique(curve$geographic_km), temporal_ma = unique(curve$temporal_ma))
      }

      pca_id <- paste("h3", region, guild, "hypothesis_pca", "geoTRUE", sep = "__")
      if (is.null(models[[pca_id]]) || !file.exists(models[[pca_id]]))
        stop("Missing H3 climate-PCA geography-on model: ", pca_id, call. = FALSE)
      pca_object <- readRDS(models[[pca_id]])
      if (!identical(pca_object$variant, "hypothesis_pca") || !isTRUE(pca_object$geo))
        stop("Expected an H3 climate-PCA geography-on model: ", pca_id, call. = FALSE)
      axis_fields <- sub("^s1[.]", "", grep(
        "^s1[.](mean|signed_change|seasonal_change)_PC1$",
        names(pca_object$gdm_table), value = TRUE))
      for (field in axis_fields) {
        metric <- sub("_PC[0-9]+$", "", field)
        pc <- sub("^.*_", "", field)
        pca_interpretation <- hypothesis_pca_extremes(cfg, region, guild, field, metric)
        curve <- hypothesis_prediction_curve(pca_object, background, field,
          metric, as.integer(n), orientation = pca_interpretation$orientation)
        pca_rows[[length(pca_rows) + 1L]] <- cbind(data.frame(continent, guild,
          variable = "climate_pca", metric, predictor = field, axis = pc,
          extreme_label = pca_interpretation$label,
          model_id = pca_id), curve)
      }
    }
  }
  curves <- do.call(rbind, rows)
  pca_curves <- do.call(rbind, pca_rows)
  background_table <- do.call(rbind, references)
  shape_summary <- hypothesis_shape_summary(curves)
  pca_shape_summary <- hypothesis_shape_summary(pca_curves)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(curves, file.path(out_dir, "h3_contrast_hypothesis_curves.csv"), row.names = FALSE)
  utils::write.csv(background_table, file.path(out_dir,
    "h3_contrast_hypothesis_background.csv"), row.names = FALSE)
  utils::write.csv(shape_summary, file.path(out_dir,
    "h3_contrast_hypothesis_shape_summary.csv"), row.names = FALSE)
  utils::write.csv(pca_curves, file.path(out_dir,
    "h3_pca_contrast_hypothesis_curves.csv"), row.names = FALSE)
  utils::write.csv(pca_shape_summary, file.path(out_dir,
    "h3_pca_contrast_hypothesis_shape_summary.csv"), row.names = FALSE)

  figures <- list()
  for (variable in names(variables)) {
    panels <- list()
    for (continent in c("Europe", "North America")) for (metric in metrics)
      panels[[length(panels) + 1L]] <- hypothesis_test_panel(curves[
        curves$variable == variable, ], metric, continent, variable)
    figure <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
      figpub_annotation(paste("H3 GDM test of conceptual", tolower(variables[[variable]]),
        "contrast hypotheses"),
        paste("H3 combined models with geography included; the y-axis shows 0.70-1.00 on the 0-1 turnover scale.",
          "Each panel varies only the named between-assemblage climate contrast.",
          "Other climate predictors, sampling, and stage duration are equalized; geographic and temporal separation are held at continent-specific all-mammal medians.",
          "These are conditional model predictions, not causal effects or uncertainty intervals.")) &
      ggplot2::theme(legend.position = "bottom")
    stem <- paste0("h3_", variable, "_contrast_hypothesis_test")
    figpub_save(out_dir, stem, figure, height = 6.5, overwrite = overwrite)
    figures[[variable]] <- figure
  }
  both_figure <- dual_axis_hypothesis_plot(curves) +
    figpub_annotation("H3 GDM test of climate-contrast hypotheses",
      paste("Bottom axes show temperature contrast; top axes show precipitation contrast.",
        "Curves are scaled to each variable's maximum contrast across the fitted groups.",
        "Precipitation is dotted and drawn over temperature; the y-axis is shared because both models predict Simpson turnover.",
        "Displayed y-range is 0.70-1.00 on the 0-1 response scale."))
  figpub_save(out_dir, "h3_temperature_precipitation_contrast_hypothesis_test",
    both_figure, height = 7.2, width = 15, overwrite = overwrite)
  figures$temperature_precipitation <- both_figure

  pca_panels <- list()
  for (continent in c("Europe", "North America")) for (metric in metrics)
    pca_panels[[length(pca_panels) + 1L]] <- hypothesis_pca_panel(
      pca_curves, metric, continent)
  pca_figure <- patchwork::wrap_plots(pca_panels, ncol = 3, guides = "collect") +
    figpub_annotation("H3 GDM test using climate PCA axes",
      paste("Only PC1 is shown, with all three guild curves in each continent-by-hypothesis panel.",
        "PC1 signs are oriented toward higher temperature; loading-based endpoint descriptions are available in the curve CSV for manual annotation.",
        "PCA inputs are standardized temperature and precipitation metrics, so annotations describe relative climate directions rather than native-unit changes.",
        "Each group's other climate axes are equalized; PC1 fits remain group-specific.",
        "Predictions are conditional model associations, not causal effects or uncertainty intervals.")) &
    ggplot2::theme(legend.position = "bottom")
  figpub_save(out_dir, "h3_climate_pca_contrast_hypothesis_test",
    pca_figure, height = 6.5, width = 7.5, overwrite = overwrite)
  figures$climate_pca <- pca_figure
  invisible(list(figures = figures, curves = curves, background = background_table,
    shape_summary = shape_summary, pca_curves = pca_curves,
    pca_shape_summary = pca_shape_summary))
}
