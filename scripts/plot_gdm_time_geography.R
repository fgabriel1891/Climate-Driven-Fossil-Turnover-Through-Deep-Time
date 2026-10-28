# Four-panel temporal and geographic GDM spline figures from saved model output.
# Source from the repository root; rendering is explicit.

source("figure_style/publication_style.R")

fig_time_geography_data <- function(cfg) {
  input_file <- file.path(cfg$tables, "all_gdm_spline_coordinates.csv")
  if (!file.exists(input_file)) {
    stop("Missing saved GDM spline coordinates: ", input_file, call. = FALSE)
  }
  splines <- utils::read.csv(input_file, stringsAsFactors = FALSE)
  required <- c("model", "resolution", "continent", "guild", "variant",
                "geo", "predictor", "x", "y")
  if (!all(required %in% names(splines))) {
    stop("Unexpected GDM spline columns in ", input_file, call. = FALSE)
  }
  keep <- splines$resolution %in% c("h3", "h2", "h1", "h0") &
    splines$variant == "combined" &
    as.logical(splines$geo) &
    splines$predictor %in% c("time_ma", "Geographic") &
    splines$guild %in% c("all", "herb", "carn") &
    splines$continent %in% c("Europe", "North America")
  x <- splines[keep, ]
  expected <- expand.grid(
    resolution = c("h3", "h2", "h1", "h0"),
    predictor = c("time_ma", "Geographic"),
    continent = c("Europe", "North America"),
    guild = c("all", "herb", "carn"),
    stringsAsFactors = FALSE
  )
  observed <- unique(x[c("resolution", "predictor", "continent", "guild")])
  if (nrow(observed) != nrow(expected) || anyNA(x$x) || anyNA(x$y)) {
    stop("Expected complete time/geography curves for every H resolution, continent, and guild.",
         call. = FALSE)
  }
  x$guild_label <- factor(
    x$guild,
    levels = c("all", "herb", "carn"),
    labels = c("All mammals", "Herbivores", "Carnivores")
  )
  x
}

fig_time_geography_panel <- function(data, resolution, predictor, continent) {
  d <- data[data$resolution == resolution & data$predictor == predictor &
              data$continent == continent, ]
  if (!nrow(d)) {
    stop("Missing spline data for ", resolution, " / ", predictor, " / ", continent,
         call. = FALSE)
  }

  # Keep axes matched between Europe and North America within each row.
  row_data <- data[data$resolution == resolution & data$predictor == predictor, ]
  xmax <- max(row_data$x, na.rm = TRUE) * 1.02
  ymax <- max(row_data$y, na.rm = TRUE) * 1.04
  x_label <- if (predictor == "time_ma") {
    "Stage midpoint age (Ma)"
  } else {
    "Geographic distance (km)"
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = x, y = y, colour = guild_label, group = model)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = figpub_colours$guild, name = NULL) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, ymax),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = continent,
      x = x_label,
      y = "GDM contribution"
    ) +
    figpub_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
}

render_fig_time_geography <- function(
    cfg = bristol_config(),
    out_dir = file.path("fig_time_geography", "output"),
    overwrite = FALSE) {
  for (package in c("ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required package: ", package, call. = FALSE)
    }
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }

  data <- fig_time_geography_data(cfg)
  resolutions <- c("h3", "h2", "h1", "h0")
  stems <- paste0("gdm_time_geography_", resolutions)
  planned <- unlist(lapply(stems, function(stem) {
    file.path(out_dir, paste0(stem, c(".png", ".pdf")))
  }), use.names = FALSE)
  if (!overwrite && any(file.exists(planned))) {
    stop("At least one time/geography figure already exists in ", out_dir,
         "; set overwrite=TRUE or choose another out_dir.", call. = FALSE)
  }

  figures <- files <- stats::setNames(vector("list", length(resolutions)), resolutions)
  for (resolution in resolutions) {
    # A-B: stage midpoint age; C-D: geographic distance.
    panels <- list(
      fig_time_geography_panel(data, resolution, "time_ma", "Europe"),
      fig_time_geography_panel(data, resolution, "time_ma", "North America"),
      fig_time_geography_panel(data, resolution, "Geographic", "Europe"),
      fig_time_geography_panel(data, resolution, "Geographic", "North America")
    )
    figure <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
      figpub_annotation(
        paste("Geological age and geographic distance in the GDM -", toupper(resolution)),
        paste(
          "Combined models with geographic distance included.",
          "Curves show fitted GDM contributions conditional on climate, sampling, duration, and the other spatial or temporal predictor."
        )
      ) &
      ggplot2::theme(legend.position = "bottom")
    figures[[resolution]] <- figure
    files[[resolution]] <- figpub_save(
      out_dir,
      stems[match(resolution, resolutions)],
      figure,
      height = 6.5,
      overwrite = overwrite
    )
  }
  invisible(list(figures = figures, data = data, files = files))
}
