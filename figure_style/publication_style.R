# Shared, print-sized visual conventions for the dedicated Bristol figures.
# All figures are designed for a 7.5-inch, two-column manuscript width.

figpub_colours <- list(
  region = c(Europe = "#0072B2", `North America` = "#D55E00"),
  climate = c(Temperature = "#D55E00", Precipitation = "#0072B2", Both = "#4D4D4D"),
  guild = c(`All mammals` = "#4D4D4D", Herbivores = "#009E73", Carnivores = "#D55E00"),
  diverging = c(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00")
)

figpub_theme <- function() {
  ggplot2::theme_bw(base_size = 10.5, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = .25),
      axis.title = ggplot2::element_text(size = 9.5),
      axis.text = ggplot2::element_text(size = 8.5, colour = "grey15"),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.text = ggplot2::element_text(size = 8.5),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 9.5, face = "bold"),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 8.5),
      plot.margin = ggplot2::margin(6, 6, 6, 6))
}

figpub_map_theme <- function() {
  ggplot2::theme_void(base_size = 10.5, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 8.5),
      axis.title = ggplot2::element_text(size = 9.5),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.text = ggplot2::element_text(size = 8.5),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(6, 6, 6, 6))
}

figpub_annotation <- function(title, caption = NULL) {
  patchwork::plot_annotation(
    tag_levels = "A", title = title, caption = caption,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(family = "sans", size = 13.5, face = "bold"),
      plot.caption = ggplot2::element_text(family = "sans", size = 8.5),
      plot.tag = ggplot2::element_text(family = "sans", size = 13, face = "bold"),
      plot.tag.position = c(0.01, 0.99)))
}

figpub_save <- function(out_dir, stem, plot, height, width = 7.5, dpi = 400,
                        overwrite = FALSE) {
  paths <- file.path(out_dir, paste0(stem, c(".png", ".pdf")))
  if (!overwrite && any(file.exists(paths)))
    stop("Figure output already exists: ", stem, "; set overwrite=TRUE")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (path in paths)
    ggplot2::ggsave(path, plot, width = width, height = height,
      units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(paths)
}
