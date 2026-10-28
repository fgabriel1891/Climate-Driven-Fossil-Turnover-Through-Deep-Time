# Fork of original status/ispline figures; one independently labelled model
# per figure. No output is produced until this function is called.
summarize_bristol_gdms <- function(cfg = bristol_config()) {
  bristol_packages(c("gdm", "ggplot2"))
  bristol_require_outputs(cfg, "gdm_catalog.rds", "fit_bristol_gdms(cfg)")
  marker <- bristol_output(cfg, "summary_complete.rds")
  catalog <- readRDS(file.path(cfg$out, "gdm_catalog.rds"))
  figures <- file.path(cfg$out, "figures")
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  message("[5/5] Exporting GDM spline curves and fit summaries")
  for (id in names(catalog)) {
    obj <- readRDS(catalog[[id]])
    curves <- gdm::isplineExtract(obj$model)
    labels <- colnames(curves$x)
    if (is.null(labels)) labels <- obj$model$predictors
    if (length(labels) != ncol(curves$x)) stop("Unexpected spline column count")
    long <- do.call(rbind, lapply(seq_len(ncol(curves$x)), function(j) {
      data.frame(predictor = labels[j], x = curves$x[, j], transformed = curves$y[, j])
    }))
    utils::write.csv(long, file.path(figures, paste0(id, "_splines.csv")), row.names = FALSE)
    plot <- ggplot2::ggplot(long, ggplot2::aes(x, transformed)) + ggplot2::geom_line() +
      ggplot2::facet_wrap(~ predictor, scales = "free") + ggplot2::theme_bw() +
      ggplot2::labs(title = id, x = "Predictor (native units or standardized PC score)", y = "GDM transformed environmental distance")
    ggplot2::ggsave(file.path(figures, paste0(id, "_splines.png")), plot, width = 12, height = 8, dpi = 200)
    predictions <- data.frame(observed = obj$model$observed, predicted = obj$model$predicted)
    utils::write.csv(predictions, file.path(figures, paste0(id, "_fitted.csv")), row.names = FALSE)
  }
  saveRDS(list(models = length(catalog), completed = Sys.time(), config = cfg,
               session = utils::sessionInfo()), marker)
  invisible(catalog)
}
