# Run from the repository root.
bristol_config <- function() list(
  root = ".",
  out = "outputs/data",
  tables = "outputs/tables",
  publication_figs = "outputs/publication_figs",
  fossil_file = "DATA/NOW_master_with_diet.csv",
  region_file = "DATA/regions_gridded/Data_RegionsGridded.shp",
  raw_dir = "DATA/climate",
  regions = list(europe = c(-10, 60, 20, 84), north_america = c(-180, -50, 20, 84)),
  stages_file = "config/stages.csv",
  spacing = c(h3 = 103.6, h2 = 274.09, h1 = 725.17, h0 = 1918.61),
  min_genera = 3L,
  min_pairs = 30L,
  max_pairs = 2000000L,
  baseline = c(stage = "Chattian", experiment = "tfgsV"),
  # A transparent analytical weighting choice, not an inverse-variance weight.
  weight_method = "richness_per_mean_duration",
  # Fixed geographic coordinates reproduce the original analysis convention.
  # See README before interpreting this as palaeogeographically aligned.
  coordinate_mode = "present_day",
  pca_variance_target = 0.90,
  pca_min_dominant_loading = 0.70,
  geo_models = c(TRUE, FALSE),
  figure_resolution = "h3",
  figure_formats = c("png", "pdf"),
  figure_dpi = 300,
  random_seed = 42L,
  prediction_resolution = "h1",
  prediction_guild = "herb",
  prediction_geo = TRUE,
  model_spec = "standard_gdm",
  overwrite = FALSE
)

read_bristol_stages <- function(cfg) {
  s <- utils::read.csv(cfg$stages_file, stringsAsFactors = FALSE)
  stopifnot(nrow(s) == 8L, identical(s$stage_id, 1:8),
            all(s$older_ma > s$younger_ma),
            all(s$younger_ma[-8] == s$older_ma[-1]),
            !anyDuplicated(s$stage), !anyDuplicated(s$experiment))
  s$duration_ma <- s$older_ma - s$younger_ma
  s$stage_midpoint_ma <- (s$older_ma + s$younger_ma)/2
  s
}

bristol_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
}

bristol_output <- function(cfg, name) {
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cfg$out, name)
  if (file.exists(path) && !cfg$overwrite) stop("Output exists: ", path, "; set cfg$overwrite=TRUE explicitly to replace")
  path
}

bristol_require_outputs <- function(cfg, names, prerequisite) {
  paths <- file.path(cfg$out, names)
  if (any(!file.exists(paths))) {
    stop("Missing prerequisite output: ", paste(paths[!file.exists(paths)], collapse = ", "),
         ". Finish ", prerequisite, " successfully before continuing. If at Browse[1]>, enter Q first.",
         call. = FALSE)
  }
  invisible(paths)
}
