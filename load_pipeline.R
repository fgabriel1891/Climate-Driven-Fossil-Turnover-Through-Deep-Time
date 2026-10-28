# Definitions only; run from the repository root.
if (!file.exists("config/stages.csv") || !file.exists("R/read_climate.R")) {
  stop("Set the working directory to the repository root first.")
}
source("config/config.R")
for (file in c("read_climate.R", "01_prepare_fossils.R", "02_prepare_climate.R",
               "03_prepare_assemblages.R", "04_fit_gdms.R", "05_summarize.R",
               "manuscript_tables.R", "publication_figures.R")) {
  source(file.path("R", file))
}

preflight_bristol <- function(cfg = bristol_config()) {
  files <- c(cfg$fossil_file, cfg$region_file, cfg$stages_file,
             bristol_pipeline_sources(cfg)$path)
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing input files: ", paste(missing, collapse = ", "))
  invisible(files)
}

run_bristol_pipeline <- function(cfg = bristol_config(), publication = TRUE, predictions = FALSE) {
  preflight_bristol(cfg)
  prepare_bristol_fossils(cfg)
  prepare_bristol_climate(cfg)
  prepare_bristol_assemblages(cfg)
  compute_bristol_turnover(cfg)
  fit_bristol_gdms(cfg)
  summarize_bristol_gdms(cfg)
  manuscript_tables(cfg)
  if (publication) publication_figures(cfg)
  if (predictions) predict_bristol_turnover(cfg)
  invisible(cfg)
}
