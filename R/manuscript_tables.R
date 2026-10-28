# Stage-aware counterparts of scripts/03_compute_turnover.R,
# 05_summary_statistics.R and summarize_fossil_turnover.R in the parent project.
compute_bristol_turnover <- function(cfg = bristol_config()) {
  bristol_packages("betapart")
  bristol_require_outputs(cfg, "assemblage_catalog.rds", "prepare_bristol_assemblages(cfg)")
  output <- bristol_output(cfg, "turnover_catalog.rds")
  catalog <- readRDS(file.path(cfg$out, "assemblage_catalog.rds"))
  result <- list()
  message("[Turnover] Saving full beta.sim matrices independently of climate completeness")
  for (key in names(catalog)) {
    obj <- readRDS(catalog[[key]])
    s <- obj$sites
    if (nrow(s) < 2L) next
    if (choose(nrow(s), 2) > cfg$max_pairs) stop("Turnover group exceeds max_pairs: ", key)
    taxa <- sort(unique(unlist(s$genera)))
    x <- matrix(0L, nrow(s), length(taxa), dimnames = list(s$site, taxa))
    for (i in seq_len(nrow(s))) x[i, match(s$genera[[i]], taxa)] <- 1L
    distance <- betapart::beta.pair(x)$beta.sim
    path <- file.path(cfg$out, paste0("turnover_", key, ".rds"))
    saveRDS(list(distance = distance, sites = s, continent = obj$continent,
                 guild = obj$guild, resolution = obj$resolution), path)
    result[[key]] <- path
  }
  saveRDS(result, output)
  invisible(result)
}

manuscript_write_table <- function(x, stem, cfg) {
  dir.create(cfg$tables, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file.path(cfg$tables, paste0(stem, ".csv")), row.names = FALSE, na = "NA")
  invisible(x)
}

manuscript_tables <- function(cfg = bristol_config()) {
  bristol_packages("dplyr")
  bristol_require_outputs(cfg, c("stage_fossils.rds", "turnover_catalog.rds"), "fossil and turnover preparation")
  output <- bristol_output(cfg, "manuscript_tables.rds")
  fossil <- readRDS(file.path(cfg$out, "stage_fossils.rds"))
  x <- fossil$occurrences
  s <- fossil$stages
  message("[Tables] Fossils, uncertainty, richness, turnover and model inventory")
  tables <- list()
  tables$summary_continent <- x |> dplyr::group_by(continent) |>
    dplyr::summarise(n_records = dplyr::n(), n_genera = dplyr::n_distinct(GENUS),
                     n_localities = dplyr::n_distinct(LIDNUM, na.rm = TRUE), .groups = "drop")
  tables$summary_guild <- x |> dplyr::group_by(continent, trophic_guild) |>
    dplyr::summarise(n_records = dplyr::n(), n_genera = dplyr::n_distinct(GENUS),
                     n_localities = dplyr::n_distinct(LIDNUM, na.rm = TRUE), .groups = "drop")
  localities <- x |> dplyr::filter(!is.na(LIDNUM)) |>
    dplyr::group_by(continent, stage_id, LIDNUM) |>
    dplyr::summarise(n_genera = dplyr::n_distinct(GENUS), .groups = "drop")
  tables$locality_richness_by_stage <- localities |> dplyr::group_by(continent, stage_id) |>
    dplyr::summarise(n_localities = dplyr::n(), mean_gen_per_loc = mean(n_genera),
                     sd_gen_per_loc = stats::sd(n_genera), median_gen_per_loc = stats::median(n_genera), .groups = "drop")
  tables$sampling_density <- localities |> dplyr::group_by(continent) |>
    dplyr::summarise(n_locality_stages = dplyr::n(), mean_gen_per_loc_stage = mean(n_genera), .groups = "drop")
  tables$Table_S1 <- dplyr::left_join(tables$summary_guild, tables$summary_continent,
                                     by = "continent", suffix = c("_guild", "_continent"))
  richness <- species_excluded <- list()
  for (continent in c("Europe", "North America")) for (guild in c("all", "herb", "carn")) for (stage in 1:8) {
    z <- x[x$continent == continent & x$stage_id == stage, ]
    if (guild != "all") z <- z[!is.na(z$trophic_guild) & z$trophic_guild == if (guild == "herb") "Herbivore" else "Carnivore", ]
    species <- if ("SPECIES" %in% names(z)) trimws(as.character(z$SPECIES)) else rep(NA_character_, nrow(z))
    identified <- !is.na(species) & nzchar(species) & !grepl("(^|[ .])(sp|spp|indet|cf|aff)([ .]|$)|[?]", species, ignore.case = TRUE)
    species_excluded[[length(species_excluded) + 1L]] <- data.frame(continent, guild, stage_id = stage, n_excluded = sum(!identified))
    richness[[length(richness) + 1L]] <- data.frame(continent, guild, stage_id = stage, stage = s$stage[stage],
      n_records = nrow(z), n_genera = length(unique(z$GENUS)), n_species = length(unique(paste(z$GENUS[identified], species[identified]))),
      n_localities = length(unique(z$LIDNUM[!is.na(z$LIDNUM)])))
  }
  tables$genus_species_richness_by_stage <- do.call(rbind, richness)
  tables$excluded_species_labels <- do.call(rbind, species_excluded)
  catalog <- readRDS(file.path(cfg$out, "turnover_catalog.rds"))
  overview <- by_stage <- by_lag <- matrices <- list()
  describe <- function(v) data.frame(n_pairs = length(v), mean_turnover = if (length(v)) mean(v) else NA_real_,
    sd_turnover = if (length(v) > 1) stats::sd(v) else NA_real_, median_turnover = if (length(v)) stats::median(v) else NA_real_,
    q25 = if (length(v)) unname(stats::quantile(v, .25)) else NA_real_, q75 = if (length(v)) unname(stats::quantile(v, .75)) else NA_real_)
  for (key in names(catalog)) {
    obj <- readRDS(catalog[[key]])
    d <- as.matrix(obj$distance)
    ids <- which(upper.tri(d), arr.ind = TRUE)
    a <- ids[, 1]; b <- ids[, 2]
    sites <- obj$sites
    val <- d[ids]
    same_cell <- sites$cell_id[a] == sites$cell_id[b]
    same_stage <- sites$stage_id[a] == sites$stage_id[b]
    meta <- data.frame(group = key, resolution = obj$resolution, continent = obj$continent, guild = obj$guild)
    selectors <- list(all_pairs = rep(TRUE, length(val)), same_stage_spatial = same_stage,
                      same_cell_temporal = same_cell, different_cell_and_stage = !same_cell & !same_stage)
    for (type in names(selectors)) overview[[length(overview) + 1L]] <- cbind(meta, comparison = type, describe(val[selectors[[type]]]))
    for (stage in 1:8) by_stage[[length(by_stage) + 1L]] <- cbind(meta, stage_id = stage, stage = s$stage[stage],
      n_assemblages = sum(sites$stage_id == stage), describe(val[same_stage & sites$stage_id[a] == stage]))
    for (lag in 1:7) by_lag[[length(by_lag) + 1L]] <- cbind(meta, stage_lag = lag,
      describe(val[same_cell & abs(sites$stage_id[a] - sites$stage_id[b]) == lag]))
    # Symmetric 8x8 stage-aggregated matrix; each off-diagonal cell uses all
    # cross-stage assemblage pairs, not just matching geographic cells.
    for (i in 1:8) for (j in i:8) {
      keep <- (sites$stage_id[a] == i & sites$stage_id[b] == j) | (sites$stage_id[a] == j & sites$stage_id[b] == i)
      row <- cbind(meta, stage1 = i, stage2 = j, describe(val[keep]))
      matrices[[length(matrices) + 1L]] <- row
      if (i != j) { row$stage1 <- j; row$stage2 <- i; matrices[[length(matrices) + 1L]] <- row }
    }
  }
  tables$turnover_comparison_summary <- do.call(rbind, overview)
  tables$turnover_same_stage_spatial <- do.call(rbind, by_stage)
  tables$turnover_same_cell_by_stage_lag <- do.call(rbind, by_lag)
  tables$turnover_stage_matrix <- do.call(rbind, matrices)
  status_path <- file.path(cfg$out, "gdm_status.csv")
  if (file.exists(status_path)) {
    status <- utils::read.csv(status_path)
    keys <- strsplit(status$group, "__", fixed = TRUE)
    status$resolution <- vapply(keys, `[`, character(1), 1)
    status$continent <- ifelse(vapply(keys, `[`, character(1), 2) == "europe", "Europe", "North America")
    status$guild <- vapply(keys, `[`, character(1), 3)
    tables$Table_S2_GDM_summary <- status
  }
  for (name in names(tables)) manuscript_write_table(tables[[name]], name, cfg)
  for (name in c("assignment_uncertainty.csv", "occurrence_filter_counts.csv", "assemblage_coverage.csv", "climate_units_audit.csv")) {
    path <- file.path(cfg$out, name)
    if (file.exists(path)) file.copy(path, file.path(cfg$tables, name), overwrite = TRUE)
  }
  saveRDS(tables, output)
  invisible(tables)
}
