# Fork of numerical_analyses/dissimilarity_modelling.R and src/gdm_functions.R.
# Standard GDM site-pair predictors; all unordered pairs within each group.
bristol_climate_names <- function() {
  unlist(lapply(c("mean", "signed_change", "seasonal_change"), function(metric) {
    paste(c("temperature", "precipitation"), metric, sep = "_")
  }), use.names = FALSE)
}

bristol_hypothesis_pca <- function(sites, cfg) {
  scores <- data.frame(row.names = seq_len(nrow(sites)))
  loadings <- list()
  models <- list()
  for (hypothesis in c("mean", "signed_change", "seasonal_change")) {
    variables <- paste(c("temperature", "precipitation"), hypothesis, sep = "_")
    x <- sites[, variables, drop = FALSE]
    varying <- vapply(x, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0, logical(1))
    if (!any(varying)) next
    if (sum(varying) == 1L) {
      name <- paste0(hypothesis, "_", if (varying[1]) "temperature_only" else "precipitation_only")
      scores[[name]] <- as.numeric(scale(x[[which(varying)]]))
      loadings[[hypothesis]] <- data.frame(hypothesis = hypothesis, axis = name,
                                          variable = variables[varying], loading = 1,
                                          explained_fraction = 1, retained = TRUE,
                                          interpretation = "Only one nonconstant variable; standardized, no PCA")
      models[[hypothesis]] <- list(variable = variables[varying], center = mean(x[[which(varying)]]), scale = stats::sd(x[[which(varying)]]))
      next
    }
    pc <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    # Orient each axis so its largest absolute loading is positive.
    for (j in seq_len(ncol(pc$rotation))) {
      dominant <- which.max(abs(pc$rotation[, j]))
      if (pc$rotation[dominant, j] < 0) {
        pc$rotation[, j] <- -pc$rotation[, j]
        pc$x[, j] <- -pc$x[, j]
      }
    }
    fraction <- pc$sdev^2/sum(pc$sdev^2)
    nretain <- which(cumsum(fraction) >= cfg$pca_variance_target)[1]
    rows <- list()
    for (j in seq_along(fraction)) {
      axis <- paste0(hypothesis, "_PC", j)
      w <- pc$rotation[, j]
      interpretable <- max(abs(w)) >= cfg$pca_min_dominant_loading
      retained <- j <= nretain && interpretable
      label <- if (abs(abs(w[1]) - abs(w[2])) <= 0.2) {
        if (prod(w) > 0) "coupled temperature/precipitation axis" else "temperature/precipitation contrast"
      } else if (abs(w[1]) > abs(w[2])) "temperature-dominated axis" else "precipitation-dominated axis"
      if (retained) scores[[axis]] <- pc$x[, j]
      rows[[j]] <- data.frame(hypothesis = hypothesis, axis = axis, variable = variables,
                              loading = as.numeric(w), explained_fraction = fraction[j],
                              retained = retained, interpretation = label)
    }
    loadings[[hypothesis]] <- do.call(rbind, rows)
    models[[hypothesis]] <- pc
  }
  if (!ncol(scores)) stop("No interpretable nonconstant PCA predictors")
  list(scores = scores, loadings = do.call(rbind, loadings), models = models)
}

bristol_pair_table <- function(sites, cfg) {
  n <- nrow(sites)
  count <- n*(n - 1)/2
  if (count < cfg$min_pairs) stop("Too few complete assemblage pairs: ", count)
  if (count > cfg$max_pairs) stop("Pair count ", count, " exceeds cfg$max_pairs; increase explicitly after checking memory. No pairs were sampled.")
  taxa <- sort(unique(unlist(sites$genera)))
  community <- matrix(0L, nrow = n, ncol = length(taxa), dimnames = list(sites$site, taxa))
  for (i in seq_len(n)) community[i, match(sites$genera[[i]], taxa)] <- 1L
  turnover <- as.matrix(betapart::beta.pair(community)$beta.sim)
  ij <- which(upper.tri(turnover), arr.ind = TRUE)
  i <- ij[, 1]
  j <- ij[, 2]
  richness_weight <- pmin(sites$richness[i], sites$richness[j])
  duration_weight <- 1/((sites$duration_ma[i] + sites$duration_ma[j])/2)
  weights <- switch(cfg$weight_method,
                    richness_per_mean_duration = richness_weight * duration_weight,
                    richness = richness_weight,
                    equal = rep(1, length(i)), stop("Unknown weight method"))
  weights <- weights/mean(weights)
  stopifnot(all(is.finite(weights)), all(weights > 0), all(is.finite(turnover[cbind(i, j)])))
  table <- data.frame(distance = turnover[cbind(i, j)], weights = weights,
                      s1.xCoord = sites$x_km[i], s1.yCoord = sites$y_km[i],
                      s2.xCoord = sites$x_km[j], s2.yCoord = sites$y_km[j])
  metadata <- data.frame(site1 = sites$site[i], site2 = sites$site[j],
                         stage1 = sites$stage[i], stage2 = sites$stage[j],
                         richness_weight = richness_weight, duration_weight = duration_weight,
                         normalized_weight = weights,
                         same_cell = sites$cell_id[i] == sites$cell_id[j],
                         same_stage = sites$stage_id[i] == sites$stage_id[j])
  list(table = table, metadata = metadata, i = i, j = j)
}

fit_bristol_gdms <- function(cfg = bristol_config()) {
  bristol_packages(c("gdm", "betapart"))
  bristol_require_outputs(cfg, "assemblage_catalog.rds", "prepare_bristol_assemblages(cfg)")
  if (cfg$model_spec != "standard_gdm") stop("Only original-concept standard GDM is supported")
  output <- bristol_output(cfg, "gdm_catalog.rds")
  assemblies <- readRDS(file.path(cfg$out, "assemblage_catalog.rds"))
  catalog <- status <- list()
  fields <- bristol_climate_names()
  message("[4/5] Fitting all-pair GDMs: combined climate and hypothesis-block PCA")
  for (key in names(assemblies)) {
    obj <- readRDS(assemblies[[key]])
    s <- obj$sites
    required <- c(fields, "time_ma", "sampling", "duration_ma", "x_km", "y_km")
    complete <- apply(as.matrix(s[, required]), 1, function(v) all(is.finite(v)))
    s <- s[complete, ]
    result <- tryCatch({
      if (anyDuplicated(s$site)) stop("Duplicate assemblage IDs")
      pairs <- bristol_pair_table(s, cfg)
      saveRDS(pairs$metadata, file.path(cfg$out, paste0("pair_metadata_", key, ".rds")))
      utils::write.csv(stats::cor(s[, fields]), file.path(cfg$out, paste0("climate_correlations_", key, ".csv")))
      pca <- tryCatch(bristol_hypothesis_pca(s, cfg), error = function(e) e)
      variants <- list(combined = s[, fields, drop = FALSE])
      if (!inherits(pca, "error")) {
        variants$hypothesis_pca <- pca$scores
        utils::write.csv(pca$loadings, file.path(cfg$out, paste0("pca_loadings_", key, ".csv")), row.names = FALSE)
        saveRDS(pca, file.path(cfg$out, paste0("pca_", key, ".rds")))
      } else {
        status[[paste0(key, "__pca")]] <- data.frame(group = key, variant = "hypothesis_pca", geo = NA,
          sites = nrow(s), pairs = nrow(pairs$table), ok = FALSE, explained = NA_real_, error = conditionMessage(pca))
      }
      for (variant in names(variants)) {
        predictors <- cbind(variants[[variant]], s[, c("time_ma", "sampling", "duration_ma")])
        varying <- vapply(predictors, function(v) length(unique(v)) >= 3L, logical(1))
        dropped <- names(predictors)[!varying]
        predictors <- predictors[, varying, drop = FALSE]
        if (!ncol(predictors)) stop("No variable predictors")
        left <- predictors[pairs$i, , drop = FALSE]
        right <- predictors[pairs$j, , drop = FALSE]
        names(left) <- paste0("s1.", names(predictors))
        names(right) <- paste0("s2.", names(predictors))
        tab <- cbind(pairs$table, left, right)
        rownames(tab) <- NULL
        for (geo in cfg$geo_models) {
          id <- paste(key, variant, paste0("geo", geo), sep = "__")
          message("  ", id)
          fit <- tryCatch(gdm::gdm(tab, geo = geo), error = function(e) e)
          ok <- !inherits(fit, "error") && !is.null(fit)
          status[[id]] <- data.frame(group = key, variant = variant, geo = geo, sites = nrow(s),
            pairs = nrow(tab), ok = ok, explained = if (ok) fit$explained else NA_real_,
            error = if (inherits(fit, "error")) conditionMessage(fit) else if (is.null(fit)) "GDM returned NULL" else "")
          if (ok) {
            path <- file.path(cfg$out, paste0("gdm_", id, ".rds"))
            saveRDS(list(model = fit, gdm_table = tab, site_data = s[, setdiff(names(s), "genera")],
                         variant = variant, geo = geo, group = key, dropped_predictors = dropped,
                         pair_metadata_file = file.path(cfg$out, paste0("pair_metadata_", key, ".rds")), config = cfg), path)
            catalog[[id]] <- path
          }
        }
      }
      TRUE
    }, error = function(e) {
      status[[key]] <<- data.frame(group = key, variant = "preparation", geo = NA, sites = nrow(s),
        pairs = nrow(s)*(nrow(s) - 1)/2, ok = FALSE, explained = NA_real_, error = conditionMessage(e))
      FALSE
    })
    utils::write.csv(do.call(rbind, status), file.path(cfg$out, "gdm_status.csv"), row.names = FALSE)
  }
  saveRDS(catalog, output)
  invisible(catalog)
}
