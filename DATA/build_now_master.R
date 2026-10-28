# Build the analysis input from the local NOW occurrence extract.
# Sourcing this file defines functions only; call build_now_master() explicitly.

build_now_master <- function(
    source_file = "DATA/NOW_filtered_by_Shan_Oct2019.csv",
    out_dir = "DATA", overwrite = FALSE) {
  if (!file.exists(source_file)) stop("Missing NOW source: ", source_file, call. = FALSE)
  if (!dir.exists(out_dir)) stop("Missing output directory: ", out_dir, call. = FALSE)
  paths <- file.path(out_dir, c("NOW_master_with_diet.csv",
    "NOW_genus_diet_audit.csv", "NOW_master_manifest.csv"))
  if (!overwrite && any(file.exists(paths)))
    stop("NOW master outputs already exist; set overwrite = TRUE", call. = FALSE)
  if (normalizePath(source_file) %in% normalizePath(paths, mustWork = FALSE))
    stop("The NOW source cannot also be a master output", call. = FALSE)

  now <- utils::read.csv(source_file, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "\\N"), check.names = FALSE,
    encoding = "UTF-8")
  required <- c("LIDNUM", "ORDER", "FAMILY", "GENUS", "SPECIES",
    "MAX_AGE", "MIN_AGE", "LAT", "LONG", "DIET_1", "DIET_2", "DIET_3")
  missing <- setdiff(required, names(now))
  if (length(missing)) stop("NOW source lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(now)) stop("NOW source contains no occurrences", call. = FALSE)
  if (basename(source_file) == "NOW_filtered_by_Shan_Oct2019.csv" &&
      identical(unname(tools::md5sum(source_file)), "f9000bd3c43e20bbb1fdd985cd7e18a5") &&
      nrow(now) != 18349L)
    stop("The archived NOW extract did not read completely", call. = FALSE)

  clean <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x %in% c("", "NA", "\\N")] <- NA_character_
    x
  }
  now$DIET_2 <- tolower(clean(now$DIET_2))
  allowed <- c("herbivore", "carnivore", "omnivore", "invert", "piscivore")
  unknown <- setdiff(unique(stats::na.omit(now$DIET_2)), allowed)
  if (length(unknown))
    stop("Unreviewed NOW DIET_2 categories: ", paste(unknown, collapse = ", "), call. = FALSE)

  genus <- clean(now$GENUS)
  invalid <- is.na(genus) | tolower(genus) %in% c("indet", "indet.", "gen.")
  taxon_key <- paste(tolower(clean(now$ORDER)), tolower(clean(now$FAMILY)),
    tolower(genus), sep = "|")
  groups <- split(seq_len(nrow(now))[!invalid], taxon_key[!invalid])
  audit <- lapply(groups, function(i) {
    diets <- sort(unique(stats::na.omit(now$DIET_2[i])))
    status <- if (!length(diets)) "missing" else if (length(diets) == 1L) "consistent" else "conflict"
    data.frame(ORDER = now$ORDER[i[1]], FAMILY = now$FAMILY[i[1]],
      GENUS = now$GENUS[i[1]], n_occurrences = length(i),
      n_diet_coded = sum(!is.na(now$DIET_2[i])),
      observed_diets = paste(diets, collapse = ";"),
      diet_status = status,
      NOW_DIET_2_GENUS = if (status == "consistent") diets else NA_character_,
      stringsAsFactors = FALSE)
  })
  audit <- do.call(rbind, audit)
  rownames(audit) <- NULL
  audit <- audit[order(audit$ORDER, audit$FAMILY, audit$GENUS), ]
  rownames(audit) <- NULL

  now$NOW_DIET_2_GENUS <- NA_character_
  now$diet_status <- ifelse(invalid, "unidentified_genus", NA_character_)
  now$diet_evidence_n <- 0L
  for (i in groups) {
    diets <- unique(stats::na.omit(now$DIET_2[i]))
    now$diet_status[i] <- if (!length(diets)) "missing" else if (length(diets) == 1L) "consistent" else "conflict"
    now$diet_evidence_n[i] <- sum(!is.na(now$DIET_2[i]))
    if (length(diets) == 1L) now$NOW_DIET_2_GENUS[i] <- diets
  }
  guild_names <- c(herbivore = "Herbivore", carnivore = "Carnivore",
    omnivore = "Omnivore", invert = "Invertivore", piscivore = "Piscivore")
  now$trophic_guild <- unname(guild_names[now$NOW_DIET_2_GENUS])
  if (anyNA(now$diet_status) || nrow(now) != sum(table(now$diet_status)))
    stop("Incomplete genus diet assignment", call. = FALSE)

  manifest <- data.frame(
    source_file = basename(source_file),
    source_md5 = unname(tools::md5sum(source_file)),
    source_description = if (basename(source_file) == "NOW_filtered_by_Shan_Oct2019.csv")
      "Local Shan October 2019 NOW extract; not a new NOW download" else
      "User-supplied NOW flat-table export; document its retrieval date and filters separately",
    source_url = "https://nowdatabase.org/now/database/",
    n_occurrences = nrow(now), n_genus_taxa = nrow(audit),
    n_consistent_genera = sum(audit$diet_status == "consistent"),
    n_conflicting_genera = sum(audit$diet_status == "conflict"),
    n_missing_genera = sum(audit$diet_status == "missing"),
    n_herbivore_occurrences = sum(now$trophic_guild == "Herbivore", na.rm = TRUE),
    n_carnivore_occurrences = sum(now$trophic_guild == "Carnivore", na.rm = TRUE),
    assignment_rule = "Unique nonmissing DIET_2 within ORDER/FAMILY/GENUS; conflict or missing means unassigned",
    stringsAsFactors = FALSE)
  utils::write.csv(now, paths[1], row.names = FALSE, na = "", fileEncoding = "UTF-8")
  utils::write.csv(audit, paths[2], row.names = FALSE, na = "", fileEncoding = "UTF-8")
  utils::write.csv(manifest, paths[3], row.names = FALSE, na = "", fileEncoding = "UTF-8")
  message("NOW master: ", nrow(now), " occurrences; ", nrow(audit),
    " genus taxa; ", manifest$n_conflicting_genera, " conflicting diet assignments")
  invisible(list(master = now, genus_audit = audit, manifest = manifest, paths = paths))
}

if (sys.nframe() == 0L) build_now_master(overwrite = TRUE)
