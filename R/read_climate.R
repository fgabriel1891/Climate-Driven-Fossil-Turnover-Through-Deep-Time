# Pipeline-owned reader. Definitions only: no downloads or analysis on source.
validate_bristol_file <- function(path) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  missing <- setdiff(c("temp_mm_1_5m", "precip_mm_srf"), names(nc$var))
  if (length(missing)) stop("Missing variables: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

read_bristol_variable <- function(path, variable) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  v <- nc$var[[variable]]
  if (is.null(v)) stop("Missing variable: ", variable)
  dn <- vapply(v$dim, function(d) tolower(d$name), character(1))
  ix <- which(grepl("^longitude(_[0-9]+)?$|^lon$", dn))
  iy <- which(grepl("^latitude(_[0-9]+)?$|^lat$", dn))
  if (length(ix) != 1L || length(iy) != 1L) stop("Ambiguous spatial dimensions: ", variable)
  other <- setdiff(seq_along(v$dim), c(ix, iy))
  if (any(vapply(v$dim[other], function(d) d$len, numeric(1)) != 1)) stop("Non-spatial dimensions must be singleton")
  unit <- function(d) tolower(trimws(d$units))
  if (!unit(v$dim[[ix]]) %in% c("degrees_east", "degree_east") ||
      !unit(v$dim[[iy]]) %in% c("degrees_north", "degree_north")) stop("Coordinates must declare geographic degrees")
  x <- v$dim[[ix]]$vals
  y <- v$dim[[iy]]$vals
  if (any(!is.finite(c(x, y))) || any(abs(y) > 90) || any(x < -180 | x > 360)) stop("Invalid coordinates")
  a <- ncdf4::ncvar_get(nc, variable, collapse_degen = FALSE)
  a <- aperm(a, c(ix, iy, other))
  dim(a) <- c(length(x), length(y))
  x <- ((x + 180) %% 360) - 180
  xo <- order(x)
  yo <- order(y, decreasing = TRUE)
  x <- x[xo]
  y <- y[yo]
  step <- function(z) {
    d <- abs(diff(z))
    if (!length(d) || any(d <= 0) || max(abs(d - d[1])) > 1e-6) stop("Grid must be regular without duplicate seam cells")
    d[1]
  }
  dx <- step(x)
  dy <- step(y)
  r <- terra::rast(nrows = length(y), ncols = length(x),
    xmin = min(x) - dx/2, xmax = max(x) + dx/2,
    ymin = min(y) - dy/2, ymax = max(y) + dy/2,
    crs = "+proj=longlat +datum=WGS84 +no_defs")
  terra::values(r) <- as.vector(a[xo, yo, drop = FALSE])
  source_units <- v$units
  if (length(source_units) != 1L || is.na(source_units) || !nzchar(source_units)) stop("Missing units: ", variable)
  u <- tolower(gsub("[[:space:]_*^]", "", source_units))
  conversion_note <- "Units interpreted as declared"
  if (variable == "temp_mm_1_5m") {
    if (u %in% c("k", "kelvin", "degreekelvin", "degreeskelvin")) {
      r <- r - 273.15
    } else if (!u %in% c("c", "degc", "celsius", "degreecelsius", "degreescelsius")) stop("Unknown temperature units: ", source_units)
    output_units <- "degC"
  } else if (variable == "precip_mm_srf") {
    attribute_text <- function(name) {
      att <- ncdf4::ncatt_get(nc, variable, name)
      if (!isTRUE(att$hasatt)) return("")
      paste(att$value, collapse = " ")
    }
    flux <- identical(tolower(trimws(attribute_text("standard_name"))), "precipitation_flux")
    title <- toupper(gsub("[[:space:]]", "", paste(attribute_text("long_name"), attribute_text("title"))))
    corrected <- u %in% c("kgm-2", "kg/m2") && flux && grepl("KG/M2/S", title, fixed = TRUE)
    if (corrected) conversion_note <- "Corrected incomplete kg m-2 to kg m-2 s-1 using precipitation_flux and explicit KG/M2/S title"
    if (u %in% c("kgm-2s-1", "kg/m2/s", "kgm-2/s", "mms-1", "mm/s") || corrected) {
      r <- r * 86400
    } else if (!u %in% c("mm/day", "mmday-1", "mmd-1")) stop("Unknown precipitation units: ", source_units)
    output_units <- "mm/day"
  } else stop("Unsupported variable: ", variable)
  names(r) <- variable
  list(raster = r, source_units = source_units, output_units = output_units,
       conversion_note = conversion_note)
}

bristol_validate_decoded <- function(decoded, path) {
  if (!is.list(decoded) || !inherits(decoded$raster, "SpatRaster")) {
    stop("Incompatible climate reader result for ", path, ". Reload pipeline/load_pipeline.R.", call. = FALSE)
  }
  for (field in c("source_units", "output_units")) {
    value <- decoded[[field]]
    if (length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Climate reader returned missing/invalid ", field, " for ", path,
           ". Reload pipeline/load_pipeline.R.", call. = FALSE)
    }
  }
  # Annotation is optional for compatibility, but required units never default.
  if (is.null(decoded$conversion_note) || !length(decoded$conversion_note)) {
    decoded$conversion_note <- "Reader supplied no conversion annotation"
  } else decoded$conversion_note <- paste(decoded$conversion_note, collapse = "; ")
  decoded
}
