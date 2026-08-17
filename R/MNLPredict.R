#' Predict source membership probabilities for unknown strains (robust)
#'
#' Reads the presence/absence matrix for unknown strains, aligns its columns
#' to the predictors used by the fitted multinomial model, fills any missing
#' predictors with zeros (absent gene), coerces to numeric, and returns the
#' class probabilities. Strain IDs from the first "ID-like" column are used as
#' row names when available.
#'
#' @param sporadic_pres_abs Character. Path to the CSV created by
#'   \code{CreateInputMNL} (usually "predict_sporadic.csv"). The separator may
#'   be semicolon or comma; it is auto-detected.
#' @param fitted_mnl A fitted \code{nnet::multinom} model (e.g. from \code{MNLFit}).
#'
#' @return A numeric matrix of class probabilities with rows = strains (IDs) and
#'   columns = sources (classes).
#'
#' @details
#' - Predictor names are taken from \code{terms(fitted_mnl)}.
#' - Any missing predictors in \code{sporadic_pres_abs} are added and set to 0.
#' - Extra columns are ignored.
#' - The function tries to detect an "ID-like" first column among:
#'   \code{"", "...1", "X", "ï..X", "row.names", "Row.names", "ID", "Id", "id",
#'   "Strain", "strain", "Sample", "sample", "Genome", "genome"}.
#'   If unique, it becomes the row names and is dropped from predictors.
#'
#' @author Laurent Guillier, \email{guillier.laurent@gmail.com}
#' @importFrom stats terms
#' @export
MNLPredict <- function(sporadic_pres_abs, fitted_mnl) {
  # --- detect separator
  first <- tryCatch(readLines(sporadic_pres_abs, n = 1, warn = FALSE), error = function(e) "")
  sep <- if (grepl(";", first, fixed = TRUE)) ";" else ","
  
  # --- read CSV robustly (keep colnames as-is)
  df <- if (identical(sep, ";")) {
    utils::read.csv2(sporadic_pres_abs, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.csv(sporadic_pres_abs,  check.names = FALSE, stringsAsFactors = FALSE)
  }
  
  # --- try to extract strain IDs from an "ID-like" first column
  id_like_names <- c("", "...1", "X", "ï..X", "row.names", "Row.names",
                     "ID", "Id", "id", "Strain", "strain", "Sample",
                     "sample", "Genome", "genome")
  id_col_idx <- which(tolower(names(df)) %in% tolower(id_like_names))
  if (length(id_col_idx) == 0) id_col_idx <- 1L  # fallback to first column
  cand <- df[[id_col_idx[1]]]
  if (is.character(cand) || is.factor(cand)) {
    cand_chr <- as.character(cand)
    if (!anyDuplicated(cand_chr)) {
      rownames(df) <- cand_chr
      df <- df[, -id_col_idx[1], drop = FALSE]
    }
  }
  
  # --- predictors expected by the model
  pred_needed <- tryCatch(attr(stats::terms(fitted_mnl), "term.labels"), error = function(e) NULL)
  if (is.null(pred_needed) || length(pred_needed) == 0) {
    stop("Could not extract predictor names from the fitted model via terms().")
  }
  
  # --- add missing predictors as zeros; keep only needed in the right order
  missing <- setdiff(pred_needed, names(df))
  for (m in missing) df[[m]] <- 0
  X <- df[, pred_needed, drop = FALSE]
  
  # --- coerce to numeric and replace NAs by 0 (absence)
  X[] <- lapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    v[is.na(v)] <- 0
    v
  })
  
  # --- predict probabilities
  probs <- predict(fitted_mnl, newdata = X, type = "probs")
  probs <- as.matrix(probs)
  
  # keep strain IDs as rownames if available
  if (!is.null(rownames(df))) rownames(probs) <- rownames(df)
  
  return(probs)
}
