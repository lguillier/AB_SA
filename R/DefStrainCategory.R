#' For each individual, this function defines the category (as text) for the
#' multinomial logistic regression from a Scoary trait file (0/1 or categorical).
#'
#' @param filenametrait Character. Path to the trait file produced by Scoary.
#'   The function auto-detects the separator (semicolon or comma).
#'
#' @return A character vector `cat_i` of the same length as the number of rows
#'   in the trait file, giving the category (source) for each strain. If a row
#'   has multiple sources marked as 1, the first matching source column (from
#'   left to right) is chosen. If no source is marked as 1 on a row, `NA` is
#'   returned for that row.
#'
#' @details
#' - The function supports two common Scoary trait file layouts:
#'   1) **Categorical format**: two columns (ID, Source-as-text). In this case,
#'      the function returns the second column as-is (coerced to `character`).
#'   2) **Binary format**: first column is the strain ID, remaining columns are
#'      sources with values in {0,1}. For each row, the function returns the
#'      *first* source whose value is 1.
#' - The first column (IDs) may have no header in some CSVs; in that case, it is
#'   temporarily named `"X"` internally. The returned vector does not include IDs.
#' - The separator is auto-detected from the first line: `;` if present, otherwise `,`.
#' - All numeric-like columns are coerced to numeric before testing equality to 1.
#'
#' @author Laurent Guillier, \email{guillier.laurent@gmail.com}
#'
#' @export
#'
DefStrainCategory <- function(filenametrait) {
  sep <- .detect_sep(filenametrait)
  df  <- readr::read_delim(filenametrait, delim = sep, show_col_types = FALSE, quote = "\"", escape_double = TRUE)
  if (ncol(df) < 1L) stop("Trait file has no columns.")
  if (!nzchar(names(df)[1])) names(df)[1] <- "X"
  
  # Case A: (ID, Source-as-text)
  if (ncol(df) == 2L && !is.numeric(df[[2]])) {
    return(as.character(df[[2]]))
  }
  # Case B: binary 0/1
  if (ncol(df) < 2L) stop("Binary trait format expected but only one column found.")
  mat <- as.data.frame(df[, -1, drop = FALSE])
  for (j in seq_along(mat)) mat[[j]] <- suppressWarnings(as.numeric(mat[[j]]))
  mm <- as.matrix(mat)
  pick <- apply(mm, 1, function(v) { w <- which(v == 1); if (length(w)) w[1] else NA_integer_ })
  cats <- colnames(mat)[pick]
  if (length(cats) != nrow(df)) stop("Categories length mismatch.")
  cats
}