#' Train/test bootstrap for multinomial logistic regression (AB_SA)
#'
#' Repeats \code{nboot} random, stratified holdout splits (train/test) of the
#' data built from \code{mnl_input}, fits a multinomial model with
#' \code{nnet::multinom} on the training set, and evaluates test accuracy and
#' per-class balanced accuracy on the test set. Returns accuracy quantiles,
#' median balanced accuracies per class, and a density object for accuracies.
#'
#' @param mnl_input Character. Path to the CSV created by \code{CreateInputMNL}
#'   (default: \code{"mnl_input_0.csv"}). Must contain a \code{Source} column
#'   (response) and numeric predictors (gene presence/absence).
#' @param percent_cross Numeric. Training proportion. Accepts \code{0–1}
#'   (e.g. \code{0.7}) or \code{1–100} (e.g. \code{70}). Default: \code{0.70}.
#' @param nboot Integer. Number of bootstrap repetitions. Default: \code{100}.
#' @param seed Integer. RNG seed for reproducibility. Default: \code{123}.
#' @param maxit Integer. Maximum iterations for \code{nnet::multinom}.
#'   Default: \code{1000}.
#' @param verbose Logical. If \code{TRUE}, prints progress every ~10% of runs.
#'   Default: \code{TRUE}.
#'
#' @details
#' \itemize{
#'   \item Stratified sampling preserves class proportions in the training set.
#'   \item Only numeric predictors are used; non-numeric columns are coerced when
#'         possible and zero-variance predictors (computed on the training set)
#'         are dropped before fitting.
#'   \item Iterations where the model fails to fit or predict are skipped; the
#'         returned summaries aggregate successful runs only.
#'   \item Balanced accuracy per class is computed one-vs-all as
#'         \eqn{(sensitivity + specificity)/2}.
#' }
#'
#' @return A \code{list} with both named and positional entries:
#' \itemize{
#'   \item \code{accuracy_quantiles} (\code{[[1]]}): named vector of test accuracy
#'         quantiles at \code{2.5\%}, \code{5\%}, \code{50\%}, \code{95\%}, \code{97.5\%}.
#'   \item \code{class_balanced_median} (\code{[[2]]}): named vector of median
#'         per-class balanced accuracies (one-vs-all).
#'   \item \code{accuracy_density} (\code{[[3]]}): \code{stats::density} object
#'         for the distribution of test accuracies.
#'   \item \code{accuracy_all}: numeric vector of accuracy values across runs.
#'   \item \code{class_balanced_all}: matrix (\code{nboot} x \code{K}) of
#'         per-class balanced accuracies across runs (rows = runs, cols = classes).
#' }
#'
#' @examples
#' \dontrun{
#' testedMNL <- MNLTrainTest("mnl_input_0.csv", percent_cross = 0.70, nboot = 100)
#' testedMNL[[1]]   # accuracy quantiles
#' testedMNL[[2]]   # median balanced accuracy per class
#' plot(testedMNL[[3]], main = "Test accuracy density")
#' }
#'
#' @seealso \code{\link{MNLFit}}, \code{\link{MNLPredict}}, \code{\link{CreateInputMNL}}
#'
#' @importFrom nnet multinom
#' @importFrom stats quantile density
#' @importFrom utils read.csv
#' @author Laurent Guillier
#' @export
#'
MNLTrainTest <- function(mnl_input = "mnl_input_0.csv",
                         percent_cross = 0.70,
                         nboot = 100,
                         seed = 123,
                         maxit = 1000,
                         verbose = TRUE) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' is required.")
  # --- read CSV robustly
  df <- tryCatch({
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::read_csv(mnl_input, show_col_types = FALSE)
    } else {
      utils::read.csv(mnl_input, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(e) utils::read.csv(mnl_input, check.names = FALSE, stringsAsFactors = FALSE))
  
  # --- fix column names: trim, replace empty, ensure unique
  nm <- trimws(names(df))
  nm[nm == ""] <- paste0("V", seq_len(sum(nm == "")))
  names(df) <- make.unique(nm)
  
  # --- find 'Source' column robustly
  src_idx <- which(tolower(trimws(names(df))) == "source")
  if (length(src_idx) < 1) stop("Column 'Source' not found in ", mnl_input)
  if (length(src_idx) > 1) src_idx <- src_idx[1]
  
  # response
  y <- factor(df[[src_idx]])
  levels_y <- levels(y)
  
  # predictors = everything except Source
  X <- df[, setdiff(seq_along(df), src_idx), drop = FALSE]
  
  # keep only numeric predictors; try coercion for 0/1-like text
  num <- vapply(X, is.numeric, logical(1))
  if (any(!num)) {
    X[!num] <- lapply(X[!num], function(v) suppressWarnings(as.numeric(v)))
    num <- vapply(X, is.numeric, logical(1))
  }
  X <- X[, num, drop = FALSE]
  if (ncol(X) == 0) stop("No numeric predictors found after cleaning.")
  
  # train proportion (accept 0–1 or 1–100)
  if (percent_cross > 1) percent_cross <- percent_cross / 100
  percent_cross <- max(min(percent_cross, 0.99), 0.01)
  
  set.seed(seed)
  
  # helpers
  zero_var_cols <- function(M) {
    vapply(as.data.frame(M), function(x) {
      ux <- unique(x)
      length(ux) <= 1 || all(is.na(x))
    }, logical(1))
  }
  stratified_indices <- function(y, p) {
    idx_train <- integer(0)
    for (lv in levels(y)) {
      idc <- which(y == lv); n <- length(idc)
      if (n == 0) next
      k <- max(1L, floor(p * n))
      if (k >= n) k <- n - 1L
      if (k > 0) idx_train <- c(idx_train, sample(idc, k))
    }
    unique(idx_train)
  }
  ba_per_class <- function(true, pred, lvls) {
    cm <- table(factor(true, lvls), factor(pred, lvls))
    N  <- sum(cm)
    out <- numeric(length(lvls)); names(out) <- lvls
    for (i in seq_along(lvls)) {
      TP <- cm[i, i]; FN <- sum(cm[i, ]) - TP; FP <- sum(cm[, i]) - TP
      TN <- N - TP - FN - FP
      sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
      out[i] <- mean(c(sens, spec), na.rm = TRUE)
    }
    out
  }
  
  acc_all <- numeric(0)
  ba_all  <- matrix(NA_real_, nrow = nboot, ncol = length(levels_y),
                    dimnames = list(NULL, levels_y))
  eff <- 0L
  
  for (b in seq_len(nboot)) {
    tr_idx <- stratified_indices(y, percent_cross)
    te_idx <- setdiff(seq_len(nrow(X)), tr_idx)
    if (length(tr_idx) < 2L || length(te_idx) < 1L) next
    
    Xtr <- X[tr_idx, , drop = FALSE]; ytr <- y[tr_idx]
    Xte <- X[te_idx, , drop = FALSE]; yte <- y[te_idx]
    
    zv <- zero_var_cols(Xtr)
    if (any(zv)) { Xtr <- Xtr[, !zv, drop = FALSE]; Xte <- Xte[, !zv, drop = FALSE] }
    if (ncol(Xtr) == 0) next
    
    fit <- try(nnet::multinom(ytr ~ ., data = data.frame(Xtr, ytr),
                              maxit = maxit, trace = FALSE, MaxNWts = 100000),
               silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    yhat <- try(predict(fit, newdata = data.frame(Xte), type = "class"), silent = TRUE)
    if (inherits(yhat, "try-error")) next
    
    acc <- mean(as.character(yhat) == as.character(yte))
    ba  <- ba_per_class(yte, yhat, levels_y)
    
    eff <- eff + 1L
    acc_all[eff]  <- acc
    ba_all[eff, ] <- ba
    
    if (verbose && (b %% max(1, floor(nboot/10)) == 0)) {
      message("Iter ", b, "/", nboot, " — acc=", sprintf("%.3f", acc))
    }
  }
  
  if (eff == 0L) stop("No successful fit across iterations. Check class balance or increase maxit.")
  
  acc_all <- acc_all[seq_len(eff)]
  ba_all  <- ba_all[seq_len(eff), , drop = FALSE]
  
  acc_q  <- stats::quantile(acc_all, probs = c(.025, .05, .50, .95, .975), na.rm = TRUE)
  ba_med <- apply(ba_all, 2, function(v) stats::median(v, na.rm = TRUE))
  dens   <- stats::density(acc_all)
  
  out <- list(
    accuracy_quantiles     = acc_q,
    class_balanced_median  = ba_med,
    accuracy_density       = dens,
    accuracy_all           = acc_all,
    class_balanced_all     = ba_all
  )
  out[[1]] <- out$accuracy_quantiles
  out[[2]] <- out$class_balanced_median
  out[[3]] <- out$accuracy_density
  out
}
