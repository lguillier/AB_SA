#' Fit a multinomial logistic model on full ABSA matrix (robust CSV handling)
#'
#' Reads "mnl_input_0.csv" (created by CreateInputMNL), cleans ghost/index
#' columns, coerces predictors to numeric (0/1), removes zero-variance
#' predictors, and fits nnet::multinom. Returns the fitted model with AIC.
#'
#' @param mnl_input Character. Path to "mnl_input_0.csv".
#' @param ref_level Optional character. Reference level for Source (e.g., "cattle").
#' @param maxit Integer. Max iterations for multinom. Default 1000.
#' @param seed Integer. RNG seed. Default 123.
#' @param verbose Logical. If TRUE, prints optimizer trace. Default FALSE.
#'
#' @return A fitted nnet::multinom object with an added element $AIC.
#' @importFrom nnet multinom
#' @author Laurent Guillier, \email{guillier.laurent@gmail.com}
#' @export
MNLFit <- function(mnl_input,
                   ref_level = NULL,
                   maxit = 1000,
                   seed = 123,
                   verbose = FALSE) {
  if (!requireNamespace("nnet",  quietly = TRUE)) install.packages("nnet")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  
  set.seed(seed)
  
  # 1) Lecture robuste (readr évite pas mal de bizarreries de noms)
  df <- readr::read_csv(mnl_input, show_col_types = FALSE)
  
  # 2) Purge des colonnes fantômes / index
  #    - noms vides
  #    - colonnes d'index fréquentes: ...1, X, X1, ï.., rowname, index, Strain
  empty_name_cols <- which(nchar(names(df)) == 0L)
  ghost_name_set  <- c("...1", "X", "X1", "X.1", "ï..", "ï..1", "rowname", "index", "Strain")
  ghost_by_name   <- which(names(df) %in% ghost_name_set)
  drop_idx        <- unique(c(empty_name_cols, ghost_by_name))
  if (length(drop_idx)) df <- df[, -drop_idx, drop = FALSE]
  
  # 3) Vérification de la cible
  if (!"Source" %in% names(df)) {
    stop("Column 'Source' not found in ", mnl_input,
         " — did CreateInputMNL write it?")
  }
  
  # 4) Préparateurs / prédicteurs
  x_cols <- setdiff(names(df), "Source")
  
  # Heuristique: retire toute colonne non-numérique illisible (et qui n'est pas 0/1)
  to_drop <- c()
  for (nm in x_cols) {
    v <- df[[nm]]
    v_num <- suppressWarnings(as.numeric(v))
    if (all(is.na(v_num))) {
      # garder si exactement "0"/"1"
      if (!all(v %in% c(0,1,"0","1"))) to_drop <- c(to_drop, nm)
    }
  }
  if (length(to_drop)) {
    df[, to_drop] <- NULL
    x_cols <- setdiff(x_cols, to_drop)
  }
  
  # Convertir tous les prédicteurs restants en numérique 0/1 (NA -> 0)
  for (nm in x_cols) {
    v_num <- suppressWarnings(as.numeric(df[[nm]]))
    v_num[is.na(v_num)] <- 0
    df[[nm]] <- v_num
  }
  
  # 5) Cible factor + niveau de référence si fourni
  df$Source <- as.factor(df$Source)
  if (!is.null(ref_level) && ref_level %in% levels(df$Source)) {
    df$Source <- stats::relevel(df$Source, ref = ref_level)
  }
  
  # Garde-fous
  if (nlevels(df$Source) < 2) stop("Source has <2 classes after cleaning.")
  if (length(x_cols) < 1)     stop("No predictor columns left after cleaning.")
  
  # 6) Retirer les colonnes à variance nulle
  zero_var <- vapply(df[x_cols], function(v) {
    uv <- unique(v)
    length(uv) <= 1 || all(is.na(v))
  }, logical(1))
  if (any(zero_var)) {
    df   <- df[, c(setdiff(x_cols, names(zero_var)[zero_var]), "Source"), drop = FALSE]
    x_cols <- setdiff(x_cols, names(zero_var)[zero_var])
  }
  
  # 7) Fit
  mod <- nnet::multinom(Source ~ ., data = df[, c(x_cols, "Source")],
                        maxit = maxit, trace = verbose, MaxNWts = 100000)
  
  # 8) Attache l'AIC au modèle pour accès direct via $AIC (utile pour ta boucle)
  mod$AIC <- stats::AIC(mod)
  
  message("Predictors used: ", length(x_cols),
          if (any(zero_var)) paste0(" | zero-variance removed: ", sum(zero_var)) else "")
  
  return(mod)
}
