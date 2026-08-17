# =========================
# helpers.R — unified utils
# =========================
# Internal (already present in your file)
.detect_sep <- function(path) {
  first <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  if (grepl(";", first)) ";" else ","
}

.normalize_gene <- function(x) {
  x <- gsub('"', "", x, fixed = TRUE)
  x <- trimws(x)
  x
}

.split_collapsed <- function(s) {
  # Split collapsed gene lists like "geneA-geneB|geneC,geneD"
  unlist(strsplit(s, "[-,;|]+"))
}

.find_pval_col <- function(nms) {
  # Prioritise 'naive p-value' if present; otherwise any p-value column
  hits <- which(grepl("naive", nms, ignore.case = TRUE) & grepl("p", nms, ignore.case = TRUE))
  if (length(hits) == 0) hits <- which(grepl("^p[-_ ]?value$", nms, ignore.case = TRUE))
  if (length(hits) == 0) hits <- which(grepl("p[-_ ]?value", nms, ignore.case = TRUE))
  if (length(hits) == 0) hits <- which(grepl("^p$", nms, ignore.case = TRUE))
  if (length(hits) == 0) NA_integer_ else hits[1]
}

# ---- Public aliases (no leading dot) for easier reuse elsewhere ----
# Keep your codebase stable even if some scripts call the non-dotted names.

#' Normalise gene names: trim + remove quotes
#' @param x character vector
#' @return character vector
if (!exists("norm_gene", mode = "function")) {
  norm_gene <- function(x) .normalize_gene(x)
}

#' Split a collapsed gene-field into tokens (commas/semicolons/pipes)
#' @param s character vector (each element may contain multiple genes)
#' @return character vector of tokens
if (!exists("split_collapse", mode = "function")) {
  split_collapse <- function(s) .split_collapsed(s)
}

#' Find the p-value column index in a Scoary results data frame
#' @param nms character vector of column names
#' @return integer index (or NA_integer_ if not found)
if (!exists("find_pval_col", mode = "function")) {
  find_pval_col <- function(nms) .find_pval_col(nms)
}

## ========= Utilitaire: prédire sur tout l'entraînement avec IDs =========
PredictOnTrainingAll <- function(model,
                                 mnl_input   = "mnl_input_0.csv",
                                 traitfile,
                                 roary_rtab) {
  # helpers
  detect_sep <- function(p){ if (grepl(";", readLines(p,1,warn=FALSE))) ";" else "," }
  clean_id   <- function(x) sub("\\.(fa|fna|fasta)(\\.gz)?$","",basename(x), ignore.case=TRUE)
  
  # 1) Lire l'input d'entraînement (mnl_input_0.csv)
  mi <- read.csv(mnl_input, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Source" %in% names(mi)) stop("Column 'Source' not found in ", mnl_input)
  X <- mi[, setdiff(names(mi), "Source"), drop = FALSE]
  # coercition sûre en numérique
  X[] <- lapply(X, function(v){ if(!is.numeric(v)) suppressWarnings(as.numeric(v)) else v })
  
  # Aligner les colonnes sur celles du modèle (au cas où)
  needed <- attr(stats::terms(model), "term.labels")
  if (length(needed) == 0) stop("No predictors found in model terms().")
  missing <- setdiff(needed, names(X))
  for (m in missing) X[[m]] <- 0
  X <- X[, needed, drop = FALSE]
  
  # 2) Prédire les probabilités in-sample
  probs <- predict(model, newdata = X, type = "probs")
  probs <- as.matrix(probs)
  
  # 3) Reconstruire les IDs exacts des souches d'entraînement (ordre Roary filtré par trait)
  traits <- read.table(traitfile, header=TRUE, sep=detect_sep(traitfile),
                       check.names=FALSE, stringsAsFactors=FALSE)
  ids_trait <- as.character(traits[[1]])
  
  rtab <- read.table(roary_rtab, header=TRUE, sep="\t",
                     check.names=FALSE, stringsAsFactors=FALSE)
  roary_ids <- setdiff(colnames(rtab), "Gene")
  
  train_ids <- roary_ids[ clean_id(roary_ids) %in% clean_id(ids_trait) ]
  if (length(train_ids) != nrow(probs)) {
    stop("Nombre d'IDs reconstruits (", length(train_ids),
         ") différent du nombre de lignes dans ", mnl_input, " (", nrow(probs),
         "). Vérifie traitfile / Rtab.")
  }
  
  rownames(probs) <- train_ids
  
  # 4) Classe prédite + proba max + assemblage
  pred_class <- colnames(probs)[ max.col(probs, ties.method = "first") ]
  max_prob   <- apply(probs, 1, max)
  
  out <- data.frame(
    Strain     = train_ids,
    TrueSource = mi$Source,
    PredClass  = pred_class,
    MaxProb    = max_prob,
    probs,
    check.names = FALSE
  )
  
  # 5) Sauvegarde
  write.csv2(out, "train_predictions_fullfit.csv", row.names = FALSE)
  message("Écrit: train_predictions_fullfit.csv (", nrow(out), " souches).")
  invisible(out)
}

