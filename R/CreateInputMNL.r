#' Build AB_SA input matrices from Scoary + Roary (robust, with alias→cluster mapping)
#'
#' This function:
#'  1) reads the Scoary trait file (ID + one 0/1 column per source) and derives the per-strain category;
#'  2) reads all Scoary `*results.csv` (one per source), detects p-value column, filters by `p_thresh`,
#'     selects top `maxGenes` per source by ascending p-value;
#'  3) maps Scoary gene fields (including merged names like "a-b-c" and functional aliases) to Roary clusters
#'     using `gene_presence_absence.csv` (alias map = cluster ID itself + "Non-unique Gene name" + first token of Annotation);
#'  4) unions the selected clusters across sources and extracts the presence/absence matrix from Roary `.Rtab`;
#'  5) writes:
#'        - "mnl_input_0.csv"      (training matrix + Source),
#'        - "predict_sporadic.csv" (unlabelled strains to predict),
#'        - optional "collapsed_map.csv" (mapping log: source, scoary_field, chosen cluster).
#'
#' @param traitfile   Character. Scoary trait file (CSV) with first column = strain ID, then one 0/1 column per source.
#' @param roaryRtab   Character. Path to Roary Rtab file (e.g., "gene_presence_absence.Rtab").
#' @param maxGenes    Integer. Maximum number of genes (per source) to keep (top by p-value).
#' @param p_thresh    Numeric. P-value threshold to keep Scoary rows (default 1e-2).
#' @param map_file    Character or NULL. If not NULL, write a CSV mapping log (default "collapsed_map.csv").
#' @param roary_csv   Character. Path to Roary CSV "gene_presence_absence.csv" (needed for alias mapping).
#'
#' @return (invisible) list with fields:
#'   \itemize{
#'     \item genes_final       vector of selected Roary clusters
#'     \item n_sources         number of sources (intersection: trait columns ∩ Scoary files)
#'     \item n_strains_train   number of labelled (training) strains
#'     \item n_strains_pred    number of unlabelled (to predict) strains
#'   }
#' @author Laurent Guillier, \email{guillier.laurent@gmail.com}
#' @export
CreateInputMNL <- function(traitfile,
                           roaryRtab,
                           maxGenes,
                           p_thresh = 1e-2,
                           map_file  = "collapsed_map.csv",
                           roary_csv = "gene_presence_absence.csv") {
  # -------- Utilities (local) -------------------------------------------------
  try_silent <- function(expr) tryCatch(expr, error = function(e) structure(list(e), class = "try-error"))
  norm_str   <- function(x) trimws(gsub('"', "", x, fixed = TRUE))
  detect_sep <- function(path) {
    first <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
    if (grepl(";", first)) ";" else ","
  }
  pcol_guess <- function(nms) {
    hits <- which(grepl("naive", nms, TRUE) & grepl("p", nms, TRUE))
    if (!length(hits)) hits <- which(grepl("^p[-_ ]?value$", nms, TRUE) |
                                       grepl("p[-_ ]?value", nms, TRUE)    |
                                       grepl("^p$", nms, TRUE))
    if (!length(hits)) NA_integer_ else hits[1]
  }
  get_src_key <- function(fname) tolower(sub("_.*$", "", basename(fname)))
  safe_read_rtab <- function(path) {
    x <- try_silent(readr::read_tsv(path, show_col_types = FALSE))
    if (inherits(x, "try-error")) {
      x <- try_silent(utils::read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE))
    }
    if (inherits(x, "try-error") || is.null(x)) stop("Cannot read Rtab: ", path)
    x
  }
  
  # -------- 1) Trait file, categories ----------------------------------------
  sep_trait <- detect_sep(traitfile)
  traits <- utils::read.csv(traitfile, header = TRUE, sep = sep_trait,
                            check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(traits) < 2) stop("Trait file must have ID + ≥1 source columns.")
  id_col <- colnames(traits)[1]
  trait_names <- tolower(colnames(traits)[-1])
  
  if (!exists("DefStrainCategory")) stop("DefStrainCategory() must be sourced before.")
  categories <- DefStrainCategory(traitfile)
  if (length(categories) != nrow(traits)) {
    stop("DefStrainCategory() does not return the same number of elements as trait rows.")
  }
  id_categories <- data.frame(id = traits[[id_col]], cat = categories, stringsAsFactors = FALSE)
  
  # -------- 2) Alias→cluster map from Roary CSV ------------------------------
  if (!file.exists(roary_csv)) stop("Roary CSV not found: ", roary_csv)
  suppressMessages({
    ro_csv <- readr::read_csv(roary_csv, show_col_types = FALSE)
  })
  needed_cols <- c("Gene","Non-unique Gene name","Annotation")
  if (!all(needed_cols %in% names(ro_csv))) {
    stop("Missing columns in Roary CSV: ", paste(setdiff(needed_cols, names(ro_csv)), collapse = ", "))
  }
  ro_csv$Gene <- norm_str(ro_csv$Gene)
  
  # identity alias = cluster
  id_map <- data.frame(alias = ro_csv$Gene, cluster = ro_csv$Gene, stringsAsFactors = FALSE)
  
  # explicit aliases from non-unique names (split ';')
  alias_exp <- tidyr::separate_rows(
    ro_csv[, c("Gene","Non-unique Gene name")],
    `Non-unique Gene name`, sep = ";", convert = FALSE
  )
  alias_exp$`Non-unique Gene name` <- norm_str(alias_exp$`Non-unique Gene name`)
  alias_exp <- subset(alias_exp, nzchar(`Non-unique Gene name`))
  alias_exp <- unique(data.frame(
    alias   = alias_exp$`Non-unique Gene name`,
    cluster = alias_exp$Gene,
    stringsAsFactors = FALSE
  ))
  
  # conservative alias from first token of Annotation
  annot_tok <- trimws(gsub("[^A-Za-z0-9_]+", "",
                           vapply(strsplit(ifelse(is.na(ro_csv$Annotation), "", ro_csv$Annotation), " "),
                                  `[`, "", 1)))
  ann_map <- unique(data.frame(alias = annot_tok, cluster = ro_csv$Gene, stringsAsFactors = FALSE))
  ann_map <- subset(ann_map, nzchar(alias))
  
  alias_map <- rbind(id_map, alias_exp, ann_map)
  alias_map <- alias_map[!duplicated(alias_map$alias), ]
  
  map_scoary_gene <- function(gfield) {
    toks <- unlist(strsplit(gfield, "[,;|\\-]+"))
    toks <- norm_str(toks[nzchar(toks)])
    hit  <- alias_map$cluster[match(toks, alias_map$alias)]
    hit  <- hit[!is.na(hit)]
    if (length(hit)) hit[1] else NA_character_
  }
  
  # -------- 3) Read Scoary files & select clusters ---------------------------
  files <- list.files(pattern = "results\\.csv$", ignore.case = TRUE)
  if (!length(files)) stop("No Scoary results.csv found in working dir.")
  
  # intersection by name
  src_from_files <- unique(vapply(files, get_src_key, character(1)))
  keep_sources <- intersect(trait_names, src_from_files)
  if (length(keep_sources) < 2) {
    stop("No overlap between trait columns and Scoary files (by name).")
  }
  files_keep <- files[vapply(files, function(f) get_src_key(f) %in% keep_sources, logical(1))]
  ord <- match(vapply(files_keep, get_src_key, character(1)), keep_sources)
  files_keep <- files_keep[order(ord)]
  
  all_enriched <- setNames(vector("list", length(keep_sources)), keep_sources)
  map_log <- list()
  
  for (f in files_keep) {
    src <- get_src_key(f)
    sc <- try_silent(readr::read_delim(f, delim = ",", quote = "\"", show_col_types = FALSE))
    if (inherits(sc, "try-error") || is.null(sc) || nrow(sc) == 0) {
      sc <- try_silent(readr::read_delim(f, delim = ";", quote = "\"", show_col_types = FALSE))
    }
    if (inherits(sc, "try-error") || is.null(sc) || nrow(sc) == 0) {
      warning("Failed to read Scoary file: ", f); next
    }
    nms <- names(sc)
    pcol_i <- pcol_guess(nms)
    if (is.na(pcol_i)) { warning("No p-value column in ", f); next }
    p <- suppressWarnings(as.numeric(sc[[pcol_i]]))
    gcol <- if ("Gene" %in% nms) "Gene" else nms[1]
    
    # filter by p-value threshold
    keep <- which(!is.na(p) & p <= p_thresh)
    if (!length(keep)) { warning("No genes <= p_thresh in ", f); next }
    sc_f <- sc[keep, , drop = FALSE]
    p_f  <- p[keep]
    
    # order by ascending p
    ordp <- order(p_f)
    sc_f <- sc_f[ordp, , drop = FALSE]
    
    # select top 'maxGenes' and map to clusters
    sel_raw <- head(sc_f[[gcol]], maxGenes)
    sel_map <- unique(vapply(sel_raw, map_scoary_gene, character(1)))
    sel_map <- sel_map[!is.na(sel_map) & nzchar(sel_map)]
    
    all_enriched[[src]] <- sel_map
    
    if (length(sel_raw)) {
      # keep log with same length; pad chosen with NA if needed
      chosen <- rep(NA_character_, length(sel_raw))
      chosen[seq_len(min(length(sel_map), length(sel_raw)))] <- sel_map
      map_log[[src]] <- data.frame(source = src,
                                   scoary_field = sel_raw,
                                   chosen = chosen,
                                   stringsAsFactors = FALSE)
    }
  }
  
  genes_final <- unique(unlist(all_enriched, use.names = FALSE))
  if (!length(genes_final)) {
    stop("No enriched genes retained (after p_thresh & mapping).")
  }
  
  # -------- 4) Read Roary Rtab & build matrices ------------------------------
  roaryRtab_df <- safe_read_rtab(roaryRtab)
  if (!"Gene" %in% names(roaryRtab_df)) stop("Column 'Gene' not found in Rtab.")
  roaryRtab_df$Gene <- norm_str(roaryRtab_df$Gene)
  
  idx <- match(genes_final, roaryRtab_df$Gene)
  missing <- is.na(idx)
  if (any(missing)) {
    cat("Example missing clusters in Rtab:\n",
        paste(head(genes_final[missing], 10), collapse = ", "), "\n")
    stop("Selected genes not found in Roary Rtab.")
  }
  mnl_input <- roaryRtab_df[idx, , drop = FALSE]
  
  # transpose: strains as rows, gene clusters as columns
  # mnl_input: rows = genes, cols = Gene + strainIDs
  mat <- as.data.frame(t(mnl_input[, -1, drop = FALSE]), stringsAsFactors = FALSE)
  colnames(mat) <- mnl_input$Gene
  rownames(mat) <- colnames(roaryRtab_df)[-1]  # strain IDs from Rtab header
  
  # match labelled strains
  row_ids <- rownames(mat)
  pos <- match(id_categories$id, row_ids)
  pos <- pos[!is.na(pos)]
  
  t_mnl_input <- mat[pos, , drop = FALSE]
  t_mnl_input$Source <- id_categories$cat
  
  utils::write.csv(t_mnl_input, file = "mnl_input_0.csv", row.names = TRUE)
  
  # prediction = all remaining strains
  t_mnl_predict <- mat[setdiff(seq_len(nrow(mat)), pos), , drop = FALSE]
  utils::write.csv2(t_mnl_predict, file = "predict_sporadic.csv", row.names = TRUE)
  
  # mapping log (optional)
  if (!is.null(map_file)) {
    if (length(map_log)) {
      log_df <- do.call(rbind, map_log)
      utils::write.csv(log_df, map_file, row.names = FALSE)
    } else {
      utils::write.csv(data.frame(source=character(), scoary_field=character(), chosen=character()),
                       map_file, row.names = FALSE)
    }
  }
  
  message(sprintf("Selected clusters (unique): %d | present in Rtab: %d | strains matched: %d | predict strains: %d",
                  length(genes_final), nrow(mnl_input), nrow(t_mnl_input), nrow(t_mnl_predict)))
  
  invisible(list(
    genes_final       = genes_final,
    n_sources         = length(keep_sources),
    n_strains_train   = nrow(t_mnl_input),
    n_strains_pred    = nrow(t_mnl_predict)
  ))
}
