#' Shared primitives for the annotation-entropy metrics
#'
#' Port of the Python package's `entropy_core.py`. Python is the reference for the
#' entropy metrics, as R is for the annotation ones, so these functions follow it
#' formula for formula and the two should agree to floating-point noise.
#'
#' Unresolvable cell types: refuse the metric, do not guess.
#' A Rao metric needs a distance for every pair of candidate cell types, and a
#' cell type with no Cell Ontology term or no row in the marker table has none.
#' The policy is the one already applied to disease conditions, where no cell
#' ontology exists at all: if any candidate cell type in the analysis cannot be
#' placed, the affected metric is not computed for that analysis.
#'
#' @keywords internal
#' @name entropy_core
NULL

#' @rdname entropy_core
ACCORDION_METRICS <- c("Shannon", "Gini", "Ontology-sp-Rao", "Ontology-Lin-Rao",
                       "Jaccard-Rao", "Overlap-Rao")

#' @rdname entropy_core
ACCORDION_ONTOLOGY_METRICS <- c("Ontology-sp-Rao", "Ontology-Lin-Rao")

#' @rdname entropy_core
ACCORDION_MARKER_METRICS <- c("Jaccard-Rao", "Overlap-Rao")

.missing_value <- function(missing) {
  switch(missing, nan = NA_real_, zero = 0, max = 1,
         stop("missing must be one of 'nan', 'zero', 'max'"))
}

# --- probability vector ---

#' @rdname entropy_core
clean_probabilities <- function(probabilities) {
  p <- as.numeric(probabilities)
  p <- p[is.finite(p) & p > 0]
  n <- length(p)
  if (n <= 1) return(list(p = NULL, n = n))
  total <- sum(p)
  if (!is.finite(total) || total <= 0) return(list(p = NULL, n = n))
  list(p = p / total, n = n)
}

#' @rdname entropy_core
normalized_shannon <- function(probabilities) {
  cp <- clean_probabilities(probabilities)
  if (is.null(cp$p)) return(0)
  # Shannon is reported normalized: raw it ranges over [0, log2(n)] and so grows
  # with the number of candidate cell types, which has nothing to do with how
  # ambiguous the annotation is.
  max(0, sum(-cp$p * log2(cp$p)) / log2(cp$n))
}

#' @rdname entropy_core
raw_gini <- function(probabilities) {
  cp <- clean_probabilities(probabilities)
  if (is.null(cp$p)) return(0)
  # Gini is reported raw, the opposite convention to Shannon, because its range
  # [0, 1 - 1/n] is already bounded by 1 for any n. Dividing by 1 - 1/n would
  # inflate small clusters: with 2 candidates the maximum is 0.5.
  max(0, 1 - sum(cp$p^2))
}

#' @rdname entropy_core
rao <- function(p, D) {
  p <- as.numeric(p)
  if (length(p) == 0 || length(D) == 0 || nrow(D) != length(p)) return(0)
  as.numeric(t(p) %*% D %*% p)
}

# --- distance matrices ---

#' @rdname entropy_core
marker_sets <- function(dt_all_marker, annot_name, require_positive_impact = TRUE) {
  if (is.null(dt_all_marker) || nrow(dt_all_marker) == 0) return(list())
  dt <- data.table::as.data.table(dt_all_marker)
  if (!(annot_name %in% names(dt)) && "annotation_per_cell" %in% names(dt)) {
    annot_name <- "annotation_per_cell"
  }
  if (!(annot_name %in% names(dt))) return(list())
  if (require_positive_impact && "gene_impact_score_per_cluster" %in% names(dt)) {
    dt <- dt[get("gene_impact_score_per_cluster") > 0]
  }
  if (nrow(dt) == 0) return(list())
  # Markers carry their sign, so a gene marking one cell type positively and
  # another negatively counts as a difference, not as shared evidence.
  sign <- ifelse(dt$marker_type == "positive", "+",
                 ifelse(dt$marker_type == "negative", "-", ""))
  signed <- paste0(ifelse(is.na(dt$marker), "", as.character(dt$marker)), sign)
  split(signed, as.character(dt[[annot_name]])) |> lapply(unique)
}

#' @rdname entropy_core
marker_distance_matrix <- function(celltypes, sets, kind = "Overlap-Rao",
                                   missing = "max") {
  celltypes <- as.character(celltypes)
  n <- length(celltypes)
  D <- matrix(0, n, n, dimnames = list(celltypes, celltypes))
  fill <- .missing_value(missing)
  if (n < 2) return(D)
  for (i in seq_len(n - 1)) {
    a <- sets[[celltypes[i]]]
    for (j in seq.int(i + 1L, n)) {
      b <- sets[[celltypes[j]]]
      if (is.null(a) || is.null(b)) {
        d <- fill
      } else {
        inter <- length(intersect(a, b))
        denom <- if (kind == "Jaccard-Rao") length(union(a, b)) else min(length(a), length(b))
        d <- if (denom > 0) 1 - inter / denom else 1
      }
      D[i, j] <- D[j, i] <- d
    }
  }
  D
}

# Cached derived structures of an ontology_index: the is_a graph and, for Lin,
# the number of subclasses of each term. Both are expensive and depend only on
# the ontology object, so they are built once per session.
.onto_cache <- new.env(parent = emptyenv())

.onto_derived <- function(ontology) {
  key <- paste0("onto_", length(ontology$id), "_", length(unlist(ontology$parents)))
  if (!is.null(.onto_cache[[key]])) return(.onto_cache[[key]])

  # The term universe: real, non-deprecated Cell Ontology terms. ontology_index
  # keeps two kinds of entry that are not terms of the ontology proper, and both
  # would land in the information-content denominator:
  #   * 179 obsolete CL terms, which are isolated (no is_a edge touches them)
  #   * 13 relation and property names swept into $id (part_of, RO:0002202,
  #     has_ontology_root_term, nine CP:...)
  # Dropping them leaves exactly the 2777 terms obonet reads from cl-basic.obo,
  # verified as a set. It matters because IC = -log((subclasses + 1) / total): a
  # wrong total shifts every IC by a constant, and Lin = 2*IC(mica)/(IC(a)+IC(b))
  # does not cancel a constant shift. With the full 2969 the Lin distances were
  # off by ~0.007 against the Python reference.
  obsolete <- unlist(ontology$obsolete)
  if (is.null(obsolete)) obsolete <- rep(FALSE, length(ontology$id))
  ids <- ontology$id[!obsolete & grepl("^CL:", ontology$id)]
  # Edges run child -> parent, as in the OBO file and in obonet. ontology_index's
  # $parents holds is_a only, which is what both metrics want: develops_from is
  # not a subsumption relation and has no business in a subclass count.
  ch <- rep(names(ontology$parents), lengths(ontology$parents))
  pa <- unlist(ontology$parents, use.names = FALSE)
  keep <- !is.na(pa) & pa %in% ids & ch %in% ids
  g <- igraph::graph_from_data_frame(
    data.frame(from = ch[keep], to = pa[keep], stringsAsFactors = FALSE),
    directed = TRUE, vertices = data.frame(name = ids, stringsAsFactors = FALSE))

  # $ancestors[[t]] is t and all its superclasses. A term's subclass count is how
  # many terms list it among their ancestors, minus itself -- counted over the
  # term universe above, as Python counts nx.ancestors over the graph's nodes.
  anc <- lapply(ontology$ancestors[ids], function(a) a[a %in% ids])
  n_sub <- table(unlist(anc, use.names = FALSE))
  subclasses <- setNames(rep(0L, length(ids)), ids)
  common <- intersect(names(n_sub), ids)
  subclasses[common] <- as.integer(n_sub[common]) - 1L

  total <- length(ids)
  IC <- setNames(-log((subclasses + 1) / total), ids)

  out <- list(graph = g, undirected = igraph::as.undirected(g, mode = "collapse"),
              IC = IC, ancestors = anc, total = total)
  .onto_cache[[key]] <- out
  out
}

#' @rdname entropy_core
ontology_distance_matrix <- function(celltypes, ontology, kind = "Ontology-sp-Rao",
                                     missing = "max") {
  celltypes <- as.character(celltypes)
  n <- length(celltypes)
  D <- matrix(0, n, n, dimnames = list(celltypes, celltypes))
  if (is.null(ontology) || n == 0) return(D)
  fill <- .missing_value(missing)

  lab2id <- setNames(names(ontology$name), unlist(ontology$name, use.names = FALSE))
  ids <- lab2id[celltypes]
  der <- .onto_derived(ontology)

  if (kind == "Ontology-sp-Rao") {
    # Edge counts on the undirected graph, divided by the matrix maximum, so the
    # values are relative to the cell types present in THIS run.
    ok <- !is.na(ids)
    if (any(ok)) {
      sp <- igraph::distances(der$undirected, v = ids[ok], to = ids[ok])
      sub <- matrix(fill, n, n)
      sub[ok, ok] <- sp
      sub[!is.finite(sub) & !is.na(sub)] <- fill
      D[] <- sub
    } else {
      D[] <- fill
    }
    diag(D) <- 0
    mx <- suppressWarnings(max(D[is.finite(D)], na.rm = TRUE))
    if (is.finite(mx) && mx > 0) D <- D / mx
    return(D)
  }

  # --- Lin ---
  IC <- der$IC; anc <- der$ancestors
  for (i in seq_len(max(0L, n - 1L))) {
    for (j in seq.int(i + 1L, n)) {
      if (is.na(ids[i]) || is.na(ids[j])) {
        d <- fill
      } else {
        common <- intersect(anc[[ids[i]]], anc[[ids[j]]])
        if (!length(common)) {
          d <- 1
        } else {
          mica <- max(IC[common], na.rm = TRUE)
          denom <- IC[[ids[i]]] + IC[[ids[j]]]
          d <- 1 - if (denom == 0) 0 else (2 * mica) / denom
        }
      }
      D[i, j] <- D[j, i] <- d
    }
  }
  diag(D) <- 0
  D
}

#' @rdname entropy_core
distance_matrix <- function(celltypes, kind, ontology = NULL, sets = NULL,
                            missing = "max") {
  if (kind %in% ACCORDION_ONTOLOGY_METRICS) {
    return(ontology_distance_matrix(celltypes, ontology, kind, missing = missing))
  }
  if (kind %in% ACCORDION_MARKER_METRICS) {
    return(marker_distance_matrix(celltypes, if (is.null(sets)) list() else sets,
                                  kind, missing = missing))
  }
  matrix(0, length(celltypes), length(celltypes))
}

# --- preconditions: can these metrics be computed at all? ---

#' @rdname entropy_core
unmapped_celltypes <- function(celltypes, ontology) {
  celltypes <- unique(as.character(celltypes))
  if (is.null(ontology)) return(sort(celltypes))
  labels <- unlist(ontology$name, use.names = FALSE)
  sort(celltypes[!(celltypes %in% labels)])
}

#' @rdname entropy_core
celltypes_without_markers <- function(celltypes, sets) {
  celltypes <- unique(as.character(celltypes))
  sort(celltypes[vapply(celltypes, function(c) length(sets[[c]]) == 0, logical(1))])
}

#' @rdname entropy_core
skip_message <- function(metrics, offenders, reason, limit = 10) {
  shown <- paste(utils::head(offenders, limit), collapse = ", ")
  if (length(offenders) > limit) {
    shown <- paste0(shown, ", ... (+", length(offenders) - limit, " more)")
  }
  paste0("WARNING: not computing ", paste(metrics, collapse = ", "),
         " for this analysis -- ", length(offenders), " cell type(s) ", reason,
         ": ", shown, ". The remaining entropy metrics are unaffected.")
}

# --- driver: entropie per cluster ---

#' Cluster-level entropies
#'
#' Port of the cluster-level block of the Python package. `info` must hold one row
#' per (cluster, candidate cell type) with columns `cluster`, `celltype`,
#' `percentage` and `celltype_impact_score`.
#'
#' `missing = "nan"` reproduces the Python path: a pair that cannot be resolved is
#' NA and propagates into Rao's Q rather than being silently treated as distance 0.
#' It can only fire for two cell types that both exist but have no connecting path,
#' because the preconditions below refuse the metric otherwise.
#'
#' @return named list of data.tables, `cluster_entropy_<metric>`
#' @keywords internal
compute_cluster_entropies <- function(info, entropy, ontology = NULL,
                                      dt_all_marker = NULL,
                                      min_percentage_celltype_entropy = 0.05,
                                      verbose = TRUE) {
  info <- data.table::as.data.table(info)
  if (identical(entropy, "all")) entropy <- ACCORDION_METRICS
  entropy <- intersect(entropy, ACCORDION_METRICS)
  if (!length(entropy) || !nrow(info)) return(list())

  info <- info[percentage >= min_percentage_celltype_entropy &
                 !grepl("unknown", celltype, ignore.case = TRUE)]
  if (!nrow(info)) return(list())

  out <- list()
  clusters <- unique(info$cluster)

  # p for one cluster: percentage x impact score, renormalized to sum 1, in the
  # table's own row order so it lines up with any distance submatrix.
  cluster_p <- function(cl) {
    d <- info[cluster == cl]
    w <- d$percentage * d$celltype_impact_score
    list(celltypes = as.character(d$celltype), p = w / sum(w))
  }

  run_metric <- function(metric, value_of) {
    vals <- vapply(clusters, function(cl) {
      cp <- cluster_p(cl); value_of(cp$celltypes, cp$p, cl)
    }, numeric(1))
    res <- data.table::data.table(cluster = clusters, v = vals)
    data.table::setnames(res, "v", paste0(metric, "_entropy"))
    out[[paste0("cluster_entropy_", metric)]] <<- res
  }

  if ("Shannon" %in% entropy) run_metric("Shannon", function(ct, p, cl) normalized_shannon(p))
  if ("Gini" %in% entropy)    run_metric("Gini",    function(ct, p, cl) raw_gini(p))

  ont_req <- intersect(ACCORDION_ONTOLOGY_METRICS, entropy)
  if (length(ont_req)) {
    all_ct <- sort(unique(as.character(info$celltype)))
    # Precondition, not a per-pair fallback: a cell type with no ontology term
    # cannot be placed at any distance, so the ontology metrics are refused for
    # the whole analysis rather than computed on a reduced set.
    unmapped <- unmapped_celltypes(all_ct, ontology)
    if (length(unmapped)) {
      if (verbose) message(skip_message(ont_req, unmapped, "have no Cell Ontology term"))
    } else {
      for (metric in ont_req) {
        D <- ontology_distance_matrix(all_ct, ontology, metric, missing = "nan")
        run_metric(metric, function(ct, p, cl) rao(p, D[ct, ct, drop = FALSE]))
      }
    }
  }

  mk_req <- intersect(ACCORDION_MARKER_METRICS, entropy)
  if (length(mk_req)) {
    # Marker sets per cluster: for the cluster being scored, each candidate's
    # markers with a positive impact score *in that cluster*. That is what
    # Jaccard and Overlap are meant to compare -- how distinguishable two
    # candidates are by the evidence present here -- and unlike picking one
    # cluster's scores for the whole run it is well defined.
    mk <- data.table::as.data.table(dt_all_marker)
    cl_col <- intersect(c("cluster", "seurat_clusters"), names(mk))
    sets_by_cluster <- list()
    if (length(cl_col)) {
      cl_col <- cl_col[1]
      for (cl in unique(as.character(mk[[cl_col]]))) {
        sets_by_cluster[[cl]] <- marker_sets(mk[as.character(get(cl_col)) == cl], "celltype")
      }
      sets_for <- function(cl) {
        s <- sets_by_cluster[[as.character(cl)]]
        if (is.null(s)) list() else s
      }
    } else {
      s_all <- marker_sets(mk, "celltype")
      sets_for <- function(cl) s_all
    }
    # Same precondition, checked per cluster since the sets are per cluster.
    no_mk <- sort(unique(unlist(lapply(clusters, function(cl)
      celltypes_without_markers(info[cluster == cl]$celltype, sets_for(cl))))))
    if (length(no_mk)) {
      if (verbose) message(skip_message(mk_req, no_mk,
                                        "have no markers left in the marker table"))
    } else {
      for (metric in mk_req) {
        local({
          met <- metric
          run_metric(met, function(ct, p, cl)
            rao(p, marker_distance_matrix(ct, sets_for(cl), met, missing = "nan")))
        })
      }
    }
  }
  out
}

# --- driver: entropie per cellula ---

#' Cell-level entropies
#'
#' Port of `calculate_entropies_from_score_matrix` in the Python package's
#' `spatial_entropy.py`, which is the path `annotation_resolution = "cell"` takes.
#' `final_dt` is the long score table, one row per (cell, candidate cell type)
#' with columns `cell`, `CL_celltype` and `diff_score`.
#'
#' Two things differ from the cluster driver, on purpose:
#'
#' * `p` is not percentage x impact score, which has no meaning for a single
#'   cell. Each cell keeps the candidates scoring above zero; when there are more
#'   than `n_top_celltypes` of them the set is reduced to those at or above the
#'   `top_cell_score_quantile_threshold` quantile *or* among the top N, and `p` is
#'   proportional to the **square** of the score. The square comes from the
#'   original formula, significance = score / sum(score) followed by
#'   p = significance * score renormalized, and is kept as such.
#' * the distance matrices are global -- one per metric, over every candidate
#'   cell type in the run -- and a pair that cannot be resolved is filled with 1
#'   for the marker metrics and 0 for the ontology ones, which is what Python
#'   does here. The preconditions are still checked over all cell types first, so
#'   the fill cannot actually fire.
#'
#' @return named list of data.tables, `cell_entropy_<metric>`
#' @keywords internal
compute_cell_entropies <- function(final_dt, cells, entropy, ontology = NULL,
                                   accordion_marker = NULL,
                                   top_cell_score_quantile_threshold = 0.90,
                                   n_top_celltypes = 5,
                                   allow_unknown = TRUE,
                                   cell_block_size = 5000,
                                   verbose = TRUE) {
  if (identical(entropy, "all")) entropy <- ACCORDION_METRICS
  entropy <- intersect(entropy, ACCORDION_METRICS)
  if (!length(entropy)) return(list())

  dt <- data.table::as.data.table(final_dt)
  if (!nrow(dt)) return(list())
  cells <- as.character(cells)

  # The healthy database calls the column CL_celltype and the disease one
  # NCIT_celltype; the metrics do not care which, only the preconditions do, and
  # those are checked against the labels themselves further down.
  ct_col <- intersect(c("CL_celltype", "NCIT_celltype", "celltype"), names(dt))
  if (!length(ct_col)) {
    stop("final_dt has no cell type column (CL_celltype, NCIT_celltype or celltype).")
  }
  if (ct_col[1] != "CL_celltype") data.table::setnames(dt, ct_col[1], "CL_celltype")

  # cells x celltypes score matrix, rows in the data's own cell order.
  wide <- data.table::dcast(dt, cell ~ CL_celltype, value.var = "diff_score")
  M <- as.matrix(wide[, -1L, with = FALSE])
  storage.mode(M) <- "double"
  pos <- match(cells, as.character(wide$cell))
  M <- M[pos, , drop = FALSE]
  M[is.na(pos), ] <- NA_real_
  celltypes <- colnames(M)
  n_ct <- ncol(M)
  n_top_celltypes <- max(1L, as.integer(n_top_celltypes))
  if (!n_ct) return(list())

  # Preconditions, checked once over every candidate cell type in the run: a
  # metric whose distances cannot all be defined is refused for the whole
  # analysis, it is not computed on a reduced set.
  dmats <- list()
  ont_req <- intersect(ACCORDION_ONTOLOGY_METRICS, entropy)
  if (length(ont_req)) {
    unmapped <- unmapped_celltypes(celltypes, ontology)
    if (length(unmapped)) {
      if (verbose) message(skip_message(ont_req, unmapped, "have no Cell Ontology term"))
      entropy <- setdiff(entropy, ont_req)
    } else {
      for (met in ont_req) {
        dmats[[met]] <- ontology_distance_matrix(celltypes, ontology, met,
                                                 missing = "zero")
      }
    }
  }

  mk_req <- intersect(ACCORDION_MARKER_METRICS, entropy)
  if (length(mk_req)) {
    # One marker set per cell type over the whole marker database: at cell
    # resolution there is no cluster whose impact scores could restrict them.
    sets <- list()
    if (!is.null(accordion_marker) && nrow(accordion_marker)) {
      mk <- data.table::as.data.table(accordion_marker)
      mk_ct <- intersect(c("CL_celltype", "NCIT_celltype", "celltype"), names(mk))
      if (length(mk_ct)) {
        mk <- unique(mk[, c(mk_ct[1], "marker", "marker_type"), with = FALSE])
        sets <- marker_sets(mk, mk_ct[1])
      }
    }
    no_mk <- celltypes_without_markers(celltypes, sets)
    if (length(no_mk)) {
      if (verbose) message(skip_message(mk_req, no_mk,
                                        "have no markers left in the marker table"))
      entropy <- setdiff(entropy, mk_req)
    } else {
      for (met in mk_req) {
        dmats[[met]] <- marker_distance_matrix(celltypes, sets, met, missing = "max")
      }
    }
  }
  if (!length(entropy)) return(list())

  # p for a block of cells, following the Python fast path step by step.
  prob_block <- function(X) {
    if (allow_unknown) {
      base <- is.finite(X) & X > 0
    } else {
      base <- is.finite(X)
      X[base & X < 0] <- 0
    }
    base[is.na(base)] <- FALSE
    counts <- rowSums(base)
    selected <- base

    red <- which(counts > n_top_celltypes)
    if (length(red)) {
      # Only the rows that need reducing. A cell with no candidate at all -- the
      # ones annotated "unknown" -- has an all-FALSE mask, and asking for a
      # quantile of nothing is both meaningless and, in Python, the source of the
      # "All-NaN slice encountered" warning.
      sm <- base[red, , drop = FALSE]
      sc <- X[red, , drop = FALSE]
      v <- sc
      v[!sm] <- NA_real_
      # Type 7, the default of both stats::quantile and numpy's "linear".
      q <- vapply(seq_len(nrow(v)), function(i)
        stats::quantile(v[i, ], probs = top_cell_score_quantile_threshold,
                        na.rm = TRUE, names = FALSE), numeric(1))
      # q has one entry per row, so recycling it against the matrix compares each
      # element with its own row's quantile.
      keep_q <- sm & sc >= q
      keep_q[is.na(keep_q)] <- FALSE

      top_n <- min(n_top_celltypes, n_ct)
      rk <- sc
      rk[!sm] <- -Inf
      idx <- vapply(seq_len(nrow(rk)), function(i)
        order(rk[i, ], decreasing = TRUE)[seq_len(top_n)], integer(top_n))
      keep_top <- matrix(FALSE, nrow(sm), ncol(sm))
      keep_top[cbind(rep(seq_len(nrow(rk)), each = top_n), as.vector(idx))] <- TRUE
      keep_top <- keep_top & sm

      selected[red, ] <- keep_q | keep_top
    }

    W <- X
    W[!selected] <- 0
    W[is.na(W)] <- 0
    W[W < 0] <- 0
    W <- W * W
    den <- rowSums(W)
    P <- matrix(0, nrow(X), ncol(X))
    ok <- is.finite(den) & den > 0
    if (any(ok)) P[ok, ] <- W[ok, , drop = FALSE] / den[ok]
    list(P = P, n = rowSums(P > 0))
  }

  metric_block <- function(metric, P, n) {
    if (metric == "Shannon") {
      # Normalized by log2 of the number of contributing cell types, which is
      # per cell here and not a property of the run.
      H <- -rowSums(ifelse(P > 0, P * log2(P), 0))
      v <- ifelse(n > 1, H / log2(n), 0)
      v[!is.finite(v)] <- 0
      return(pmax(0, v))
    }
    if (metric == "Gini") {
      # Raw 1 - sum(p^2), as at cluster level: already bounded by 1 for any n.
      v <- ifelse(n > 1, 1 - rowSums(P * P), 0)
      v[!is.finite(v)] <- 0
      return(pmin(1, pmax(0, v)))
    }
    D <- dmats[[metric]]
    if (is.null(D)) return(rep(0, nrow(P)))
    # Rao's Q row by row: diag(P D P') == rowSums((P %*% D) * P).
    v <- rowSums((P %*% D) * P)
    v[!is.finite(v)] <- 0
    pmax(0, v)
  }

  cell_block_size <- max(1L, as.integer(cell_block_size))
  starts <- seq(1L, length(cells), by = cell_block_size)
  vals <- lapply(entropy, function(met) numeric(length(cells)))
  names(vals) <- entropy
  for (s in starts) {
    e <- min(s + cell_block_size - 1L, length(cells))
    pb <- prob_block(M[s:e, , drop = FALSE])
    for (met in entropy) vals[[met]][s:e] <- metric_block(met, pb$P, pb$n)
  }

  out <- list()
  for (met in entropy) {
    res <- data.table::data.table(cell = cells, v = vals[[met]])
    data.table::setnames(res, "v", paste0(met, "_entropy"))
    out[[paste0("cell_entropy_", met)]] <- res
  }
  out
}
