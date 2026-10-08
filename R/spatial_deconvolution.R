#' Spatial deconvolution and spot-level entropy
#'
#' Port of the Python package's `accordion_spatial_deconvolution` and
#' `plot_spatial_entropy`. Python is the reference here, as it is for the entropy
#' metrics: the deconvolution exists only there, so these follow it step by step.
#'
#' @keywords internal
#' @name spatial_deconvolution
NULL

#' Keep only the `n` largest entries of each row, zeroing the rest.
#'
#' @keywords internal
.accordion_keep_top_n <- function(P, n) {
  n <- min(as.integer(n), ncol(P))
  if (n < 1L || !nrow(P)) return(P)
  idx <- vapply(seq_len(nrow(P)),
                function(i) order(P[i, ], decreasing = TRUE)[seq_len(n)],
                integer(n))
  keep <- matrix(FALSE, nrow(P), ncol(P))
  keep[cbind(rep(seq_len(nrow(P)), each = n), as.vector(idx))] <- TRUE
  P * keep
}

#' Per-spot cell type composition
#'
#' `final_dt` is the long score table, one row per (spot, candidate cell type)
#' with columns `cell`, `CL_celltype` and `diff_score`: the same table the
#' annotation and the cell-level entropies are built from.
#'
#' The steps and their order are the ones in Python, and the order matters: the
#' threshold is applied **before** the top-N cap, so a spot can end up with fewer
#' than `n_top_celltypes` contributions, and the renormalisation comes last, so
#' what is reported are fractions of the kept signal and not of the original.
#'
#' @return a spots x celltypes matrix whose rows sum to 1 (a row with no positive
#'   signal at all stays all zero)
#' @keywords internal
.accordion_deconvolution_proportions <- function(final_dt, cells,
                                                 deconvolution_threshold = 0.05,
                                                 n_top_celltypes = 5,
                                                 deconvolution_empty_fallback = FALSE) {
  dt <- data.table::as.data.table(final_dt)
  cells <- as.character(cells)
  wide <- data.table::dcast(dt, cell ~ CL_celltype, value.var = "diff_score")
  M <- as.matrix(wide[, -1L, with = FALSE])
  storage.mode(M) <- "double"
  pos <- match(cells, as.character(wide$cell))
  M <- M[pos, , drop = FALSE]
  M[is.na(pos), ] <- 0
  dimnames(M) <- list(cells, colnames(M))

  # Only positive evidence contributes to a composition: a negative score means
  # the spot argues *against* that cell type, not that it holds a negative
  # fraction of it.
  M[!is.finite(M) | M < 0] <- 0
  rs <- rowSums(M)
  P <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
  ok <- rs > 0
  if (any(ok)) P[ok, ] <- M[ok, , drop = FALSE] / rs[ok]
  # The pre-threshold fractions, kept for the fallback below.
  P_norm <- P

  if (!is.null(deconvolution_threshold) && deconvolution_threshold) {
    P[P < as.numeric(deconvolution_threshold)] <- 0
  }

  n_cap <- if (is.null(n_top_celltypes)) 0L else as.integer(n_top_celltypes)
  # A spot is a mixture of a few cell types, not of every candidate.
  if (n_cap > 0L && ncol(P) > n_cap) P <- .accordion_keep_top_n(P, n_cap)

  # A spot whose signal is spread thin over many cell types has no single
  # contribution above the threshold and would come out empty despite having
  # positive signal. With the fallback on, it keeps its top types from the
  # pre-threshold fractions instead.
  emptied <- rowSums(P) == 0 & rs > 0
  if (isTRUE(deconvolution_empty_fallback) && any(emptied)) {
    n_fb <- if (n_cap > 0L) min(n_cap, ncol(P_norm)) else ncol(P_norm)
    P[emptied, ] <- .accordion_keep_top_n(P_norm[emptied, , drop = FALSE], n_fb)
  }

  rs2 <- rowSums(P)
  ok2 <- rs2 > 0
  if (any(ok2)) P[ok2, ] <- P[ok2, , drop = FALSE] / rs2[ok2]
  P[!ok2, ] <- 0
  P
}

#' Run the Accordion spatial deconvolution
#'
#' Thin wrapper around [accordion()] for spatial data, where each observation is
#' a spot to be deconvolved into candidate cell types rather than a single cell.
#' It is the R equivalent of the Python `accordion_spatial_deconvolution`, and
#' calls `accordion(..., annotation_resolution = "cell", do_deconvolution = TRUE)`.
#'
#' @param data Seurat object (or raw counts matrix) of spatial data.
#' @param entropy Entropy metric(s) to compute per spot, or "all". Default
#'   "Overlap-Rao". NULL or FALSE skips the entropy.
#' @param deconvolution_threshold Contributions below this fraction are dropped
#'   from a spot's composition. Default 0.05.
#' @param min_percentage_celltype_entropy Passed through to [accordion()].
#' @param annotation_name Prefix of the output columns and of the proportions
#'   assay. Default "accordion".
#' @param ... Further arguments for [accordion()], notably `assay` (often
#'   "Spatial" for Visium), `tissue`, `species`, `n_top_celltypes` and
#'   `deconvolution_empty_fallback`.
#'
#' @return The input object with, for a Seurat object:
#'   \describe{
#'     \item{\code{<annotation_name>_per_cell}}{assigned cell type per spot.}
#'     \item{\code{<annotation_name>_per_cell_score}}{its score.}
#'     \item{\code{<annotation_name>_per_cell_percentage}}{percentage (0-100) of
#'       the assigned cell type in that spot.}
#'     \item{assay \code{<annotation_name>_proportions}}{cell types x spots, each
#'       spot's column summing to 1. Plots straight out of
#'       \code{Seurat::SpatialFeaturePlot(data, features = "<a cell type>")}
#'       after \code{DefaultAssay(data) <- "<annotation_name>_proportions"}.}
#'     \item{\code{<annotation_name>_<metric>_entropy}}{one column per entropy
#'       metric requested.}
#'   }
#'   With a raw matrix in input the same tables come back in the output list.
#'
#' @export
accordion_spatial_deconvolution <- function(data,
                                            entropy = "Overlap-Rao",
                                            deconvolution_threshold = 0.05,
                                            min_percentage_celltype_entropy = 0.05,
                                            annotation_name = "accordion",
                                            ...) {
  accordion(data,
            annotation_resolution = "cell",
            do_deconvolution = TRUE,
            entropy = entropy,
            deconvolution_threshold = deconvolution_threshold,
            min_percentage_celltype_entropy = min_percentage_celltype_entropy,
            annotation_name = annotation_name,
            ...)
}

#: Colormap per metric, the same ones the Python plot uses. They are
#: RColorBrewer palette names, which is what ggplot2's *_distiller scales take.
ACCORDION_ENTROPY_PALETTES <- c(
  "Shannon" = "Greens",
  "Gini" = "Oranges",
  "Ontology-sp-Rao" = "Reds",
  "Ontology-Lin-Rao" = "Purples",
  "Jaccard-Rao" = "Blues",
  "Overlap-Rao" = "YlOrBr"
)

#' Spot coordinates for the spatial plot
#'
#' @keywords internal
.accordion_spot_coordinates <- function(data, coordinates = NULL) {
  if (!is.null(coordinates)) {
    co <- as.data.frame(coordinates)
    if (ncol(co) < 2) stop("coordinates must have at least two columns (x and y)")
    out <- data.frame(x = as.numeric(co[[1]]), y = as.numeric(co[[2]]))
    rownames(out) <- rownames(co)
    return(out)
  }
  co <- NULL
  if (length(methods::slotNames(data)) && length(data@images)) {
    co <- try(Seurat::GetTissueCoordinates(data), silent = TRUE)
    if (inherits(co, "try-error")) co <- NULL
  }
  if (is.null(co)) {
    stop("No spatial coordinates found. Pass them explicitly with ",
         "`coordinates = <data.frame with x and y and the spots as row names>`.")
  }
  co <- as.data.frame(co)
  # Seurat v5 returns x/y/cell, v4 imagerow/imagecol.
  nm <- tolower(names(co))
  xi <- if ("x" %in% nm) which(nm == "x") else which(nm == "imagecol")
  yi <- if ("y" %in% nm) which(nm == "y") else which(nm == "imagerow")
  if (!length(xi) || !length(yi)) { xi <- 1L; yi <- 2L }
  rn <- if ("cell" %in% nm) as.character(co[[which(nm == "cell")]]) else rownames(co)
  out <- data.frame(x = as.numeric(co[[xi[1]]]), y = as.numeric(co[[yi[1]]]))
  rownames(out) <- rn
  out
}

#' Plot the spot-level entropy on the tissue
#'
#' R equivalent of the Python `plot_spatial_entropy`: the spots drawn at their
#' tissue coordinates and coloured by entropy, next to the distribution of the
#' same values. Run [accordion_spatial_deconvolution()] (or [accordion()] with
#' `entropy =` at cell resolution) first.
#'
#' @param data Seurat object annotated by the Accordion.
#' @param entropy_type Which metric to draw, e.g. "Shannon", "Gini",
#'   "Overlap-Rao". Default "Overlap-Rao".
#' @param annotation_name The prefix used when running the Accordion.
#' @param coordinates Optional data.frame of spot coordinates (first two columns
#'   x and y, spots as row names). By default they are read from the object's
#'   image slot, which is where Visium data keeps them.
#' @param palette RColorBrewer palette name. NULL (the default) picks the
#'   metric's own: Greens for Shannon, Oranges for Gini, Reds for
#'   Ontology-sp-Rao, Purples for Ontology-Lin-Rao, Blues for Jaccard-Rao,
#'   YlOrBr for Overlap-Rao.
#' @param spot_size Point size. Raise it on sparse tissues, lower it on dense
#'   ones. Default 1.5.
#' @param title,xlabel,ylabel,legend_label Text of the spatial panel.
#' @param hist_title,hist_xlabel,hist_ylabel Text of the distribution panel.
#' @param bins Histogram bins. Default 30.
#' @param save Optional file path to write the figure to.
#' @param width,height,dpi Size and resolution used when saving.
#'
#' @return A list with the assembled plot, the two panels and the plotted data.
#' @export
plot_spatial_entropy <- function(data,
                                 entropy_type = "Overlap-Rao",
                                 annotation_name = "accordion",
                                 coordinates = NULL,
                                 palette = NULL,
                                 spot_size = 1.5,
                                 title = NULL,
                                 xlabel = "Spatial X",
                                 ylabel = "Spatial Y",
                                 legend_label = NULL,
                                 hist_title = "Distribution",
                                 hist_xlabel = NULL,
                                 hist_ylabel = "Frequency",
                                 bins = 30,
                                 save = NULL,
                                 width = 12,
                                 height = 5,
                                 dpi = 300) {
  if (is.null(palette)) {
    # Single brackets, not double: an unknown metric should fall back to the
    # default palette, not error out with "subscript out of bounds".
    palette <- unname(ACCORDION_ENTROPY_PALETTES[entropy_type])
    if (is.na(palette)) palette <- "Blues"
  }
  ent_col <- paste0(annotation_name, "_", entropy_type, "_entropy")
  md <- data@meta.data
  if (!(ent_col %in% colnames(md))) {
    avail <- grep("_entropy$", colnames(md), value = TRUE)
    stop("Column '", ent_col, "' not found in the metadata. ",
         if (length(avail)) paste0("Available: ", paste(avail, collapse = ", "), ".")
         else "Run the annotation with entropy enabled at cell resolution first.")
  }

  co <- .accordion_spot_coordinates(data, coordinates)
  common <- intersect(rownames(co), rownames(md))
  if (!length(common)) stop("The coordinates and the metadata share no spot.")
  pd <- data.frame(cell = common,
                   x = co[common, "x"], y = co[common, "y"],
                   entropy = as.numeric(md[common, ent_col]),
                   stringsAsFactors = FALSE)
  pd$entropy[!is.finite(pd$entropy)] <- 0

  resolved_title <- if (is.null(title)) paste(entropy_type, "entropy") else title
  resolved_legend <- if (is.null(legend_label)) resolved_title else legend_label
  resolved_hist_x <- if (is.null(hist_xlabel)) resolved_title else hist_xlabel
  # The colour scale is fixed to [0, 1], the range of every metric, so two
  # tissues or two metrics stay comparable instead of each being stretched over
  # its own range.
  p_space <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = y, colour = entropy)) +
    ggplot2::geom_point(size = spot_size) +
    ggplot2::scale_colour_distiller(palette = palette, direction = 1,
                                    limits = c(0, 1), name = resolved_legend,
                                    guide = ggplot2::guide_colourbar(
                                      barheight = grid::unit(0.45, "npc"))) +
    # Image coordinates grow downwards, so the y axis is reversed to show the
    # tissue the right way up.
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = resolved_title, x = xlabel, y = ylabel) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15),
                   axis.title = ggplot2::element_text(face = "bold"),
                   # Framed like the Python panel, where the tissue sits inside a
                   # box rather than against two bare axes.
                   panel.border = ggplot2::element_rect(colour = "black",
                                                        fill = NA, linewidth = 0.6),
                   legend.title = ggplot2::element_text(face = "bold", size = 11))

  hist_color <- RColorBrewer::brewer.pal(9, palette)[7]
  mean_val <- mean(pd$entropy, na.rm = TRUE)
  p_hist <- ggplot2::ggplot(pd, ggplot2::aes(x = entropy)) +
    ggplot2::geom_histogram(bins = bins, fill = hist_color, colour = "white",
                            linewidth = 0.3, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = mean_val, colour = "#8B1A1A",
                        linetype = "dashed", linewidth = 0.7) +
    # The label goes on whichever side of the mean line has room: with a
    # left-skewed metric like Shannon the mean sits near the right edge and the
    # text would otherwise be cut off.
    ggplot2::annotate("text", x = mean_val, y = Inf,
                      label = paste0("Mean: ", format(round(mean_val, 3), nsmall = 3)),
                      hjust = if (mean_val > mean(range(pd$entropy, na.rm = TRUE))) 1.1 else -0.1,
                      vjust = 1.8, size = 3.6, colour = "#8B1A1A") +
    ggplot2::labs(title = hist_title, x = resolved_hist_x, y = hist_ylabel) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
                   axis.title = ggplot2::element_text(face = "bold"))

  combined <- cowplot::plot_grid(p_space, p_hist, ncol = 2, rel_widths = c(3, 1.4))
  if (!is.null(save)) {
    ggplot2::ggsave(save, combined, width = width, height = height, dpi = dpi)
  }
  list(plot = combined, spatial = p_space, distribution = p_hist, data = pd)
}
