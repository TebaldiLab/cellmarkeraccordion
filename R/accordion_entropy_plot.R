#' Where the entropy of a run lives, and at which resolution
#'
#' At cell resolution it is a metadata column; at cluster resolution it is a
#' table in `misc`, one row per cluster, which is spread back over the cells so
#' that both can be drawn on the same embedding.
#'
#' @keywords internal
.accordion_entropy_values <- function(data, entropy_type, annotation_name) {
  md <- data@meta.data
  cell_col <- paste0(annotation_name, "_", entropy_type, "_entropy")
  if (cell_col %in% colnames(md)) {
    return(list(resolution = "cell",
                values = stats::setNames(as.numeric(md[[cell_col]]), rownames(md)),
                cluster_column = NULL))
  }

  info <- data@misc[[annotation_name]][["cluster_resolution"]][["detailed_annotation_info"]]
  tab <- info[[paste0("cluster_entropy_", entropy_type)]]
  if (!is.null(tab)) {
    tab <- as.data.frame(tab)
    cl_col <- names(tab)[1]
    val_col <- paste0(entropy_type, "_entropy")
    if (!(cl_col %in% colnames(md))) {
      stop("The cluster entropy is keyed by '", cl_col,
           "', which is not a column of the metadata.")
    }
    v <- tab[[val_col]][match(as.character(md[[cl_col]]), as.character(tab[[cl_col]]))]
    return(list(resolution = "cluster",
                values = stats::setNames(as.numeric(v), rownames(md)),
                cluster_column = cl_col))
  }

  avail_cell <- sub("_entropy$", "",
                    sub(paste0("^", annotation_name, "_"), "",
                        grep(paste0("^", annotation_name, "_.*_entropy$"),
                             colnames(md), value = TRUE)))
  avail_cl <- sub("^cluster_entropy_", "",
                  grep("^cluster_entropy_", names(info), value = TRUE))
  avail <- unique(c(avail_cell, avail_cl))
  stop("No '", entropy_type, "' entropy found for '", annotation_name, "'. ",
       if (length(avail)) paste0("Available: ", paste(avail, collapse = ", "), ".")
       else paste0("Run the annotation with entropy = \"all\" (or a specific ",
                   "metric) first."))
}

#' Which embedding to draw on
#'
#' @keywords internal
.accordion_pick_reduction <- function(data, reduction = NULL) {
  have <- tryCatch(SeuratObject::Reductions(data), error = function(e) character(0))
  if (!is.null(reduction)) {
    if (!(reduction %in% have)) {
      stop("Reduction '", reduction, "' not found. ",
           if (length(have)) paste0("Available: ", paste(have, collapse = ", "), ".")
           else "The object has none: run RunPCA/RunUMAP first.")
    }
    return(reduction)
  }
  for (r in c("umap", "tsne", "pca")) if (r %in% have) return(r)
  if (length(have)) return(have[1])
  NULL
}

#' Plot the annotation entropy of a single-cell dataset
#'
#' The companion of [plot_spatial_entropy()] for data that has no tissue
#' coordinates: the cells on their embedding, coloured by entropy, next to how
#' those values are distributed. It reads both resolutions -- the cell-level
#' entropy from the metadata, the cluster-level one from
#' \code{misc[[annotation_name]]} -- so it works whichever way the annotation was
#' run.
#'
#' @param data Seurat object annotated by the Accordion with `entropy` enabled.
#' @param entropy_type Which metric to draw: "Shannon", "Gini",
#'   "Ontology-sp-Rao", "Ontology-Lin-Rao", "Jaccard-Rao", "Overlap-Rao".
#' @param annotation_name The prefix used when running the Accordion.
#' @param reduction Embedding to draw on. NULL (the default) takes the first of
#'   umap, tsne, pca that the object has.
#' @param group_by Metadata column to split the right-hand panel by, typically
#'   the clusters or the Accordion annotation. NULL (the default) draws a
#'   histogram at cell resolution and one bar per cluster at cluster resolution.
#' @param palette RColorBrewer palette. NULL picks the metric's own, the same
#'   ones [plot_spatial_entropy()] uses.
#' @param point_size Point size on the embedding. Default 0.6.
#' @param max_groups With `group_by`, at most this many groups are drawn, the
#'   largest ones. Default 30.
#' @param title,legend_label Text of the embedding panel.
#' @param dist_title,dist_xlabel,dist_ylabel Text of the distribution panel.
#' @param bins Histogram bins. Default 30.
#' @param save Optional file path to write the figure to.
#' @param width,height,dpi Size and resolution used when saving.
#'
#' @return A list with the assembled plot, the two panels separately and the
#'   plotted data.
#' @export
plot_entropy <- function(data,
                         entropy_type = "Shannon",
                         annotation_name = "accordion",
                         reduction = NULL,
                         group_by = NULL,
                         palette = NULL,
                         point_size = 0.6,
                         max_groups = 30,
                         title = NULL,
                         legend_label = NULL,
                         dist_title = NULL,
                         dist_xlabel = NULL,
                         dist_ylabel = NULL,
                         bins = 30,
                         save = NULL,
                         width = 12,
                         height = 5,
                         dpi = 300) {
  if (is.null(palette)) {
    palette <- unname(ACCORDION_ENTROPY_PALETTES[entropy_type])
    if (is.na(palette)) palette <- "Blues"
  }
  ent <- .accordion_entropy_values(data, entropy_type, annotation_name)
  md <- data@meta.data
  pd <- data.frame(cell = rownames(md), entropy = as.numeric(ent$values),
                   stringsAsFactors = FALSE)
  pd$entropy[!is.finite(pd$entropy)] <- NA_real_

  if (is.null(group_by) && ent$resolution == "cluster") group_by <- ent$cluster_column
  if (!is.null(group_by)) {
    if (!(group_by %in% colnames(md))) {
      stop("group_by = '", group_by, "' is not a column of the metadata.")
    }
    pd$group <- as.character(md[[group_by]])
  }

  red <- .accordion_pick_reduction(data, reduction)
  resolved_title <- if (is.null(title)) paste(entropy_type, "entropy") else title
  resolved_legend <- if (is.null(legend_label)) resolved_title else legend_label
  resolved_dist_x <- if (is.null(dist_xlabel)) resolved_title else dist_xlabel

  p_emb <- NULL
  if (!is.null(red)) {
    emb <- as.data.frame(SeuratObject::Embeddings(data, red)[, 1:2, drop = FALSE])
    names(emb) <- c("dim1", "dim2")
    pd$dim1 <- emb[pd$cell, "dim1"]
    pd$dim2 <- emb[pd$cell, "dim2"]
    axis_lab <- toupper(red)
    p_emb <- ggplot2::ggplot(pd[!is.na(pd$entropy), ],
                             ggplot2::aes(x = dim1, y = dim2, colour = entropy)) +
      ggplot2::geom_point(size = point_size) +
      # Fixed to [0, 1], the range every metric lives in, so two datasets or two
      # metrics stay comparable instead of each being stretched over its own range.
      ggplot2::scale_colour_distiller(palette = palette, direction = 1,
                                      limits = c(0, 1), name = resolved_legend,
                                      guide = ggplot2::guide_colourbar(
                                        barheight = grid::unit(0.45, "npc"))) +
      ggplot2::labs(title = resolved_title,
                    x = paste0(axis_lab, " 1"), y = paste0(axis_lab, " 2")) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15),
                     axis.title = ggplot2::element_text(face = "bold"),
                     panel.border = ggplot2::element_rect(colour = "black",
                                                          fill = NA, linewidth = 0.6),
                     legend.title = ggplot2::element_text(face = "bold", size = 11))
  } else {
    warning("The object has no embedding, drawing only the distribution. ",
            "Run RunPCA/RunUMAP first for the left panel.")
  }

  dat <- pd[!is.na(pd$entropy), ]
  hist_color <- RColorBrewer::brewer.pal(9, palette)[7]

  if (is.null(group_by)) {
    mean_val <- mean(dat$entropy)
    p_dist <- ggplot2::ggplot(dat, ggplot2::aes(x = entropy)) +
      ggplot2::geom_histogram(bins = bins, fill = hist_color, colour = "white",
                              linewidth = 0.3, alpha = 0.85) +
      ggplot2::geom_vline(xintercept = mean_val, colour = "#8B1A1A",
                          linetype = "dashed", linewidth = 0.7) +
      ggplot2::annotate("text", x = mean_val, y = Inf,
                        label = paste0("Mean: ", format(round(mean_val, 3), nsmall = 3)),
                        hjust = if (mean_val > mean(range(dat$entropy))) 1.1 else -0.1,
                        vjust = 1.8, size = 3.6, colour = "#8B1A1A") +
      ggplot2::labs(title = if (is.null(dist_title)) "Distribution" else dist_title,
                    x = resolved_dist_x,
                    y = if (is.null(dist_ylabel)) "Frequency" else dist_ylabel)
  } else {
    # Only the largest groups: with a per-cell annotation the tail is made of
    # cell types carrying two or three cells, whose distribution says nothing.
    keep <- names(sort(table(dat$group), decreasing = TRUE))
    if (length(keep) > max_groups) {
      message("Drawing the ", max_groups, " largest groups out of ", length(keep), ".")
      keep <- keep[seq_len(max_groups)]
    }
    dat <- dat[dat$group %in% keep, , drop = FALSE]
    med <- stats::aggregate(entropy ~ group, dat, stats::median)
    dat$group <- factor(dat$group, levels = med$group[order(med$entropy)])
    # At cluster resolution every cell of a cluster carries the same number, so
    # there is nothing to draw a distribution of: one bar per group is the honest
    # shape, and a violin would collapse to a flat tick.
    one_value_per_group <- all(vapply(split(dat$entropy, dat$group),
                                      function(v) length(unique(v)) <= 1,
                                      logical(1)))
    p_dist <- if (one_value_per_group) {
      agg <- med
      agg$group <- factor(agg$group, levels = levels(dat$group))
      ggplot2::ggplot(agg, ggplot2::aes(x = group, y = entropy)) +
        ggplot2::geom_col(ggplot2::aes(fill = entropy), width = 0.75) +
        ggplot2::scale_fill_distiller(palette = palette, direction = 1,
                                      limits = c(0, 1), guide = "none")
    } else {
      ggplot2::ggplot(dat, ggplot2::aes(x = group, y = entropy)) +
        ggplot2::geom_violin(fill = hist_color, colour = NA, alpha = 0.55,
                             scale = "width", trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.25, fill = "white",
                              colour = "grey25", linewidth = 0.35)
    }
    p_dist <- p_dist +
      # Cell type names run long ("naive thymus-derived CD4-positive,
      # alpha-beta T cell"), so they are wrapped rather than left to push the
      # panel out of the figure.
      ggplot2::scale_x_discrete(labels = function(z) stringr::str_wrap(z, width = 28)) +
      ggplot2::coord_flip() +
      ggplot2::ylim(0, max(1, max(dat$entropy, na.rm = TRUE))) +
      ggplot2::labs(title = if (is.null(dist_title)) paste("By", group_by) else dist_title,
                    x = NULL,
                    y = if (is.null(dist_ylabel)) resolved_dist_x else dist_ylabel)
  }

  p_dist <- p_dist +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
                   axis.title = ggplot2::element_text(face = "bold"))

  # The right panel gets more room when its labels are long, so that wrapping
  # them does not squeeze the plotting area to nothing.
  w_dist <- 2
  if (!is.null(group_by)) {
    lab_n <- suppressWarnings(max(nchar(as.character(dat$group)), 0))
    w_dist <- 2 + min(1.6, lab_n / 28)
  }
  combined <- if (is.null(p_emb)) p_dist else
    cowplot::plot_grid(p_emb, p_dist, ncol = 2, rel_widths = c(3, w_dist))
  if (!is.null(save)) {
    ggplot2::ggsave(save, combined, width = width, height = height, dpi = dpi)
  }
  list(plot = combined, embedding = p_emb, distribution = p_dist, data = pd,
       resolution = ent$resolution)
}
