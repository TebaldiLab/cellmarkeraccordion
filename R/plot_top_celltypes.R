#' The detailed annotation tables of a run, whatever the input was
#'
#' @keywords internal
.accordion_detailed_info <- function(data, annotation_name, resolution = "cluster") {
  slot <- paste0(resolution, "_resolution")
  info <- if (inherits(data, "Seurat")) {
    data@misc[[annotation_name]][[slot]][["detailed_annotation_info"]]
  } else {
    data[[annotation_name]][[slot]][["detailed_annotation_info"]]
  }
  if (is.null(info) || !length(info)) {
    stop("No detailed annotation info for '", annotation_name, "' at ", resolution,
         " resolution. Run the annotation with ",
         "include_detailed_annotation_info = TRUE.")
  }
  info
}

#' First of `candidates` that is a column of `df`
#'
#' @keywords internal
.accordion_first_col <- function(df, candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[1] else NULL
}

#: Rows of the entropy block, in the order they are drawn, with the palette each
#: metric keeps across every plot in the package.
.ACCORDION_ENTROPY_ROWS <- c("Shannon", "Gini", "Ontology-sp-Rao",
                             "Ontology-Lin-Rao", "Jaccard-Rao", "Overlap-Rao")

#' Dot plot of the top cell types per cluster, with the entropy underneath
#'
#' The R version of the Python package's `plot_top_celltypes_per_cluster`.
#' Clusters on the x axis, candidate cell types on the y, dot area growing with
#' the impact score, and the cell type that wins a cluster ringed in black. Under
#' the plot, the size of each cluster and one heatmap row per entropy metric
#' found in the annotation, so the question "how contested was this label?" is
#' read off the same columns as the labels themselves.
#'
#' Only the entropy metrics actually present are drawn, so the figure adapts to
#' having all of them, some, or none.
#'
#' @param data Seurat object (or the list returned for a matrix input) annotated
#'   with `include_detailed_annotation_info = TRUE`.
#' @param annotation_name The prefix used when running the Accordion.
#' @param resolution "cluster" (the default). The cell-resolution tables have no
#'   per-cluster entropy, so this is the resolution the plot is about.
#' @param cluster_colors Named vector of colours, one per cluster. NULL (the
#'   default) uses the same hue palette Seurat gives `DimPlot`, so the columns
#'   match the clusters on the UMAP.
#' @param show_entropy Logical, draw the entropy rows. Default TRUE.
#' @param entropy_types Which metrics to draw, in order. NULL takes every one
#'   present.
#' @param show_entropy_values Logical, write the value inside each tile.
#' @param show_n_cells Logical, draw the "Number of cells" strip.
#' @param size_range Smallest and largest dot size, in mm. Default c(2, 11).
#' @param legend_quantiles Which quantiles of the impact scores to show in the
#'   size legend. Default c(0.25, 0.5, 0.75, 0.95), as in the Python plot.
#' @param title,xlabel,ylabel,legend_title Text overrides.
#' @param celltype_fontsize,ax_fontsize,value_fontsize Font sizes of the cell
#'   type labels, the axes, and the numbers inside the entropy tiles.
#' @param entropy_row_height Height of one entropy row relative to the dot plot.
#'   Default 0.055.
#' @param save Optional file path to write the figure to.
#' @param width,height,dpi Size and resolution used when saving. NULL sizes the
#'   figure from the number of clusters and cell types.
#'
#' @return A list with the assembled plot, the single panels and the data drawn.
#' @export
plot_top_celltypes_per_cluster <- function(data,
                                           annotation_name = "accordion",
                                           resolution = "cluster",
                                           cluster_colors = NULL,
                                           show_entropy = TRUE,
                                           entropy_types = NULL,
                                           show_entropy_values = TRUE,
                                           show_n_cells = TRUE,
                                           size_range = c(2, 11),
                                           legend_quantiles = c(0.25, 0.5, 0.75, 0.95),
                                           title = "Top cell types per cluster",
                                           xlabel = "Cluster",
                                           ylabel = NULL,
                                           legend_title = "Impact score",
                                           celltype_fontsize = 8,
                                           ax_fontsize = 10,
                                           value_fontsize = 2.6,
                                           entropy_row_height = 0.055,
                                           save = NULL,
                                           width = NULL,
                                           height = NULL,
                                           dpi = 300) {
  info <- .accordion_detailed_info(data, annotation_name, resolution)
  if (is.null(info[["top_celltypes"]])) {
    stop("No 'top_celltypes' table in the detailed annotation info. ",
         "Available: ", paste(names(info), collapse = ", "), ".")
  }
  df <- as.data.frame(info[["top_celltypes"]])

  score_col <- "celltype_impact_score"
  if (!(score_col %in% names(df))) {
    stop("Column '", score_col, "' not found in top_celltypes. ",
         "Available: ", paste(names(df), collapse = ", "), ".")
  }
  ct_col <- .accordion_first_col(df, c(paste0(annotation_name, "_per_cluster"),
                                       annotation_name, "annotation_per_cluster",
                                       "CL_celltype", "celltype"))
  cl_col <- .accordion_first_col(df, setdiff(names(df),
                                             c(ct_col, score_col, "ncell_tot_cluster",
                                               "perc_celltype_cluster", "percentage")))
  if (is.null(ct_col) || is.null(cl_col)) {
    stop("Could not tell which columns hold the cluster and the cell type. ",
         "Available: ", paste(names(df), collapse = ", "), ".")
  }

  df <- df[!is.na(df[[cl_col]]) & !is.na(df[[ct_col]]) & !is.na(df[[score_col]]), ,
           drop = FALSE]
  df$cluster <- as.character(df[[cl_col]])
  df$celltype <- as.character(df[[ct_col]])
  df$score <- as.numeric(df[[score_col]])

  # Clusters left to right: numerically when they are numbers, which is what
  # anyone reading "0, 1, 2, 10" expects, and alphabetically otherwise.
  clusters <- unique(df$cluster)
  num <- suppressWarnings(as.numeric(clusters))
  clusters <- if (!anyNA(num)) clusters[order(num)] else sort(clusters)

  # Rows: each cell type sits in the cluster where it scores highest, and inside
  # that block the order is by score, so a cluster's winner heads its own block
  # instead of being buried in it.
  best <- do.call(rbind, lapply(split(df, df$celltype), function(d) d[which.max(d$score), ]))
  best$cx <- match(best$cluster, clusters)
  best <- best[order(best$cx, -best$score), ]
  celltypes <- best$celltype

  df$cluster <- factor(df$cluster, levels = clusters)
  df$celltype <- factor(df$celltype, levels = rev(celltypes))
  df$winner <- unsplit(lapply(split(df$score, df$cluster, drop = TRUE),
                              function(s) s == max(s)), df$cluster, drop = TRUE)

  if (is.null(cluster_colors)) {
    cluster_colors <- stats::setNames(scales::hue_pal()(length(clusters)), clusters)
  } else {
    miss <- setdiff(clusters, names(cluster_colors))
    if (length(miss)) {
      cluster_colors <- c(cluster_colors,
                          stats::setNames(scales::hue_pal()(length(miss)), miss))
    }
  }

  # A negative impact score means the cluster argues against that cell type; it
  # has no size to draw, and feeding it to the sqrt scale would make the point
  # disappear with a NaN warning. Clipped to zero, as the Python plot does.
  df$size_value <- pmax(df$score, 0)

  # Dot area grows with the square root of the score, as in the Python plot: the
  # impact scores span an order of magnitude and a linear area would leave the
  # small ones invisible. ggplot's size scale already squares, so the sqrt
  # transform here is what gives area ~ sqrt(score).
  brk <- unique(round(stats::quantile(df$score[df$score > 0],
                                      probs = legend_quantiles, na.rm = TRUE)))
  # A break rounded past the largest score would fall outside the scale and be
  # silently dropped, taking the biggest dot out of the legend.
  brk <- sort(brk[brk > 0 & brk <= max(df$size_value, na.rm = TRUE)])

  p_main <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = celltype)) +
    ggplot2::geom_point(data = df[!df$winner, , drop = FALSE],
                        ggplot2::aes(size = size_value, colour = cluster), alpha = 0.75) +
    # The winner is the label the cluster actually got: ringed, not recoloured,
    # so it reads as "this one" without leaving the cluster's colour.
    ggplot2::geom_point(data = df[df$winner, , drop = FALSE],
                        ggplot2::aes(size = size_value, fill = cluster),
                        shape = 21, colour = "black", stroke = 1.1, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = cluster_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = cluster_colors, guide = "none") +
    ggplot2::scale_size_continuous(trans = "sqrt", range = size_range,
                                   breaks = brk, name = legend_title) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    # The order is set explicitly, not left to the factor levels: the two point
    # layers draw different subsets of the rows, and letting the scale collect
    # its values from them put all the winners at the top instead of each at the
    # head of its own cluster's block.
    ggplot2::scale_y_discrete(limits = levels(df$celltype), drop = FALSE) +
    ggplot2::labs(title = title, x = xlabel,
                  y = if (is.null(ylabel)) NULL else ylabel) +
    ggplot2::theme_bw(base_size = ax_fontsize) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ax_fontsize * 1.35,
                                         hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      axis.text.y = ggplot2::element_text(size = celltype_fontsize),
      axis.text.x = ggplot2::element_text(size = ax_fontsize),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = ax_fontsize * 0.9),
      legend.key = ggplot2::element_blank())

  # Everything below the dot plot shares its x scale and hides its own, so the
  # strips line up column by column under the clusters.
  strip_theme <- function() {
    ggplot2::theme_minimal(base_size = ax_fontsize) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(angle = 0, hjust = 1, vjust = 0.5,
                                             size = ax_fontsize * 0.85),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(1, 1, 1, 1))
  }

  p_n <- NULL
  if (show_n_cells && "ncell_tot_cluster" %in% names(df)) {
    nd <- unique(df[, c("cluster", "ncell_tot_cluster")])
    nd <- nd[!duplicated(nd$cluster), ]
    p_n <- ggplot2::ggplot(nd, ggplot2::aes(x = cluster, y = 1)) +
      ggplot2::geom_text(ggplot2::aes(label = format(ncell_tot_cluster, big.mark = "")),
                         size = value_fontsize * 1.15) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::labs(y = "Number of cells") +
      strip_theme()
  }

  ent_plots <- list()
  if (show_entropy) {
    present <- sub("^cluster_entropy_", "",
                   grep("^cluster_entropy_", names(info), value = TRUE))
    order_rows <- c(.ACCORDION_ENTROPY_ROWS,
                    setdiff(present, .ACCORDION_ENTROPY_ROWS))
    rows <- if (is.null(entropy_types)) intersect(order_rows, present)
            else intersect(entropy_types, present)
    for (met in rows) {
      tab <- as.data.frame(info[[paste0("cluster_entropy_", met)]])
      val_col <- paste0(met, "_entropy")
      if (!(val_col %in% names(tab))) next
      tcl <- .accordion_first_col(tab, c(cl_col, "cluster", "seurat_clusters",
                                         setdiff(names(tab), val_col)))
      if (is.null(tcl)) next
      ed <- data.frame(cluster = factor(clusters, levels = clusters),
                       value = tab[[val_col]][match(clusters, as.character(tab[[tcl]]))])
      pal <- unname(ACCORDION_ENTROPY_PALETTES[met])
      if (is.na(pal)) pal <- "Blues"
      p <- ggplot2::ggplot(ed, ggplot2::aes(x = cluster, y = 1, fill = value)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4)
      if (show_entropy_values) {
        p <- p + ggplot2::geom_text(
          # White on the dark end of the ramp, black on the light end: the same
          # 0.45 cut the Python plot uses.
          ggplot2::aes(label = formatC(value, format = "f", digits = 2),
                       colour = value > 0.45),
          size = value_fontsize, show.legend = FALSE) +
          ggplot2::scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"),
                                       guide = "none", na.value = "black")
      }
      p <- p +
        # Fixed to [0, 1], the range of every metric, so the rows are comparable
        # with each other and across datasets. The metric's name is the colour
        # bar's own title rather than an axis label, which is what puts it to the
        # left of the bar as in the Python figure.
        ggplot2::scale_fill_distiller(palette = pal, direction = 1,
                                      limits = c(0, 1), breaks = c(0, 1),
                                      name = paste(met, "entropy"),
                                      na.value = "grey92",
                                      guide = ggplot2::guide_colourbar(
                                        direction = "horizontal",
                                        title.position = "left", title.hjust = 1,
                                        title.vjust = 1,
                                        barwidth = grid::unit(2.4, "cm"),
                                        barheight = grid::unit(0.3, "cm"),
                                        ticks = FALSE,
                                        label.position = "bottom")) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        strip_theme() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       legend.position = "left",
                       legend.margin = ggplot2::margin(0, 6, 0, 0),
                       legend.title = ggplot2::element_text(size = ax_fontsize * 0.85),
                       legend.text = ggplot2::element_text(size = ax_fontsize * 0.7))
      ent_plots[[met]] <- p
    }
  }

  panels <- c(list(p_main), if (!is.null(p_n)) list(p_n), ent_plots)
  heights <- c(1, if (!is.null(p_n)) entropy_row_height * 0.8,
               rep(entropy_row_height, length(ent_plots)))
  combined <- cowplot::plot_grid(plotlist = panels, ncol = 1,
                                 align = "v", axis = "lr",
                                 rel_heights = heights)

  # Sized so that every cell type gets a readable line and the strips below keep
  # their height: too wide and the dot plot flattens, too short and the labels
  # collide.
  if (is.null(width)) width <- max(9, 0.5 * length(clusters) + 7)
  if (is.null(height)) {
    height <- max(6.5, 0.30 * length(celltypes) + 2.2 +
                    0.5 * (length(ent_plots) + as.integer(!is.null(p_n))))
  }
  if (!is.null(save)) {
    ggplot2::ggsave(save, combined, width = width, height = height, dpi = dpi,
                    limitsize = FALSE, bg = "white")
  }
  list(plot = combined, main = p_main, n_cells = p_n, entropy = ent_plots,
       data = df, width = width, height = height)
}
