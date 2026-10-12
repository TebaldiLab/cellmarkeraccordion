#' The object the Accordion plotting functions return
#'
#' `plot_top_celltypes_per_cluster()`, `plot_entropy()` and
#' `plot_spatial_entropy()` all return a list: `$plot` is the assembled figure
#' and the other elements are its pieces, so that a panel can be restyled or
#' reused on its own.
#'
#' A bare list, though, does not behave like a plot: printing it at the console
#' dumps the whole structure, and `ggplot2::ggsave()` stops with
#' *no applicable method for 'grid.draw' applied to an object of class "list"*.
#' Giving the list a class and three small methods fixes both, while `$plot`,
#' `$data` and the rest keep working exactly as before.
#'
#' @param x the object returned by one of the Accordion plotting functions.
#' @param ... passed on to the underlying print/plot method.
#' @param recording passed on to [grid::grid.draw()].
#' @return `print()` and `plot()` return `x` invisibly, after drawing it.
#' @name accordion_plot_object
NULL

#' @keywords internal
#' @noRd
.as_accordion_plot <- function(x) {
  class(x) <- c("accordion_plot", class(x))
  x
}

#' @rdname accordion_plot_object
#' @method print accordion_plot
#' @export
print.accordion_plot <- function(x, ...) {
  print(x$plot, ...)
  invisible(x)
}

#' @rdname accordion_plot_object
#' @method plot accordion_plot
#' @export
plot.accordion_plot <- function(x, ...) {
  print(x$plot, ...)
  invisible(x)
}

#' @rdname accordion_plot_object
#' @method grid.draw accordion_plot
#' @export
#' @importFrom grid grid.draw
grid.draw.accordion_plot <- function(x, recording = TRUE) {
  grid::grid.draw(x$plot, recording = recording)
}
