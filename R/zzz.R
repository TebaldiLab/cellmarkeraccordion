owlready2 <- NULL
.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  owlready2 <<- reticulate::import("owlready2", delay_load = TRUE)
}
