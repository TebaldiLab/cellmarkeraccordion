.onLoad <- function(libname, pkgname) {
  library(reticulate)

  # Check if Python is available; if not, configure Miniconda automatically
  if (!py_available(initialize = FALSE)) {
    # Automatically install Miniconda if needed
    if (!miniconda_exists()) {
      install_miniconda()
    }
    use_miniconda(required = TRUE)
  }
}
