.onLoad <- function(libname, pkgname) {
  library(reticulate)
  print("ciao")
  # Suppress reticulate prompt
  options(reticulate.miniconda.prompt = FALSE)

  # Check or configure Python environment
  if (!py_available(initialize = FALSE)) {
    if (!miniconda_exists()) {
      install_miniconda()
    }
    use_miniconda(required = TRUE)
  }

}
