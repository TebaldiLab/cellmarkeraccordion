#' List cell types available in the Cell Marker Accordion database
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available cell types.
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#' @param tissue Character string or character string vector specifying the tissue
#'   for which to extract the associate list of available cell types.
#'    If NULL, information from all tissues are retrieved.
#' @return List of cell types available in the Cell Marker Accordion database
#' @import data.table
#' @export
list_celltypes<-function(species = c("Human","Mouse"),
               tissue = NULL
){
  data("accordion_marker", package = "cellmarkeraccordion",envir = environment())

  input_species<-species
  if(!is.null(tissue)){
    output_table<-accordion_marker[species %in% input_species & Uberon_tissue %in% tissue]
  } else{
    output_table<-accordion_marker[species %in% input_species]
  }
  return(unique(output_table$CL_celltype))
}

