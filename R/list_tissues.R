#' List tissues available in the Cell Marker Accordion database
#' @docType function
#' @name list_tissues
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available tissues.
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#' @param celltype Character string or character string vector specifying the
#'   celltype for which to extract the associate list of available tissues.
#'    If NULL, information from all cell types are retrieved.
#' @return List of tissues available in the Cell Marker Accordion database
#' @import data.table
#' @export
list_tissues<-function(species = c("Human","Mouse"),
                       celltype = NULL
){
  data(accordion_marker, package = "cellmarkeraccordion")
  input_species<-species
  if(!is.null(celltype)){
    output_table<-accordion_marker[species %in% input_species & CL_celltype %in% celltype]
  } else{
    output_table<-accordion_marker[species %in% input_species]
  }
  return(unique(output_table$Uberon_tissue))
}

