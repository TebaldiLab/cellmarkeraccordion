#' List tissues available in the Cell Marker Accordion database
#'
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available tissues.
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#'  @param celltype Character string or character string vector specifying the
#'   celltype for which to extract the associate list of available tissues.
#'    If NULL, information from all cell types are retrieved.
#' @return List of tissues available in the Cell Marker Accordion database

list_tissues<-function(species = c("Human","Mouse"),
                       celltype = NULL
){
  load("accordion_marker_v2.0.rda")
  input_species<-species
  if(!is.null(celltype)){
    output_table<-accordion_marker_v2.0[species %in% input_species & CL_celltype %in% celltype]
  } else{
    output_table<-accordion_marker_v2.0[species %in% input_species]
  }
  return(unique(output_table$Uberon_tissue))
}

