#' List tissues available in the Cell Marker Accordion disease database
#'
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available tissues.
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#'  @param disease Character string or character string vector specifying
#'   diseases to consider. If NULL, information from all diseases are considered.
#'   Default is NULL.
#'  @param celltype Character string or character string vector specifying the
#'   celltype for which to extract the associate list of available tissues.
#'    If NULL, information from all cell types are retrieved.
#' @return List of tissues available in the Cell Marker Accordion disease
#' database

list_disease_tissues<-function(species = c("Human","Mouse"),
                               disease = NULL,
                       aberrant_celltype = NULL
){
  load("accordion_disease_Ncit.rda")
  input_species<-species
  if(!is.null(celltype)){
    output_table<-accordion_disease[species %in% input_species &  NCIT_celltype %in% aberrant_celltype]
  }
  if(!is.null(disease)){
    output_table<-accordion_disease[species %in% input_species &  DO_diseasetype %in% disease]
  }
  if(is.null(celltype) & is.null(disease)){
    output_table<-accordion_disease[species %in% input_species]
  }
  return(unique(output_table$Uberon_tissue))
}

