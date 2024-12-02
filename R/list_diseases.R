#' List diseases available in the Cell Marker Accordion disease database
#'
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available diseases
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#'  @param tissue Character string or character string vector specifying the tissue
#'   for which to extract the associate list of available diseases.
#'    If NULL, information from all tissues are retrieved.
#'  @param aberrant_celltype Character string or character string vector
#'  specifying the aberrant celltype for which to extract the associate list of
#'  available tissues. If NULL, information from all aberrant cell types are
#'  retrieved.
#'  @return List of diseases available in the Cell Marker Accordion disease
#'  database
#' @import data.table
#' @export


list_diseases<-function(species = c("Human","Mouse"),
                       tissue = NULL,
                       aberrant_celltype = NULL
){
  data(accordion_disease)
  input_species<-species
  if(!is.null(tissue)){
    output_table<-accordion_disease[species %in% input_species & Uberon_tissue %in% tissue]
  }
  if(!is.null(aberrant_celltype)){
    output_table<-accordion_disease[species %in% input_species & NCIT_celltype %in% aberrant_celltype]
  }
  if( is.null(tissue) & is.null(aberrant_celltype)){
    output_table<-accordion_disease[species %in% input_species]
  }

  return(unique(output_table$DO_diseasetype))
}

