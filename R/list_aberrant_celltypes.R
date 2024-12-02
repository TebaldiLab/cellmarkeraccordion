#' List aberrant cell types available in the Cell Marker Accordion disease database
#'
#' @param species Character string or character string vector specifying the
#'   species for which to extract the associate list of available cell types.
#'   Currently, either “Human” and/or “Mouse” are supported. Default is
#'   c("Mouse",“Human”).
#'  @param disease Character string or character string vector specifying
#'   diseases to consider. If NULL, information from all diseases are considered.
#'   Default is NULL.
#'  @param tissue Character string or character string vector specifying the tissue
#'   for which to extract the associate list of available cell types.
#'    If NULL, information from all tissues are retrieved.
#' @return List of aberrant cell types available in the Cell Marker Accordion
#' disease database
#' @import data.table
#' @export
list_aberrant_celltypes<-function(species = c("Human","Mouse"),
                                  disease=NULL,
                         tissue = NULL
){
  data(disease_accordion_marker)
  input_species<-species
  if(!is.null(tissue)){
    output_table<-disease_accordion_marker[species %in% input_species & Uberon_tissue %in% tissue]
  }
  if(!is.null(disease)){
    output_table<-disease_accordion_marker[species %in% input_species & DO_diseasetype %in% disease]
  }
  if(is.null(tissue) & is.null(disease)){
    output_table<-disease_accordion_marker[species %in% input_species]
  }
  return(unique(output_table$NCIT_celltype))
}
