#' Detect the gene nomenclature used by the input dataset
#'
#' Internal helper. Marker genes of the Accordion database are stored with
#' human (upper case, e.g. "CD3E") or mouse (title case, e.g. "Cd3e") symbols. 
#' When markers of more than one species are merged together, they have to be 
#' converted to a single nomenclature, i.e. the one used by the input dataset. 

#' @param genes_names Character vector with the gene names of the input dataset.
#' @param markers Optional character vector with the marker symbols to be converted.
#' @return Either "Humna" or "Mouse".
#' @keywords internal
#' @noRd

detect_gene_nomenclature <- function(genes_names, markers = NULL){

    gene_names <- unique(as.character(gene_names))
    gene_names <- gene_names[!is.na(gene_names) & nzchar(gene_names)]

    if (length(gene_names) == 0){
        return("Human")
    }
    
    if (!is.null(markers) && length(markers) > 0){
        markers <- unique(as.character(markers))
        markers <- markers[!is.na(markers) & nzchar(markers)]
        n_human <- sum(unique(toupper(markers)) %in% gene_names)
        n_mouse <- sum(unique(str_to_title(markers)) %in% gene_names)
        if(n_human > 0 || n_mouse > 0){
            return(if(n_human >= n_mouse) "Human" else "Mouse")
        }
    }

    symbols <- gene_names[grepl("[[:alpha:]]", gene_names)]
    if(length(symbols) == 0){
        return("Human")
    }
    frac_upper <- mean(symbols == toupper(symbols))
    if(frac_upper >= 0.5) "Human" else "Mouse"
}

utils::globalVariables(c("invalid_species", ".input_species"))