#' Internal helpers for Suerat v4 (slot) / v5 (layer) compatibility
#'
#' These helpers centralize that logic so the rest of the package never has to 
#' branch on the version again. 
#' @keywords internal 
#' @noRd
NULL

.use_layers <- function(){
    utils::packageVersion("SeuratObject") >= "5.0.0"
}

#' Get an assay layer/slot independently of the Seurat version 
#' @keywords internal 
#' @noRd

.accordion_get_layer <- function(object, assay, layer){
    if(.use_layers()){
        tryCatch(
            SeuratObject::LayerData(object, assay = assay, layer = layer),
            error = function(e) NULL
        )
    } else{
        SeuratObject::SetAssayData(
            object, assay = assay, slot = layer, new.data = new.data
        )
    }
}

#' Join split layers (Seurat v5)
#' @keywords internal 
#' @noRd

.accordion_join_layers <- function(object, assay){
    if(!.use_layers()){
        return(object)
    }
    lyrs <- tryCatch(
        SeuratObject::Layers(object[[assay]]), 
        error = function(e) NULL
    )
    if(length(lyrs) > 1 &&
        (sum(grepl("^counts", lyrs)) > 1 || sum(grepl("^data", lyrs)) > 1)){
        object <- tryCatch(
            SeuratObject::JoinLayers(object, assay = assay), 
            error = function(e) object
        )
    }
    object
}

utils::globalVariables(c("new.data", ".accordion_set_layer", "accordion_join_layers"))
