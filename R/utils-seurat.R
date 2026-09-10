# Internal helpers for Seurat v4 (slot) / v5 (layer) compatibility

#' @keywords internal
#' @noRd

.use_layers <- function() {
    utils::packageVersion("SeuratObject") >= "5.0.0"
}

#' Get an assay layer/slot idenpendently of the Seurat version

#' @keywords internal
#' @noRd

.accordion_get_layer <- function(objectm assay, layer){
    if (.use_layers()){
        tryCatch(
            SeuratObject::LayerData(object, assay = assay, layer = layer), 
            error = function(e) NULL
        )
    } else{
        SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    }
}

#' Set an assay layer/slot independently of the Seurat version
#' @keywords internal
#' @noRd

.accordion_set_layer <- function(object, assay, layer, new.data){
    if (.use_layers()){
        SeuratObject::LayerData(object, assay = assay, layer = layer) <- new.data
        object
    } else {
        SeuratObject::SetAssayData(
            object, assay = assay, slot = layer, new.data = new.data
        )
    }
}

#' Join split layers (Seurat v5) so single-matrix accessors work
#' @keywords internal
#' @noRd

.accordion_join_layers <- function(object, assay){
    if (!.use_layers()){
        return(object)
    }
    lyrs <- tryCatch(
        SeuratObject::Layers(object[[assay]]), 
        error = function(e) NULL
    )
    if (length(lyrs) > 1 &&
        (sum(grepl("^counts", lyrs)) > 1 || sum(grepl("^data", lyrs)) > 1)){
            object <- tryCatch(
                SeuratObject::JoinLayers(object, assay = assay), 
                error = function(e) object
            )
        }
        object
}