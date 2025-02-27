#' Automatically identify and interpreting cell cycle state of single-cell
#' populations
#'
#' This function identifies cell cycle states exploiting the collection of
#' marker genes associated to each phase, including G0. It takes in input either
#' a Seurat object or a raw or normalized count matrix and return in output the
#' cell cycle assignment and the detailed informations of the annotation results
#' (added to the Seurat object or as a list).
#'
#' @param data Either a  Seurat object (version 4 or 5) or a raw or normalized
#'   count matrix with genes on rows and cells on columns. If raw counts are
#'   provided, data are log-normalized exploiting the NormalizeData() function
#'   from the Seurat package.
#' @param cluster_info in case \code{data} is a Seurat object,
#'   \code{cluster_info} should be need to be a character string specifying the
#'   name of the column in the metadata that contains cluster ids; if
#'   \code{data} is a count matrix, \code{cluster_info} should be need to be a
#'   data frame or data table containing cluster identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “cluster”, which specifies
#'   the clustering id’s for each cell. This parameter is necessary only when
#'   the input is a count matrix and only if the \code{annotation_resolution}
#'   parameter is set to “cluster”. Default is “seurat_clusters”.
#' @param assay Character string specifying the Assay of the Seurat object. This
#'   parameter is necessary only  in case \code{data} is a Seurat object.
#'   Default is “RNA”.
#' @param species Character string or character string vector specifying the
#'   species. Currently, either “Human” and/or “Mouse” are supported. If
#'   multiple species are selected, marker genes are merged together. Default is
#'   “Human”.
#' @param annotation_resolution Character string or character string vector
#'   specifying the resolution of the annotation. Either “cluster” and/or “cell”
#'   are supported. Default is “cell”.
#' @param annotation_name Character string specifying the name of the column
#'   in either the metadata of the input Seurat object or in the input
#'   \code{cluster_info} where the annotation will be stored. Per cluster and
#'   per cell annotation results will be stored in the
#'   \code{annotation_name}_per_cluster and \code{annotation_name}_per_cell
#'   columns respectively.
#'   If \code{include_detailed_annotation_info} parameter is set to TRUE, the
#'   detailed information the stored in a list named \code{annotation_name}.
#'   Default is “accordion_cell_cycle”.
#' @param cluster_score_quantile_threshold numeric value in (0,1) specifying the
#'   cluster score quantile threshold. For each cell a score specific for each
#'   cell type is computed. To annotate a cluster cl, for each cell type the
#'   \code{cluster_score_quantile_threshold} is computed across cells belonging
#'   to that cluster and the cell type with the maximum score is then assigned
#'   to the cluster cl. Default is 0.75.
#' @param allow_unknown Logical value indicating whether to allow cells or
#'   clusters to be labeled as “unknown”. If it is set to TRUE, cells or
#'   clusters with negative scores are assigned to the “unknown” category.
#'   Default is TRUE.
#' @param include_detailed_annotation_info Logical value indicating whether to
#'   store information on the top cell types and markers in the output. If TRUE,
#'   a nested list named \code{annotation_name} is created. If
#'   \code{resolution_annotation} is set to “cluster” and/or “cell, sublists
#'   named “cluster_resolution” and/or “cell_resolution” are then added. Inside
#'   the sublist “detailed_annotation_info” the \code{n_top_markers} markers,
#'   group by \code{group_markers_by} and the \code{n_top_celltypes} cell types
#'   are then included. If a Seurat object is provided as input the list is
#'   stored in the misc slot of the object (object@misc@\code{annotation_name}).
#'   If the input is a count matrix, the list is returned in the final output.
#'   Default is FALSE.
#' @param condition_group_info in case \code{data} is a Seurat object,
#'  \code{condition_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains condition ids for each cell;
#'  if \code{data} is a count matrix, \code{condition_group_info} should be need to be a
#'   data frame or data table containing condition identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “condition”, which specifies
#'   the condition id’s for each cell.  Default is NULL.
#' @param celltype_group_info in case \code{data} is a Seurat object,
#'  \code{celltype_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains cell types ids for each cell;
#'  if \code{data} is a count matrix, \code{celltype_group_info} should be need to be a
#'   data frame or data table containing cell types identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “cell_type”, which specifies
#'   the cell types for each cell.  Default is NULL.
#' @param group_markers_by Character string or character string vector
#'   specifying the classification of marker genes. It possible to retrieve
#'   \code{n_top_markers} marker genes for each cell type identified with
#'   cluster ("celltype_cluster") or cell ("celltype_cell") resolution;
#'   \code{n_top_markers} marker genes per cluster ("cluster") or per cell
#'   ("cell") can be also obtained. Either "celltype_cluster", "celltype_cell",
#'   "cluster" and/or "cell". Default is "celltype_cell".
#' @param n_top_celltypes Integer value specifying the number of the top cell
#'   types to be included in the output for each cluster and cell depending on
#'   the selected \code{annotation_resolution} parameter Default is 5.
#' @param n_top_markers Integer value specifying the number of the top markers
#'   to be included in the output for each cell type, cluster or cell depending
#'   on the selected \code{annotation_resolution} and \code{group_markers_by}
#'   parameters. Default is 5.
#' @param top_marker_score_quantile_threshold numeric value in (0,1] specifying
#'   the marker score quantile threshold. For each marker a score specific for
#'   each cell is computed. To identify the \code{n_top_markers} for a cluster
#'   cl or a cell type ct, the \code{top_marker_score_quantile_threshold} is
#'   computed across cells belonging to that cluster or labeled as ct, and the
#'   \code{n_top_markers} with the maximum score are reported. Default is 0.75.
#' @param plot Logical value indicating whether to store plots displaying
#'   detailed annotation information.  This parameter can be set to TRUE only
#'   when \code{include_detailed_annotation_info} is set to TRUE. If TRUE,
#'   lollipop plots displaying the top \code{n_top_markers} group by
#'   \code{group_markers_by} and \code{n_top_celltypes} for each
#'   \code{annotation_resolution} together with the cell types hierarchies based
#'   on the cell ontology structure are stored in the \code{annotation_name}
#'   list. Default is FALSE.
#'
#' @return A Seurat object or a list
#' @details If a Seurat object was provided in input, the function returns the
#' Seurat object with markers-based scaled data in the scale.data slot and cell
#' types annotation results in the metadata. If
#' \code{include_detailed_annotation_info} and \code{plot} were set to TRUE, a
#' list named \code{annotation_name} containing cell types and markers
#' information, together with ggplot objects, is stored in the “misc” slot. If a
#' count matrix was provided in input, the function returns a list containing
#' the following elements:
#'
#' \describe{
#' \item{"scaled_matrix":}{normalized and scaled expression matrix;}
#' }
#' If \code{annotation_resolution} is set to “cell”:
#' \describe{
#' \item{"cell_annotation":}{data table containing cell types annotation results for each cell;}
#' }
#' If \code{annotation_resolution} is set to “cluster”:
#' \describe{
#' \item{"cluster_annotation":}{data table containing cell types annotation results for each cell;}
#' }
#' If \code{include_detailed_annotation_info} is set to TRUE:
#' \describe{
#' \item{"\code{annotation_name}":}{list containing detailed information of cell types annotation.}
#' }
#' @import scales
#' @import plyr
#' @import data.table
#' @import Seurat
#' @import ggplot2
#' @import stringr
#' @import knitr
#' @importFrom methods as
#' @importFrom stats quantile
#' @export
accordion_cellcycle<-function(data,
                              cluster_info = "seurat_clusters",
                              assay = "RNA",
                              species = "Human",
                              annotation_resolution = "cell",
                              annotation_name = "accordion_cell_cycle",
                              cluster_score_quantile_threshold = 0.75,
                              allow_unknown = FALSE,
                              include_detailed_annotation_info = FALSE,
                              condition_group_info = NULL,
                              celltype_group_info = NULL,
                              group_markers_by = "celltype_cell",
                              n_top_celltypes = 5,
                              n_top_markers = 5,
                              top_marker_score_quantile_threshold = 0.75,
                              plot = FALSE
                              ){


  #count matrix  data
  #check the type of input (Seurat data or raw count matrix)
  if(!(inherits(data, "Seurat"))){

    #check that is not an empty count matrix
    if(sum(dim(data)) == 0){
      stop("Count matrix is empty")
    }
    data_type<-"matrix"
    #check if the first column is the gene columns
    if(inherits(data[,1],"character")){
      setDF(data)
      # Set the barcodes as the row names
      rownames(data) <- data[[1]]
      data[[1]] <- NULL
    }
    if(inherits(data,"dgCMatrix")){
      data <- as(as.matrix(data), "sparseMatrix")

    }
    data <- CreateSeuratObject(counts = data)
    accordion_output<-list() #set output list
    #check that cluster_info is present if cluster is in annotation_resolution
    #if both cluster and cell resolution are set if the cluster_info is not provided or is not correct, only the per cell annotation is performed
    if("cluster" %in% annotation_resolution & "cell" %in% annotation_resolution){
      if(is.null(cluster_info)){
        warning("cluster_info not found. Please provide a data table or data frame specifying cell clusters to perform per cluster annotation.Cell types annotation will be perform only with per cell resolution.")
      } else {
        if(!(inherits(cluster_info, "data.table")) | !(inherits(cluster_info, "data.frame"))){
          warning("Invalid input type. cluster_info needs to be a data table or data frame specifying cell clusters to perform per cluster annotation. Cell types annotation will be perform only with per cell resolution.")
        } else { #if exists check that contain columns name
          if(!("cell" %in% colnames(cluster_info))){
            warning("cell column not found in cluster_info. Please provide a data table or data frame with a column named cell contaning cell ids. Cell types annotation will be perform only with per cell resolution.")
          }
          if(!("cluster" %in% colnames(cluster_info))){
            warning("cluster column not found in cluster_info. Please provide a data table or data frame with a column named cluster contaning cluster ids. Cell types annotation will be perform only with per cell resolution.")
          } else if("cell" %in% colnames(cluster_info) & "cluster" %in% colnames(cluster_info)){
            cluster_table<-as.data.table(cluster_info)[,c("cell","cluster")]
          }
        }
      }
      #if only per cluster resolution is set if the cluster_info is not provided or is not correct stop
    } else if ("cluster" %in% annotation_resolution & !("cell" %in% annotation_resolution)){
      if(is.null(cluster_info)){
        stop("cluster_info not found. Please provide a data table or data frame specifying cell clusters to perform per cluster annotation, or set cell in the annotation_resolution parameter to perform annotation with per cell resolution.")
      } else {
        if(!(inherits(cluster_info, "data.table")) | !(inherits(cluster_info, "data.frame"))){
          stop("Invalid input type. cluster_info needs to be a data table or data frame specifying cell clusters to perform per cluster annotation, or set cell in the annotation_resolution parameter to perform annotation with per cell resolution.")
        } else{ #if exists check that contain columns name
          if(!("cell" %in% colnames(cluster_info))){
            stop("cell column not found in cluster_info. Please provide a data table or data frame with a column named cell contaning cell ids.")
          }
          if(!("cluster" %in% colnames(cluster_info))){
            stop("cluster column not found in cluster_info. Please provide a data table or data frame with a column named \"cluster\" contaning cluster ids.")
          } else if("cell" %in% colnames(cluster_info) & "cluster" %in% colnames(cluster_info)){
            cluster_table<-as.data.table(cluster_info)[,c("cell","cluster")]
            colnames(cluster_table)<-c("cell","seurat_clusters")
          }
        }
      }
    }

    #Seurat data

  } else{
    data_type<-"seurat"

    #check assay
    if(is.null(data@assays[[assay]])){
      stop("Invalid assay provided")
    } else{

      DefaultAssay(data)<-assay
      #check that the Seurat data not contain an empty count matrix
      if (assay != "integrated"){
        if(sum(dim(GetAssayData(data, assay=assay, slot='counts')))==0){
          stop("Count matrix is empty")
        }
      }
      #check that the cluster column is present in the data
      if("cluster" %in% annotation_resolution & "cell" %in% annotation_resolution){
        if(!(inherits(cluster_info, "character"))){
          warning("Invalid input type: cluster_info needs to be a character string specifying the name of the column in the meta data containing cluster id's. Cell types annotation will be perform only with per cell resolution.")
        } else if (!cluster_info %in% colnames(data@meta.data)){
          warning(paste0(eval(cluster_info), " meta data column not found. Please provide a valid character string specifying the name of the column in the meta data containing cluster id's. Cell types annotation will be perform only with per cell resolution."))
        } else if (cluster_info %in% colnames(data@meta.data)){
          seurat_clusters<-cluster_info
        }
      } else if ("cluster" %in% annotation_resolution & !("cell" %in% annotation_resolution)){
        if(!(inherits(cluster_info, "character"))){
          stop("Invalid input type: cluster_info needs to be a character string specifying the name of the column in the meta data containing cluster id's.")
        } else if (!cluster_info %in% colnames(data@meta.data)){
          stop(paste0(eval(cluster_info), " meta data column not found. Please provide a valid character string specifying the name of the column in the meta data containing cluster id's."))
        } else if(cluster_info %in% colnames(data@meta.data)){
          cluster_table<-as.data.table(data@meta.data)[,cell:=rownames(data@meta.data)]
          col<-c("cell",eval(cluster_info))
          cluster_table<-cluster_table[, ..col]
          colnames(cluster_table)<-c("cell","seurat_clusters")
        }
      }
    }
  }

  data("cell_cycle_markers", package = "cellmarkeraccordion",envir = environment())

  if((length(species) ==1 & !(species %in% c("Human","Mouse"))) | (length(species) == 2 & setequal(species, c("Human","Mouse")))){
    warning("Invalid species type")
    if(all(rownames(data)[1:10] %in% toupper(rownames(data)[1:10]))){
      cell_cycle_markers<-cell_cycle_markers[species %in% "Human"]
      warning("The dataset might be human. Human markers are indeed used.")
    } else{
      cell_cycle_markers<-cell_cycle_markers[species %in% "Mouse"]
      warning("The dataset might be mouse. Mouse markers are indeed used.")
    }
  }

  # check group_markers_by input
  if(!(group_markers_by %in% c("cluster","celltype_cluster","cell","celltype_cell","score_cell"))){
    warning("invalid group_by. Please select \"cluster\",\"celltype_cluster\", \"cell\", \"celltype_cell\" or \"score_cell\"")
    if("cluster" %in% annotation_resolution){
      group_markers_by<-"celltype_cluster"
    } else if("cell" %in% annotation_resolution){
      group_markers_by<-"celltype_cell"
    }
  }
  if("cluster" %in% annotation_resolution){
    if(!(group_markers_by %in% c("cluster","celltype_cluster"))){
      group_markers_by<-"celltype_cluster"
    }
  }
  if("cell" %in% annotation_resolution){
    if(!(group_markers_by %in% c("cell","celltype_cell","score_cell"))){
      group_markers_by<-"celltype_cell"
    }
  }

  #avoid warnings
  suppressWarnings({
    if(sum(dim(GetAssayData(data, assay=assay, slot='counts')))!=0){
      #perform data normalization if not already performed
      if(identical(GetAssayData(data, assay=assay, slot='counts'), GetAssayData(data, assay=assay, slot='data')) | sum(dim(GetAssayData(data, assay=assay, slot='data')))==0){
        data <- NormalizeData(data)
      }
    }
  })

    # subselect genes only found in data
    cell_cycle_markers<-cell_cycle_markers[marker %in% rownames(data)]
    cell_cycle_markers[,weight:=1]
    #Evidence consistency score log-transformed
    cell_cycle_markers[,weight_scaled := log10(weight)+1]

    #compute SPs for positive and negative markers
    mark_spec<-ddply(cell_cycle_markers,.(marker,marker_type),nrow)
    colnames(mark_spec)<-c("marker","marker_type","SPs")
    cell_cycle_markers<-merge(cell_cycle_markers,mark_spec,by=c("marker","marker_type"),all.x = TRUE)

    length_ct_pos<-uniqueN(cell_cycle_markers[marker_type=="positive"]$cell_type)
    length_ct_neg<-uniqueN(cell_cycle_markers[marker_type=="negative"]$cell_type)

    #scale and log transforme SPs
    cell_cycle_markers<-cell_cycle_markers[marker_type=="positive",SPs_reg := scales::rescale(as.numeric(SPs), to = c(1,length_ct_pos),from = c(length_ct_pos,1))
    ][marker_type=="negative",SPs_reg := scales::rescale(as.numeric(SPs), to = c(1,length_ct_neg),from = c(length_ct_neg,1))
    ][,c("cell_type","marker","marker_type","SPs","SPs_reg","weight_scaled","weight")]
    cell_cycle_markers[,SPs_reg:=log10(SPs_reg)+1]

    setkey(cell_cycle_markers,marker,cell_type)

    # merge Z_scaled_dt and accordion table
    cell_cycle_markers[,combined_score := SPs_reg * weight_scaled]

    # store original scale.data slot if present
    if(sum(dim(GetAssayData(data, assay=assay, slot='scale.data')))!=0){
      orig.scale_data<-GetAssayData(data, assay=assay, slot='scale.data')
    }

    # scale data based on markers used for the annotation
    data<-ScaleData(data, features = unique(cell_cycle_markers$marker))
    Zscaled_data<-GetAssayData(data, assay=assay, slot='scale.data')
    Zscaled_data<-as.data.table(as.data.frame(Zscaled_data),keep.rownames = "marker")
    setkey(Zscaled_data, marker)
    Zscaled_m_data<-melt.data.table(Zscaled_data,id.vars = c("marker"))
    colnames(Zscaled_m_data)<-c("marker","cell","expr_scaled")

    # compute the score for each cell
    dt_score<-merge.data.table(Zscaled_m_data,cell_cycle_markers, by="marker",allow.cartesian = TRUE)
    dt_score[,score := expr_scaled * combined_score]
    dt_score_ct <- unique(dt_score[, c("cell_type", "cell")])
    setkey(dt_score, cell_type, cell, marker_type)
    sum_dt <- dt_score[data.table("cell_type" = rep(dt_score_ct$cell_type, each = 2),
                                  "cell" = rep(dt_score_ct$cell, each = 2),
                                  "marker_type" = c("positive", "negative")),
                       .(score= (sum(score)/(sqrt((sum(weight_scaled * SPs_reg)))))), by = .EACHI]


    sum_dt<-unique(sum_dt)
    sum_dt[is.na(score), score := 0]
    final_dt <- sum_dt[marker_type == "positive"
    ][, diff_score := score - sum_dt[marker_type == "negative", score]
    ][, marker_type := NULL][,score := NULL]


    # annotation per cluster
    if ("cluster" %in% annotation_resolution){
      setkey(cluster_table, cell)
      final_dt_cluster<-merge.data.table(final_dt, cluster_table, by="cell")
      final_dt_cluster[, quantile_score_cluster:= quantile(diff_score,probs = cluster_score_quantile_threshold, na.rm=TRUE), by=c("seurat_clusters","cell_type")]
      if (allow_unknown == T){
        final_dt_cluster[quantile_score_cluster < 0, annotation_per_cell:= "unknown"][quantile_score_cluster > 0, annotation_per_cell := cell_type]
      } else {
        final_dt_cluster[, annotation_per_cell := cell_type]
      }
      # add the annotation results in the metadata of the Seurat data
      anno_dt_cell<-final_dt_cluster[order(-diff_score)][,head(.SD, 1),"cell"]
      anno_dt_cell_ptc<-anno_dt_cell[,ncell_celltype_cluster:= .N,by=c("seurat_clusters","annotation_per_cell")]
      anno_dt_cell_ptc[,ncell_tot_cluster:= .N, by="seurat_clusters"]
      anno_dt_cell_ptc[,perc_celltype_cluster:= round((ncell_celltype_cluster/ncell_tot_cluster)*100, digits = 2)]
      anno_dt_cl<-anno_dt_cell_ptc[order(-quantile_score_cluster)][,head(.SD, 1),"seurat_clusters"][,-c("cell","cell_type","diff_score","ncell_tot_cluster","ncell_celltype_cluster")]
      anno_dt_cl<-anno_dt_cl[,c("seurat_clusters","annotation_per_cell","quantile_score_cluster","perc_celltype_cluster")]
      colnames(anno_dt_cl)<-c("seurat_clusters","annotation_per_cluster","quantile_score_cluster","percentage")

      #if less than 10% of cells are labeled as the top cell type assigned as unknown
      if (allow_unknown == T){
        anno_dt_cl[percentage < 10, annotation_per_cluster:= "unknown"]
      }

      name<-paste0(annotation_name,"_per_cluster")
      name_score<-paste0(annotation_name,"_per_cluster_score")

      if(data_type == "seurat"){
        data@meta.data[,name] = ""
        data@meta.data[,name_score] = ""

        for (cl in unique(anno_dt_cl$seurat_clusters)){
          data@meta.data[which(data@meta.data$seurat_clusters == cl),name]<- anno_dt_cl[seurat_clusters==cl]$annotation_per_cluster
          data@meta.data[which(data@meta.data$seurat_clusters == cl),name_score]<- anno_dt_cl[seurat_clusters==cl]$quantile_score_cluster

        }
      } else {
        cluster_table<-merge(cluster_table,anno_dt_cl[,c("seurat_clusters","annotation_per_cluster")], by="seurat_clusters")
        cluster_table<-cluster_table[,c("cell","seurat_clusters","annotation_per_cluster")]
        colnames(cluster_table)<-c("cell","cluster",eval(name))

        accordion_output<-list(data@assays[[assay]]@scale.data, cluster_table)
        names(accordion_output)<-c("scaled_matrix","cluster_annotation")
      }

    }

    # annotation per cell
    if ("cell" %in% annotation_resolution){
      name<-paste0(annotation_name,"_per_cell")
      name_score<-paste0(annotation_name,"_per_cell_score")

      data@meta.data[,name] = ""
      data@meta.data[,name_score] = ""
      anno_dt_cell<-final_dt[order(-diff_score)][,head(.SD, 1),"cell"]
      if (allow_unknown == T){
        anno_dt_cell[diff_score < 0, annotation_per_cell:= "allow"][diff_score > 0, annotation_per_cell := cell_type]

      } else {
        anno_dt_cell[, annotation_per_cell := cell_type]

      }

      # add the annotation result to the data metadata
      if(!identical(colnames(data),anno_dt_cell$cell)){
        anno_dt_cell<-anno_dt_cell[order(match(anno_dt_cell$cell,colnames(data))),]
      }
      if(data_type == "seurat"){
        data@meta.data[,name]<-anno_dt_cell$annotation_per_cell
        data@meta.data[,name_score]<-anno_dt_cell$diff_score
      } else{
        cell_table<-anno_dt_cell[,c("cell","annotation_per_cell","diff_score")]
        colnames(cell_table)<-c("cell",eval(name), eval(name_score))

        if(!is_empty(accordion_output)){
          accordion_output<-append(accordion_output,cell_table)
          names(accordion_output)<-c(names(accordion_output), "cell_annotation")
        } else {
          accordion_output<-list(GetAssayData(data, assay=assay, slot='scale.data'), cell_table)
          names(accordion_output)<-c("scaled_matrix","cell_annotation")
        }

      }

    }

    if(include_detailed_annotation_info == T){
      if(data_type == "seurat"){
        data <- include_detailed_annotation_info_helper(data,
                                                        data_type,
                                                        annotation_resolution,
                                                        final_dt_cluster,
                                                        anno_dt_cl,
                                                        dt_score,
                                                        annotation_name,
                                                        group_markers_by,
                                                        dt_top_marker,
                                                        cluster_info,
                                                        final_dt,
                                                        anno_dt_cell,
                                                        n_top_celltypes,
                                                        n_top_markers,
                                                        top_marker_score_quantile_threshold,
                                                        condition_group_info,
                                                        celltype_group_info)
      } else{
        accordion_output<-include_detailed_annotation_info_helper(accordion_output,
                                                                  data_type,
                                                                  annotation_resolution,
                                                                  final_dt_cluster,
                                                                  anno_dt_cl,
                                                                  dt_score,
                                                                  annotation_name,
                                                                  group_markers_by,
                                                                  dt_top_marker,
                                                                  cluster_info,
                                                                  final_dt,
                                                                  anno_dt_cell,
                                                                  n_top_celltypes,
                                                                  n_top_markers,
                                                                  top_marker_score_quantile_threshold,
                                                                  condition_group_info,
                                                                  celltype_group_info)
      }

    }

    #re-assigned the original scale.data slot
    if(exists("orig.scale_data")){
      accordion_scale.data<-list()
      accordion_scale.data[["accordion_scale.data"]]<-GetAssayData(object = data, assay = assay, slot = "scale.data")
      data@misc[[annotation_name]]<-append(data@misc[[annotation_name]], accordion_scale.data)
      data[[assay]]$scale.data <- orig.scale_data
    }

    if(include_detailed_annotation_info==T & plot == T){
      if(data_type == "seurat"){
        data<-accordion_plot(data, info_to_plot = annotation_name, resolution = annotation_resolution, group_markers_by = group_markers_by)
      } else{
        accordion_output<-accordion_plot(accordion_output, info_to_plot = annotation_name, resolution = annotation_resolution, group_markers_by = group_markers_by)

      }
    }
    if(data_type == "seurat"){
      return(data)
    } else{
      return(accordion_output)
    }

  }

