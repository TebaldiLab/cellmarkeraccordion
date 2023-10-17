#' Automatically annotating single-cell populations with custom marker genes
#' sets
#'
#' This function performs cell types or signatures/pathways annotation based on
#' cusom marker genes set. It takes in input either a Seurat object or a raw or
#' normalized count matrix and a table of marker genes associated to cell types
#' or even to pathways and return in output the cell types/pathways assignment
#' (added to the Seurat object or as a list).
#'
#'@param data Either a  Seurat object (version 4.9) or a raw or normalized count
#'  matrix with genes on rows and cells on columns. If raw counts are provided,
#'  data are log-normalized exploiting the NormalizeData() function from the
#'  Seurat package.
#'@param marker_table Data table or data frame containing cell type
#'  markers. The table needs to have at least two columns, the
#'  \code{category_column},  which specifies cell types or categories, and the
#'  \code{marker_column}, which specifies the corresponding markers on each row.
#'  Columns indicating the marker type (either positive or negative), and the
#'  marker weight can be optionally included.
#'@param category_column String characters specifying the name of the
#'  \code{marker_table} column containing cell types or categories.
#'  Default is “cell_type”.
#'@param marker_column String characters specifying the name of the
#'  \code{marker_table} column containing markers. Default is “marker”.
#'@param marker_type_column Optional string characters specifying the name of
#'  the \code{marker_table} column containing string characters
#'  indicating the type of markers, either “positive” or “negative”. If no
#'  \code{marker_type_column} is found in the \code{marker_table} all
#'  markers are considered “positive”. Default is “marker_type”.
#'@param weight_column  Optional string characters specifying the name of the
#'  \code{marker_table} column containing numeric value indicating the
#'  weight for each marker. If no \code{weight_column} is found in the
#'  \code{marker_table} all markers are equally weighted as 1. Default
#'  is “weight”.
#'@param cluster_info in case \code{object} is a Seurat object,
#'  \code{cluster_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains cluster ids; if
#'  \code{object} is a count matrix,
#'@param assay Character string specifying the Assay of the Seurat object. This
#'  parameter is necessary only  in case \code{data} is a Seurat object. Default
#'  is “RNA”.
#'@param min_n_marker Integer value specifying the minimum number of markers to
#'  keep for each cell type. Only cell types with a number of markers >= this
#'  threshold are kept.  Default is 5.
#'@param max_n_marker Integer value specifying the maximum number of markers to
#'  keep for each cell type. For the selection, markers are ranked according to
#'  their combined score, obtained by multiplying evidence consistency score and
#'  specificity score. If  NULL, no filter is applied. Default is NULL.
#'@param annotation_resolution Character string or character string vector
#'  specifying the resolution of the annotation. Either “cluster” and/or “cell”
#'  are supported. Default is “cluster”.
#'@param cluster_score_quantile_threshold numeric value in [0,1] specifying the
#'  cluster score quantile threshold. For each cell a score specific for each
#'  cell type is computed. To annotate a cluster cl, for each cell type the
#'  \code{cluster_score_quantile_threshold} is computed across cells belonging
#'  to that cluster and the cell type with the maximum score is then assigned to
#'  the cluster cl. Default is 0.75.
#'@param allow_unknown Logical value indicating whether to allow cells or
#'  clusters to be labeled as “unknown”. If it is set to TRUE, cells or clusters
#'  with negative scores are assigned to the “unknown” category. Default is
#'  TRUE.
#'@param annotation_name Character string specifying the name of the column in
#'  either the metadata of the input Seurat object or in the input
#'  \code{cluster_info} where the annotation will be stored. Per cluster and per
#'  cell annotation results will be stored in the
#'  \code{annotation_name}_per_cluster and \code{annotation_name}_per_cell
#'  columns respectively.
#'  If \code{include_detailed_annotation_info} parameter is set to TRUE, the
#'  detailed information the stored in a list named \code{annotation_name}.
#'  Default is “accordion_custom”.
#'@param include_detailed_annotation_info Logical value indicating whether to
#'  store information on the top cell types and markers in the output. If TRUE,
#'  a nested list named \code{annotation_name} is created. If
#'  \code{resolution_annotation} is set to “cluster” and/or “cell, sublists
#'  named “cluster_resolution” and/or “cell_resolution” are then added. Inside
#'  the sublist “detailed_annotation_info” the \code{n_top_markers} markers,
#'  group by \code{group_markers_by} and the \code{n_top_celltypes} cell types
#'  are then included. If a Seurat object is provided as input the list is
#'  stored in the misc slot of the object (object@misc@\code{annotation_name}). If the input
#'  is a count matrix, the list is returned in the final output. Default is
#'  TRUE.
#'   @param condition_group_info in case \code{object} is a Seurat object,
#'  \code{condition_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains condition ids for each cell;
#'  if \code{object} is a count matrix, \code{condition_group_info} should be need to be a
#'   data frame or data table containing condition identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “condition”, which specifies
#'   the condition id’s for each cell.  Default is NULL.
#'  @param cell_type_group_info in case \code{object} is a Seurat object,
#'  \code{cell_type_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains cell types ids for each cell;
#'  if \code{object} is a count matrix, \code{cell_type_group_info} should be need to be a
#'   data frame or data table containing cell types identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “cell_type”, which specifies
#'   the cell types for each cell.  Default is NULL.
#'@param group_markers_by Character string or character string vector specifying
#'  the classification of marker genes. It possible to retrieve
#'  \code{n_top_markers} marker genes for each cell type identified with cluster
#'  ("celltype_cluster") or cell ("celltype_cell") resolution;
#'  \code{n_top_markers} marker genes per cluster ("cluster") or per cell
#'  ("cell") can be also obtained. Either "celltype_cluster", "celltype_cell",
#'  "cluster" and/or "cell". Default is "celltype_cluster".
#'@param n_top_celltypes Integer value specifying the number of the top cell
#'  types to be included in the output for each cluster and cell depending on
#'  the selected \code{annotation_resolution} parameter Default is 5.
#'@param n_top_markers Integer value specifying the number of the top markers to
#'  be included in the output for each cell type, cluster or cell depending on
#'  the selected \code{annotation_resolution} and \code{group_markers_by}
#'  parameters. Default is 5.
#'@param top_marker_score_quantile_threshold numeric value in (0,1] specifying
#'  the marker score quantile threshold. For each marker a score specific for
#'  each cell is computed. To identify the \code{n_top_markers} for a cluster cl
#'  or a cell type ct, the \code{top_marker_score_quantile_threshold} is
#'  computed across cells belonging to that cluster or labeled as ct, and the
#'  \code{n_top_markers} with the maximum score are reported. Default is 0.75.

#' @return A Seurat object or a list
#' @details If a Seurat object was provided in input, the function returns the
#' Seurat object with markers-based scaled data in the scale.data slot and cell
#' types annotation results in the metadata. If
#' \code{include_detailed_annotation_info} and \code{plot} were set to TRUE, a
#' list containing cell types and markers information, together with ggplot
#' objects, is stored in the “misc@\code{annotation_name}” slot. If a count matrix was
#' provided in input, the function returns a list containing the following
#' elements:
#'
#' \itemize{
#' \item{"scaled_matrix":}{normalized and scaled expression matrix;}
#' }
#' If \code{annotation_resolution} is set to “cell”:
#' \itemize{
#' \item{"cell_annotation":}{data table containing cell types annotation results for each cell;}
#' }
#' If \code{annotation_resolution} is set to “cluster”:
#' \itemize{
#' \item{"cluster_annotation":}{data table containing cell types annotation results for each cell;}
#' }
#' If \code{include_detailed_annotation_info} is set to TRUE:
#' \itemize{
#' \item{"\code{annotation_name}":}{list containing detailed information of cell types annotation.}
#' }
#' @import scales
#' @import plyr
#' @import data.table
#' @import Seurat
#' @import ggplot2
#' @import stringr
#' @export
accordion_custom_annotation<-function(data,
                           marker_table,
                           category_column = "cell_type",
                           marker_column = "marker",
                           marker_type_column = "marker_type",
                           weight_column = "weight",
                           cluster_info = "seurat_clusters",
                           assay = "RNA",
                           min_n_marker = 5,
                           max_n_marker = NULL,
                           annotation_resolution = "cluster",
                           cluster_score_quantile_threshold = 0.75,
                           allow_unknown = TRUE,
                           annotation_name = "accordion_custom",
                           include_detailed_annotation_info = TRUE,
                           condition_group_info = NULL,
                           cell_type_group_info = NULL,
                           group_markers_by = "celltype_cluster",
                           n_top_celltypes = 5,
                           n_top_markers = 5,
                           top_marker_score_quantile_threshold = 0.75
){
  #count matrix  data
  #check the type of input (Seurat object or raw count matrix)
  if(class(data) != "Seurat"){

    #check that is not an empty count matrix
    if(sum(dim(data)) == 0){
      stop("Count matrix is empty")
    }
    data_type<-"matrix"
    #check if the first column is the gene columns
    if(class(data[,1]) == "character"){
      setDF(data)
      # Set the barcodes as the row names
      rownames(data) <- data[[1]]
      data[[1]] <- NULL
    }
    if(class(data) != "dgCMatrix"){
      data <- as(as.matrix(data), "sparseMatrix")

    }
    data <- CreateAssaydata(counts = data)
    #check that cluster_info is present if cluster is in annotation_resolution
    #if both cluster and cell resolution are set if the cluster_info is not provided or is not correct, only the per cell annotation is performed
    if("cluster" %in% annotation_resolution & "cell" %in% annotation_resolution){
      if(is.null(cluster_info)){
        warning("cluster_info not found. Please provide a data table or data frame specifying cell clusters to perform per cluster annotation.Cell types annotation will be perform only with per cell resolution.")
      } else {
        if(class(cluster_info) != "data.table" | class(cluster_info) != "data.frame"){
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
        if(class(cluster_info) != "data.table" | class(cluster_info) != "data.frame"){
          stop("Invalid input type. cluster_info needs to be a data table or data frame specifying cell clusters to perform per cluster annotation, or set cell in the annotation_resolution parameter to perform annotation with per cell resolution.")
        } else{ #if exists check that contain columns name
          if(!("cell" %in% colnames(cluster_info))){
            stop("cell column not found in cluster_info. Please provide a data table or data frame with a column named cell contaning cell ids.")
          }
          if(!("cluster" %in% colnames(cluster_info))){
            stop("cluster column not found in cluster_info. Please provide a data table or data frame with a column named cluster contaning cluster ids.")
          } else if("cell" %in% colnames(cluster_info) & "cluster" %in% colnames(cluster_info)){
            cluster_table<-as.data.table(cluster_info)[,c("cell","cluster")]
            colnames(cluster_table)<-c("cell","seurat_clusters")
          }
        }
      }
    }

    #Seurat object

  } else{
    data_type<-"seurat"

    #check assay
    if(is.null(data@assays[[assay]])){
      stop("Invalid assay provided")
    } else{

      DefaultAssay(data)<-assay
      #check that the Seurat object not contain an empty count matrix
      if (assay != "integrated"){
        if(sum(dim(data@assays[[assay]]@counts))==0){
          stop("Count matrix is empty")
        }
      }
      #check that the cluster column is present in the data
      if("cluster" %in% annotation_resolution & "cell" %in% annotation_resolution){
        if(class(cluster_info) != "character"){
          warning("Invalid input type: cluster_info needs to be a character string specifying the name of the column in the meta data containing cluster id's. Cell types annotation will be perform only with per cell resolution.")
        } else if (!cluster_info %in% colnames(data@meta.data)){
          warning(paste0(eval(cluster_info), " meta data column not found. Please provide a valid character string specifying the name of the column in the meta data containing cluster id's. Cell types annotation will be perform only with per cell resolution."))
        } else if (cluster_info %in% colnames(data@meta.data)){
          seurat_clusters<-cluster_info
        }
      } else if ("cluster" %in% annotation_resolution & !("cell" %in% annotation_resolution)){
        if(class(cluster_info) != "character"){
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

  #check input tabel
  if(nrow(marker_table) < 2){
    stop("Insufficient number of columns. The marker table must contains at least \"cell_type\"  and \"marker\" columns ")
  } else {
    if(category_column %in% colnames(marker_table) & marker_column %in% colnames(marker_table)){
      if(!marker_type_column %in% colnames(marker_table)){
        marker_table[,marker_type:="positive"]
      }
      if(!weight_column %in% colnames(marker_table)){
        marker_table[,weight:=1]
      }
      col_vec<-c(category_column, marker_column, marker_type_column, weight_column)
      marker_table<-marker_table[,..col_vec]
      colnames(marker_table)<-c("cell_type","marker","marker_type","weight")
    }

  }
  # keep only the max_n_marker genes for each cell type
  if(!is.null(max_n_marker)){
    if(!is.numeric(max_n_marker) | !(max_n_marker %in% 1 == 0)){
      warning("Invalid max_n_marker type. Parameter max_n_marker must be an integer value. No filter is applied")
    } else {
      accordion_marker<-accordion_marker[order(-combined_score)][,head(.SD, max_n_marker), by="cell_type"]
    }
  }

  # number of markers for each cell type
  accordion_marker[,length:= .N, by="cell_type"]
  if(!is.null(min_n_marker)){
    if(!is.numeric(min_n_marker) | !(min_n_marker %in% 1 == 0)){
      warning("Invalid min_n_marker type. Parameter min_n_marker must be an integer value. No filter is applied")
    } else{
      accordion_marker<-accordion_marker[length >= min_n_marker]
    }
  }
  #check input group_markers_by
  if("cell" %in% annotation_resolution & !("cluster" %in% annotation_resolution)){
    if(!(group_markers_by %in% c("celltype_cell","cell"))){
      group_markers_by<-"celltype_cell"
    }
  }

  if(sum(dim(data@assays[[assay]]@counts))!=0){
    #perform data normalization if not already performed
    if(identical(data@assays[[assay]]@counts, data@assays[[assay]]@data)){
      data <- NormalizeData(data)
    }
  }
  # subselect genes only found in data
  marker_table<-marker_table[marker %in% rownames(data)]

  #Evidence consistency score log-transformed
  marker_table[,weight_scaled := log10(weight)+1]

  #compute specificity for positive and negative markers
  mark_spec<-ddply(marker_table,.(marker,marker_type),nrow)
  colnames(mark_spec)<-c("marker","marker_type","specificity")
  marker_table<-merge(marker_table,mark_spec,by=c("marker","marker_type"),all.x = TRUE)

  length_ct_pos<-uniqueN(marker_table[marker_type=="positive"]$cell_type)
  length_ct_neg<-uniqueN(marker_table[marker_type=="negative"]$cell_type)

  #scale and log transforme specificity
  marker_table<-marker_table[marker_type=="positive",specificity_scaled := scales::rescale(as.numeric(specificity), to = c(1,length_ct_pos),from = c(length_ct_pos,1))
  ][marker_type=="negative",specificity_scaled := scales::rescale(as.numeric(specificity), to = c(1,length_ct_neg),from = c(length_ct_neg,1))
  ][,c("cell_type","marker","marker_type","specificity","specificity_scaled","weight_scaled","weight")]
  marker_table[,specificity_scaled:=log10(specificity_scaled)+1]

  setkey(marker_table,marker,cell_type)

  # merge Z_scaled_dt and accordion table
  marker_table[,combined_score := specificity_scaled * weight_scaled]


  # scale data based on markers used for the annotation
  data<-ScaleData(data, features = unique(marker_table$marker))
  scale_data_mat<-data@assays[[assay]]@scale.data

  Zscaled_data<-setDT(as.data.frame(scale_data_mat))
  Zscaled_data[,marker:=rownames(scale_data_mat)]

  setkey(Zscaled_data, marker)
  Zscaled_m_data<-melt.data.table(Zscaled_data,id.vars = c("marker"))
  colnames(Zscaled_m_data)<-c("marker","cell","expr_scaled")

  # compute the score for each cell
  dt_score<-merge.data.table(Zscaled_m_data,marker_table, by="marker",allow.cartesian = TRUE)
  dt_score[,score := expr_scaled * combined_score]
  dt_score_ct <- unique(dt_score[, c("cell_type", "cell")])
  setkey(dt_score, cell_type, cell, marker_type)
  sum_dt <- dt_score[data.table("cell_type" = rep(dt_score_ct$cell_type, each = 2),
                                "cell" = rep(dt_score_ct$cell, each = 2),
                                "marker_type" = c("positive", "negative")),
                     .(score= (sum(score)/(sqrt((sum(weight_scaled * specificity_scaled)))))), by = .EACHI]


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
    anno_dt_cl<-final_dt_cluster[order(-quantile_score_cluster)][,head(.SD, 1),"seurat_clusters"][,-c("cell","diff_score")]
    anno_dt_cl<-anno_dt_cl[,c("seurat_clusters","annotation_per_cell","quantile_score_cluster")]
    colnames(anno_dt_cl)<-c("seurat_clusters","annotation_per_cluster","quantile_score_cluster")

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
      anno_dt_cell[diff_score < 0, annotation_per_cell:= "unknown"][diff_score > 0, annotation_per_cell := cell_type]

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
        accordion_output<-list(data@assays[[assay]]@scale.data, cell_table)
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
                                                      cell_type_group_info)
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
                                                                cell_type_group_info)
    }

  }
    if(data_type == "seurat"){
      return(data)
    } else{
      return(accordion_output)
    }
  }
