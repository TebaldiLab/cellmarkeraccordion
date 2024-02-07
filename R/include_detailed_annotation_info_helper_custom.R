#' Helper function to include detailed annotation information
#' @noRd
#' @param data
#' @param data_type
#' @param annotation_resolution
#' @param final_dt_cluster
#' @param anno_dt_cl
#' @param dt_score
#' @param annotation_name
#' @param group_markers_by
#' @param dt_top_marker
#' @param cluster_info
#' @param final_dt
#' @param anno_dt_cell
#' @param n_top_celltypes
#' @param n_top_markers
#' @param top_marker_score_quantile_threshold
#' @param top_cell_score_quantile_threshold
#' @param condition_group_info
#' @param cell_type_group_info
#'
#' @return
#' @export
#'
#' @examples
include_detailed_annotation_info_helper_custom<-function(data,
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
                                                  top_cell_score_quantile_threshold,
                                                  condition_group_info,
                                                  cell_type_group_info){
  if(data_type == "matrix"){
    info_list<-list()
  }

  if(resolution =="cell"){
    if(!(group_markers_by %in% c("cell","celltype_cell","score_cell"))){
      group_markers_by<-"celltype_cell"
    }
  } else if(resolution =="cluster"){
    if(!(group_markers_by %in% c("cluster","celltype_cluster"))){
      group_markers_by<-"celltype_cluster"
      }
  }

  if ("cluster" %in% annotation_resolution){
    anno_cluster<-merge(final_dt_cluster, anno_dt_cl[,-"quantile_score_cluster"], by=c("seurat_clusters"))
    anno_dt_cell<-anno_cluster[order(-diff_score)][,head(.SD, 1),"cell"]

    anno_dt_cell_ptc<-anno_dt_cell[,ncell_celltype_cluster:= .N,by=c("seurat_clusters","annotation_per_cell")]
    anno_dt_cell_ptc[,ncell_tot_cluster:= .N, by="seurat_clusters"]
    anno_dt_cell_ptc[,perc_celltype_cluster:= round((ncell_celltype_cluster/ncell_tot_cluster)*100, digits = 2)]
    dt_top_marker<-unique(merge.data.table(dt_score,anno_dt_cell_ptc, by=c("cell_type","cell")))
    anno_dt_cl_rank<-unique(anno_dt_cell_ptc[,-c("cell","diff_score")])[order(-quantile_score_cluster)][,head(.SD, n_top_celltypes),"seurat_clusters"][,c("seurat_clusters","annotation_per_cell","quantile_score_cluster","ncell_tot_cluster","perc_celltype_cluster")]
    name<-paste0(annotation_name,"_per_cluster")
    colnames(anno_dt_cl_rank)<-c("seurat_clusters",name,"celltype_impact_score", "ncell_tot_cluster","perc_celltype_cluster")

    if(data_type == "seurat"){
      colnames(anno_dt_cl_rank)[colnames(anno_dt_cl_rank) == "seurat_clusters"] <- cluster_info
    } else{
      colnames(anno_dt_cl_rank)[colnames(anno_dt_cl_rank) == "seurat_clusters"] <- "cluster"
    }
    cluster_res_detailed_annotation_info<-list()
    cluster_res_detailed_annotation_info[["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes"]] <- as.data.table(anno_dt_cl_rank)


    if("cluster" %in% group_markers_by){

      # for each cluster retrieves first N markers
      dt_top<- unique(dt_top_marker[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","seurat_clusters")][,c("cell_type","marker","seurat_clusters","quantile_score_marker","weight","specificity","marker_type","annotation_per_cluster")])
      dt_top_marker_by_cl<-dt_top[order(seurat_clusters,-quantile_score_marker)][,head(.SD, n_top_markers),by=c("seurat_clusters")]
      dt_top_marker_by_cl<-as.data.table(dt_top_marker_by_cl)[,c("seurat_clusters","annotation_per_cluster","marker","marker_type","weight","specificity","quantile_score_marker")]
      colnames(dt_top_marker_by_cl)<-c("seurat_clusters",eval(name), "marker","marker_type","weight","specificity","gene_impact_score_per_cluster")

      #rename seurat_clusters columns using the original name in input
      if(data_type == "seurat"){
        colnames(dt_top_marker_by_cl)[colnames(dt_top_marker_by_cl) == "seurat_clusters"] <- cluster_info
      } else{
        colnames(dt_top_marker_by_cl)[colnames(dt_top_marker_by_cl) == "seurat_clusters"] <- "cluster"
      }

      cluster_res_detailed_annotation_info[["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_cluster"]] <- as.data.table(dt_top_marker_by_cl)
    }
    if("celltype_cluster" %in% group_markers_by){
      # for each cluster retrieves first N markers
      dt_top<- unique(dt_top_marker[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cluster")][,c("cell_type","marker","seurat_clusters","quantile_score_marker","weight","specificity","marker_type","annotation_per_cluster")])
      dt_top_marker_by_cl<-dt_top[order(seurat_clusters,-quantile_score_marker)][,head(.SD, n_top_markers),by=c("annotation_per_cluster")]
      dt_top_marker_by_cl<-as.data.table(dt_top_marker_by_cl)[,c("seurat_clusters","annotation_per_cluster","marker","marker_type","weight","specificity","quantile_score_marker")]
      colnames(dt_top_marker_by_cl)<-c("seurat_clusters",eval(name), "marker","marker_type","weight","specificity","gene_impact_score_per_celltype_cluster")

      #rename seurat_clusters columns using the original name in input
      if(data_type == "seurat"){
        colnames(dt_top_marker_by_cl)[colnames(dt_top_marker_by_cl) == "seurat_clusters"] <- cluster_info
      } else{
        colnames(dt_top_marker_by_cl)[colnames(dt_top_marker_by_cl) == "seurat_clusters"] <- "cluster"
      }

      cluster_res_detailed_annotation_info[["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cluster"]] <- as.data.table(dt_top_marker_by_cl)
    }

    if(data_type == "seurat"){
      data@misc[[annotation_name]]<-cluster_res_detailed_annotation_info
    } else{
      info_list[[annotation_name]]<-cluster_res_detailed_annotation_info
      data<-append(data,info_list)
    }
  }
  if ("cell" %in% annotation_resolution){

    dt_top_ct_by_cell<-final_dt[order(-diff_score)][,head(.SD, n_top_celltypes),cell]
    dt_top_ct_per_cell<-as.data.table(dt_top_ct_by_cell)[,c("cell","cell_type","diff_score")]

    name<-paste0(annotation_name,"_per_cell")
    name_score<-paste0(annotation_name,"_per_cell_score")

    colnames(dt_top_ct_per_cell)<-c("cell",eval(name),eval(name_score))

    cell_res_detailed_annotation_info<-list()
    cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_celltypes"]] <- as.data.table(dt_top_ct_per_cell)

    dt_top_marker<-unique(merge.data.table(dt_score,anno_dt_cell, by=c("cell_type","cell")))

    if("cell" %in% group_markers_by){
      dt_top <- unique(dt_top_marker[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","cell")][,c("cell_type","marker","marker_type","quantile_score_marker","cell","score","weight","specificity","annotation_per_cell", "diff_score")])
      #for each cell retrieves first N cell type and first N markers
      dt_top_marker_by_cell<-dt_top[order(-quantile_score_marker)][,head(.SD, n_top_markers),cell]

      dt_top_marker_per_cell<-as.data.table(dt_top_marker_by_cell)[,c("cell","annotation_per_cell","diff_score","marker","marker_type","weight","specificity","score")]
      colnames(dt_top_marker_per_cell)<-c("cell",eval(name),eval(name_score), "marker","marker_type","weight","specificity", "gene_impact_score_per_cell")

      cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_cell"]] <- as.data.table(dt_top_marker_per_cell)

    }
    if("celltype_cell" %in% group_markers_by){
      #for each cell retrieves first N cell type and first N markers
      if(!is.null(condition_group_info) & is.null(cell_type_group_info)){
        if(data_type == "seurat"){
          condition_table<-data@meta.data
          condition_table<-as.data.table(condition_table)[,cell:=rownames(condition_table)]
          col_vec<-c("cell",condition_group_info)
          condition_table<-condition_table[,..col_vec]
          colnames(condition_table)<-c("cell","condition")
          condition_table<-condition_table[,c("cell","condition")]

        } else{
          col_vec<-c("cell",condition_group_info, cell_type_group_info)
          condition_table<-as.data.table(condition_group_info)[,..col_vec]
          colnames(condition_table)<-c("cell","condition","celltype")
        }
        dt_top_marker_condition<-merge(dt_top_marker, condition_table, by="cell")
        dt_top <- unique(dt_top_marker_condition[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell","condition")][,c("condition","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top<-unique(dt_top[,c("annotation_per_cell","condition","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top[order(-quantile_score_marker)][,head(.SD, n_top_markers),c("annotation_per_cell","condition")]

        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), eval(condition_group_info),"marker","marker_type","gene_impact_score_per_celltype_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell"]] <- as.data.table(dt_top_marker_by_cell)
      }
      if(!is.null(condition_group_info) & !(is.null(cell_type_group_info))){
        if(data_type == "seurat"){
          condition_table<-data@meta.data
          condition_table<-as.data.table(condition_table)[,cell:=rownames(condition_table)]
          col_vec<-c("cell",condition_group_info, cell_type_group_info)
          condition_table<-condition_table[,..col_vec]
          colnames(condition_table)<-c("cell","condition","celltype")
          condition_table<-condition_table[,c("cell","condition","celltype")]
        } else{
          col_vec<-c("cell",condition_group_info, cell_type_group_info)
          condition_table<-as.data.table(condition_group_info)[,..col_vec]
          colnames(condition_table)<-c("cell","condition","celltype")

        }
        dt_top_marker_condition<-merge(dt_top_marker, condition_table, by="cell")
        dt_top <- unique(dt_top_marker_condition[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell","condition","celltype")][,c("condition","celltype","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top<-unique(dt_top[,c("annotation_per_cell","condition","celltype","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top[order(-quantile_score_marker)][,head(.SD, n_top_markers),c("annotation_per_cell","condition","celltype")]

        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), eval(condition_group_info),eval(cell_type_group_info),"marker","marker_type","gene_impact_score_per_celltype_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell"]] <- as.data.table(dt_top_marker_by_cell)
      } else {
        dt_top <- unique(dt_top_marker[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell")][,c("cell","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-unique(dt_top[,c("annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top_marker_by_cell[order(-quantile_score_marker)][,head(.SD, n_top_markers),annotation_per_cell]
        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), "marker","marker_type","gene_impact_score_per_celltype_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell"]] <- as.data.table(dt_top_marker_by_cell)
      }
    }
    if ("score_cell" %in% group_markers_by){
      if(!is.null(condition_group_info) & is.null(cell_type_group_info)){
          if(data_type == "seurat"){
            condition_table<-data@meta.data
            condition_table<-as.data.table(condition_table)[,cell:=rownames(condition_table)]
            col_vec<-c("cell",condition_group_info)
            condition_table<-condition_table[,..col_vec]
            colnames(condition_table)<-c("cell","condition")
            condition_table<-condition_table[,c("cell","condition")]

          } else{
            col_vec<-c("cell",condition_group_info, cell_type_group_info)
            condition_table<-as.data.table(condition_group_info)[,..col_vec]
            colnames(condition_table)<-c("cell","condition","celltype")
          }
        dt_top_marker_condition<-merge(dt_top_marker, condition_table, by="cell")
        # consider only cells with a score greater than the 90 percentile (as default)
        dt_top_marker_condition[,quantile_diff_score := quantile((unique(.SD, by = c("cell")))$diff_score, probs = top_cell_score_quantile_threshold, na.rm=TRUE),by="cell_type"]
        dt_top_marker_condition<-dt_top_marker_condition[diff_score >= quantile_diff_score]

        dt_top <- unique(dt_top_marker_condition[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell","condition")][,c("condition","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top<-unique(dt_top[,c("annotation_per_cell","condition","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top[order(-quantile_score_marker)][,head(.SD, n_top_markers),c("annotation_per_cell","condition")]

        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), eval(condition_group_info),"marker","marker_type","gene_impact_score_per_score_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_score_cell"]] <- as.data.table(dt_top_marker_by_cell)
      } else if(!is.null(condition_group_info) & !(is.null(cell_type_group_info))){
        if(data_type == "seurat"){
          condition_table<-data@meta.data
          condition_table<-as.data.table(condition_table)[,cell:=rownames(condition_table)]
          col_vec<-c("cell",condition_group_info, cell_type_group_info)
          condition_table<-condition_table[,..col_vec]
          colnames(condition_table)<-c("cell","condition","celltype")
          condition_table<-condition_table[,c("cell","condition","celltype")]
        } else{
          col_vec<-c("cell",condition_group_info, cell_type_group_info)
          condition_table<-as.data.table(condition_group_info)[,..col_vec]
          colnames(condition_table)<-c("cell","condition","celltype")

        }
        dt_top_marker_condition<-merge(dt_top_marker, condition_table, by="cell")
        dt_top_marker_condition[,quantile_diff_score := quantile((unique(.SD, by = c("cell")))$diff_score, probs = top_cell_score_quantile_threshold, na.rm=TRUE),by="cell_type"]
        dt_top_marker_condition<-dt_top_marker_condition[diff_score >= quantile_diff_score]

        dt_top <- unique(dt_top_marker_condition[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell","condition","celltype")][,c("condition","celltype","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top<-unique(dt_top[,c("annotation_per_cell","condition","celltype","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top[order(-quantile_score_marker)][,head(.SD, n_top_markers),c("annotation_per_cell","condition","celltype")]

        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), eval(condition_group_info),eval(cell_type_group_info),"marker","marker_type","gene_impact_score_per_score_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_score_cell"]] <- as.data.table(dt_top_marker_by_cell)
      } else if(is.null(condition_group_info) & is.null(cell_type_group_info)){
        dt_top_marker[,quantile_diff_score := quantile((unique(.SD, by = c("cell")))$diff_score, probs = top_cell_score_quantile_threshold, na.rm=TRUE),by="cell_type"]
        dt_top_marker<-dt_top_marker[diff_score >= quantile_diff_score]
        dt_top <- unique(dt_top_marker[, quantile_score_marker := quantile(score,probs = top_marker_score_quantile_threshold, na.rm=TRUE), by=c("marker","marker_type","annotation_per_cell")][,c("cell","annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-unique(dt_top[,c("annotation_per_cell","marker","marker_type","quantile_score_marker","weight","specificity")])
        dt_top_marker_by_cell<-dt_top_marker_by_cell[order(-quantile_score_marker)][,head(.SD, n_top_markers),annotation_per_cell]
        name<-paste0(annotation_name,"_per_cell")
        name_score<-paste0(annotation_name,"_per_cell_score")

        colnames(dt_top_marker_by_cell)<-c(eval(name), "marker","marker_type","gene_impact_score_per_score_cell","weight","specificity")

        cell_res_detailed_annotation_info[["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_score_cell"]] <- as.data.table(dt_top_marker_by_cell)
      }
    }
    if(data_type == "seurat"){
      if(is_empty(data@misc[[annotation_name]])){
        data@misc[[annotation_name]]<-cell_res_detailed_annotation_info
      } else {
        data@misc[[annotation_name]]<-append(data@misc[[annotation_name]], cell_res_detailed_annotation_info)
      }
    } else{
      if(is_empty(info_list[[annotation_name]])){
        info_list[[annotation_name]]<-cell_res_detailed_annotation_info
        data<-append(data,info_list)

      } else{
        info_list<-append(info_list,cell_res_detailed_annotation_info)
        data<-append(data,info_list)
      }
    }
  }
  return(data)

}
