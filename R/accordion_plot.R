#' Interpreting annotation results
#'
#' This function generates lollipop plots displaying the detailed annotation
#' results obtained with the \code{accordion_annotation},
#' \code{accordion_disease_annotation} and \code{accordion_custom_annotation}
#' functions.
#'
#' top cell types (or pathways) and top markers It takes in input either a
#' Seurat object or a raw or normalized count matrix and a table of marker genes
#' associated to cell types or even to pathways and return in output the cell
#' types/pathways assignment and the detailed informations of the annotation
#' results (added to the Seurat object or as a list).
#'
#' @param data A Seurat object or a list containing “detailed_annotation_info”,
#'   from either accordion(), accordion_disease() or accordion_custom()
#'   functions, run with include_detailed_annotation_info parameter set to TRUE.
#' @param info_to_plot Character string or character string vector specifying
#'   the list from which extract the detailed annotation information, either
#'   “accordion”, “accordion_disease” or “accordion_custom”, for which returns
#'   the plot, either “accordion”, “accordion_disease” or “accordion_custom”.
#'   Default is “accordion”.
#' @param resolution Character string or character string vector specifying the
#'   annotation resolution for which provided the plots. Either “cluster” and/or
#'   “cell” are supported. Default is “cluster”.
#' @param group_markers_by Character string or character string vector
#'   specifying the classification of marker genes. It is possible to retrieve
#'   top marker genes for each cell type identified with cluster
#'   ("celltype_cluster") or cell (“celltype_cell”) resolution; top  marker
#'   genes per cluster ("cluster") or per cell ("cell") can be also obtained.
#'   Additionally, by setting \code{group_markers_by}
#'   to "score_cell", the \code{n_top_markers} marker genes only for
#'   cells with a score greater than \code{top_cell_score_quantile_threshold} are
#'   retrieved. Either "celltype_cluster", "celltype_cell",
#'   "cluster", "cell" or "score_cell". Default is "celltype_cluster".
#'
#' @return A Seurat object or a list.
#' @details If a Seurat object was provided in input, the function returns the
#' Seurat object with a list of ggplot objects added to the "misc" slot in the
#' \code{info_to_plot} list. If a list was provided in input, the function
#' returns the same list with the addition of the ggplot objects.
#' @import ontologyIndex
#' @import igraph
#' @import ontologyPlot
#' @import ontoProc
#' @import ggraph
#' @import cowplot
#' @import ggnewscale
#' @import purrr
#' @export
accordion_plot<-function(data,
                         info_to_plot = "accordion",
                         resolution = "cluster",
                         group_markers_by = "celltype_cluster"
                         ){


  bs<-22
  # check group_markers_by input
  if(!(resolution %in% c("cluster","cell"))){
    stop("invalid resolution Please select \"cluster\" or \"cell\"")
  } else{
    resolution_slot <- paste0(resolution,"_resolution")
  }

  # check group_markers_by input
  if(!(group_markers_by %in% c("cluster","celltype_cluster","cell","celltype_cell","score_cell"))){
    stop("invalid group_by. Please select \"cluster\",\"celltype_cluster\", \"cell\" or \"celltype_cell\"")
  } else if(resolution %in% "cell"){
    if(!(group_markers_by %in% c("cell","celltype_cell","score_cell"))){
      group_markers_by<-"celltype_cell"
    }
  }
  top_markers<-paste0("top_markers_per_",group_markers_by)


  # check di input data
  if(class(data) != "Seurat"){
    data_type <- "matrix"

    if(is_empty(data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]])){
      stop(paste0("detailed annotation info list for ", resolution, " resolution is empty"))
    } else {
      top_cell_type_dt<-data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][["top_celltypes"]]
      top_marker_dt<-data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[top_markers]]
    }
  } else {
    data_type <- "seurat"
    if(is_empty(data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]])){
      stop(paste0("detailed annotation info slot for ", resolution_slot, " resolution is empty"))
    } else{
      top_cell_type_dt<-data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][["top_celltypes"]]
      top_marker_dt<-data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[top_markers]]
    }
  }

  marker_slot_plot<-paste0(top_markers,"_plot")
  celltype_slot_plot<-"top_celltypes_plot"


    if(info_to_plot == "accordion"){
      data(cell_onto)
      ontology_celltype<-as.data.frame(cell_onto[["name"]])
      colnames(ontology_celltype)<-"cell_type"
      ontology_celltype$cell_ID<-rownames(ontology_celltype)
      ontology_celltype<-as.data.table(ontology_celltype)
      if("cluster" %in% resolution){
        top_celltypes<-merge(top_cell_type_dt,ontology_celltype,  by.x="accordion_per_cluster", by.y="cell_type")
      } else if ("cell" %in% resolution){
        top_celltypes<-merge(top_cell_type_dt,ontology_celltype,  by.x="accordion_per_cell", by.y="cell_type")
      }
    } else{
      top_celltypes<-top_cell_type_dt
    }
  cluster_column_name<-colnames(top_cell_type_dt)[1]
  info_to_plot_per_cluster<-paste0(info_to_plot, "_per_cluster")
  info_to_plot_per_cell<-paste0(info_to_plot, "_per_cell")

    if("cluster" %in% group_markers_by){
      group<-as.vector(unique(top_marker_dt[, get(cluster_column_name)]))
    } else if ("celltype_cluster" %in% group_markers_by){
      group<-as.vector(unique(top_marker_dt[, get(info_to_plot_per_cluster)]))
    } else if ("cell" %in% group_markers_by){
      group<-unique(top_marker_dt$cell)
    } else if( "celltype_cell" %in% group_markers_by){
      group<-as.vector(unique(top_marker_dt[, get(info_to_plot_per_cell)]))
    } else if("score_cell" %in% group_markers_by){
      group<-as.vector(unique(top_marker_dt[, get(info_to_plot_per_cell)]))
    }


      for (gr in group){
      #lolliplot with top N markers per cluster
        #lolliplot with top N markers per cluster
        if("cluster" %in% group_markers_by){
          top_dt_cl<-top_marker_dt[get(cluster_column_name) == gr]
          colnames(top_dt_cl)[colnames(top_dt_cl) == "gene_impact_score_per_cluster"] <- "impact_score"
          name<-paste0(gr,"_", as.vector(unique(top_dt_cl[, get(info_to_plot_per_cluster)])))

        } else if("celltype_cluster" %in% group_markers_by){
          top_dt_cl<-top_marker_dt[get(info_to_plot_per_cluster) == gr]
          colnames(top_dt_cl)[colnames(top_dt_cl) == "gene_impact_score_per_celltype_cluster"] <- "impact_score"
          name<-gr

        } else if("cell" %in% group_markers_by){
          top_dt_cl<-top_marker_dt[cell == gr]
          colnames(top_dt_cl)[colnames(top_dt_cl) == "gene_impact_score_per_cell"] <- "impact_score"
          name<-paste0(gr,"_", as.vector(unique(top_dt_cl[, get(info_to_plot_per_cell)])))


        } else if("celltype_cell" %in% group_markers_by){
          top_dt_cl<-top_marker_dt[get(info_to_plot_per_cell) == gr]
          colnames(top_dt_cl)[colnames(top_dt_cl) == "gene_impact_score_per_celltype_cell"] <- "impact_score"
          name<-gr

        } else if("score_cell" %in% group_markers_by){
          top_dt_cl<-top_marker_dt[get(info_to_plot_per_cell) == gr]
          colnames(top_dt_cl)[colnames(top_dt_cl) == "gene_impact_score_per_score_cell"] <- "impact_score"
          name<-gr
        }

        top_dt_cl<-top_dt_cl[order(impact_score)]
        top_dt_cl<-top_dt_cl[impact_score > 0]

        top_dt_cl[,marker:=factor(marker,levels = unique(marker))]
        top_dt_cl<-top_dt_cl[, marker_type:=factor(marker_type, levels=c("positive","negative"))]

        #filtering out negative quantile

        colfunc_pos <- colorRampPalette(c("#8B1A1A", "#f1d9d4"))
        vec_pos<-colfunc_pos(5)

        colfunc_neg <- colorRampPalette(c("#104E8B", "#cfdbe7"))
        vec_neg<-colfunc_neg(5)

        top_dt_cl[,specificity_ratio:=as.character(round(1/specificity, digits = 2))]

        top_dt_cl[specificity_ratio == 1,specificity_range:= "1"]
        top_dt_cl[specificity_ratio == 0.50,specificity_range:= "0.50"]
        top_dt_cl[specificity_ratio == 0.33,specificity_range:= "0.33"]
        top_dt_cl[specificity_ratio == 0.25,specificity_range:= "0.25"]
        top_dt_cl[specificity_ratio < 0.25,specificity_range:= "<0.25"]
        top_dt_cl[marker_type =="positive",specificity_positive := specificity_range]
        top_dt_cl[marker_type =="negative",specificity_negative := specificity_range]

        #add colors
        top_dt_cl[,color_combo:=as.character()]
        top_dt_cl[specificity_ratio == 1, color_combo:=ifelse(marker_type=="positive", vec_pos[1], vec_neg[1])]
        top_dt_cl[specificity_ratio == 0.50, color_combo:=ifelse(marker_type=="positive", vec_pos[2], vec_neg[2])]
        top_dt_cl[specificity_ratio == 0.33, color_combo:=ifelse(marker_type=="positive", vec_pos[3], vec_neg[3])]
        top_dt_cl[specificity_ratio == 0.25, color_combo:=ifelse(marker_type=="positive", vec_pos[4], vec_neg[4])]
        top_dt_cl[specificity_ratio < 0.25, color_combo:=ifelse(marker_type=="positive", vec_pos[5], vec_neg[5])]
        top_dt_cl<-top_dt_cl[order(-specificity_ratio)]

        vec_pos<-top_dt_cl[marker_type=="positive"]$color_combo
        names(vec_pos)<-top_dt_cl[marker_type=="positive"]$specificity_range
        vec_pos<-vec_pos[!duplicated(vec_pos)]

        vec_neg<-top_dt_cl[marker_type=="negative"]$color_combo
        names(vec_neg)<-top_dt_cl[marker_type=="negative"]$specificity_range
        vec_neg<-vec_neg[!duplicated(vec_neg)]

        if("accordion" %in% info_to_plot | "accordion_disease" %in% info_to_plot){
          top_dt_cl[, EC_score_range:= EC_score][EC_score > 5, EC_score_range:=5]
          pl <- ggplot(top_dt_cl, aes(impact_score, marker)) +
            geom_vline(xintercept = 0, linetype = 2) +
            geom_segment(aes(x = 0, xend = impact_score, y = marker, yend = marker, color=specificity_range),linewidth = bs/10, show.legend = F) +
            geom_point(aes(size=EC_score_range,color=specificity_range), alpha= 1, shape = 16) +
            theme_bw(base_size = bs) +
            scale_size(range=c(4,13), guide_legend(title="EC score"),breaks = c(1,2,3,4,5), limits=c(1,5))+
            scale_color_manual(values=vec_pos, breaks = names(vec_pos),  guide_legend(title="Specificity\n(positive)"), limits = force)+

            if(!is.na(unique(top_dt_cl$specificity_negative))){
              pl <- pl + new_scale("color") +
                new_scale("size") +
                geom_vline(xintercept = 0, linetype = 2) +
                geom_segment(aes(x = 0, xend = impact_score, y = marker, yend = marker, color=specificity_range),data = subset(top_dt_cl, !is.na(specificity_negative)),linewidth = bs/10, show.legend = F) +
                geom_point(aes(size=weight_range,color=specificity_range),data = subset(top_dt_cl, !is.na(specificity_negative)),alpha= 1, shape = 16) + #, stroke = NA
                theme_bw(base_size = bs) +
                scale_size(range=c(4,13), guide_legend(title="EC score"),breaks = c(1,2,3,4,5), limits=c(1,5))+
                scale_color_manual(values=vec_neg, breaks = names(vec_neg),  guide_legend(title="Specificity\n(negative)"))
            }

          pl <- pl + theme(panel.border = element_blank()) +
            labs(x = "Gene impact score", y = "") +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_text(size=bs),
                  text = element_text(size = bs),
                  axis.text.y = element_text(size = bs)) +
            theme(legend.position = "right", legend.margin = margin(10,0,0,0), legend.box.margin = margin(-5,-5,-5,5)) +
            theme(legend.text = element_text(margin = margin(l = -7, unit = "pt")), legend.key.size = unit(1.2,"line")) +
            ggtitle(name)



        } else {
          top_dt_cl[, weight_range:= weight][weight > 5, weight_range:=5]
          pl <- ggplot(top_dt_cl, aes(impact_score, marker)) +
            geom_vline(xintercept = 0, linetype = 2) +
            geom_segment(aes(x = 0, xend = impact_score, y = marker, yend = marker, color=specificity_range),linewidth = bs/10, show.legend = F) +
            geom_point(aes(size=weight_range,color=specificity_range), alpha= 1, shape = 16) +
            theme_bw(base_size = bs) +
            scale_size(range=c(4,13), guide_legend(title="Weight"),breaks = c(1,2,3,4,5), limits=c(1,5))+
            scale_color_manual(values=vec_pos, breaks = names(vec_pos),  guide_legend(title="Specificity\n(positive)"))

            if(!is.na(unique(top_dt_cl$specificity_negative))){
            pl <- pl + new_scale("color") +
                      new_scale("size") +
                      geom_vline(xintercept = 0, linetype = 2) +
                      geom_segment(aes(x = 0, xend = impact_score, y = marker, yend = marker, color=specificity_range),data = subset(top_dt_cl, !is.na(specificity_negative)),linewidth = bs/10, show.legend = F) +
                      geom_point(aes(size=weight_range,color=specificity_range),data = subset(top_dt_cl, !is.na(specificity_negative)),alpha= 1, shape = 16) + #, stroke = NA
                      theme_bw(base_size = bs) +
                      scale_size(range=c(4,13), guide_legend(title="Weight"),breaks = c(1,2,3,4,5), limits=c(1,5))+
                      scale_color_manual(values=vec_neg, breaks = names(vec_neg),  guide_legend(title="Specificity\n(negative)"))
            }

            pl <- pl + theme(panel.border = element_blank()) +
            labs(x = "Gene impact score", y = "") +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_text(size=bs),
                  text = element_text(size = bs),
                  axis.text.y = element_text(size = bs)) +
            theme(legend.position = "right", legend.margin = margin(10,0,0,0), legend.box.margin = margin(-5,-5,-5,5)) +
            theme(legend.text = element_text(margin = margin(l = -7, unit = "pt")), legend.key.size = unit(1.2,"line")) +
            ggtitle(name)
        }

        if(uniqueN(top_dt_cl$specificity_range) == 1){
          if(unique(top_dt_cl$specificity_range) == 1){
            pl <- pl + guides(color="none")
          }
        }
        if(uniqueN(top_dt_cl$weight_range) == 1){
          if(unique(top_dt_cl$weight_range == 1)){
            pl <- pl + guides(size="none")
          }
        }

        if("condition" %in% colnames(top_dt_cl)){
          col<-hue_pal()(uniqueN(top_dt_cl))
          pl<- pl + facet_grid(condition ~ .)
        }

        if(data_type == "seurat"){
            data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[marker_slot_plot]][[name]]<-pl
        } else {
          data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[marker_slot_plot]][[name]]<-pl

        }
      }
        if(!"cell" %in% group_markers_by){
          #global plot markers
          if("cluster" %in% group_markers_by){
            colnames(top_marker_dt)[colnames(top_marker_dt) == "gene_impact_score_per_cluster"] <- "impact_score"
            colnames(top_marker_dt)[colnames(top_marker_dt) == eval(cluster_column_name)] <- "group"

          } else if("celltype_cluster" %in% group_markers_by){
            colnames(top_marker_dt)[colnames(top_marker_dt) == "gene_impact_score_per_celltype_cluster"] <- "impact_score"
            colnames(top_marker_dt)[colnames(top_marker_dt) == eval(info_to_plot_per_cluster)] <- "group"

            } else if ("score_cell" %in% group_markers_by){
              colnames(top_marker_dt)[colnames(top_marker_dt) == "gene_impact_score_per_score_cell"] <- "impact_score"
              colnames(top_marker_dt)[colnames(top_marker_dt) == eval(info_to_plot_per_cell)] <- "group"

            } else if ("celltype_cell" %in% group_markers_by){
              colnames(top_marker_dt)[colnames(top_marker_dt) == "gene_impact_score_per_score_cell"] <- "impact_score"
              colnames(top_marker_dt)[colnames(top_marker_dt) == eval(info_to_plot_per_cell)] <- "group"
            }

          top_marker_dt<-top_marker_dt[order(group)]
          top_marker_dt<-top_marker_dt[,marker:=factor(marker,levels=unique(marker))]

          if(! ("condition" %in% colnames(top_marker_dt))){

              dotplot<- ggplot(top_marker_dt, aes(x=group, y = marker, color = group, size = impact_score)) +
                geom_point() +
                theme_bw(base_size = bs) +
                scale_size(range=c(4,10))+
                guides(colour="none")+
                theme(panel.border = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(size=bs),
                      text = element_text(size = bs),
                      axis.text.y = element_text(size = bs),
                      legend.position = "right", legend.margin = margin(10,0,0,0), legend.box.margin = margin(-5,-5,-5,5),
                      legend.text = element_text(margin = margin(l = -7, unit = "pt")), legend.key.size = unit(1.2,"line"),
                      axis.text.x = element_text(angle = 45, hjust=1))+
                scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "_"),
                                                               width = 20))
          } else {

            dotplot<- ggplot(top_marker_dt, aes(x=group, y = marker, color = group, alpha= condition, size = impact_score, group=condition)) +
              geom_point(position=position_dodge(width=0.3)) +
              scale_alpha_discrete(range = c(0.6, 1))  +
              theme_bw(base_size = bs) +
              scale_size(range=c(4,10))+
              guides(colour="none")+
              theme(panel.border = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.ticks.y = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_text(size=bs),
                    text = element_text(size = bs),
                    axis.text.y = element_text(size = bs),
                    legend.position = "right", legend.margin = margin(10,0,0,0), legend.box.margin = margin(-5,-5,-5,5),
                    legend.text = element_text(margin = margin(l = -7, unit = "pt")), legend.key.size = unit(1.2,"line"),
                    axis.text.x = element_text(angle = 45, hjust=1))+
              scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "_"),
                                                             width = 20))
          }

          if(data_type == "seurat"){
            data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[marker_slot_plot]][["global"]]<-dotplot
          } else {
            data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[marker_slot_plot]][["global"]]<-dotplot

          }
        }

        if("cluster" %in% resolution){
          group<-as.vector(unique(top_celltypes[, get(cluster_column_name)]))
          #global cell types
          colnames(top_celltypes)[colnames(top_celltypes) == "celltype_impact_score"] <- "impact_score"
          colnames(top_celltypes)[colnames(top_celltypes) == eval(cluster_column_name)] <- "group"
          colnames(top_celltypes)[colnames(top_celltypes) == eval(info_to_plot_per_cluster)] <- "cell_type"

          top_celltypes<-top_celltypes[order(group)]
          top_celltypes<-top_celltypes[,cell_type:=factor(cell_type,levels=unique(cell_type))]

          dotplot_ct<- ggplot(top_celltypes, aes(x=group, y = cell_type, color = group, size = impact_score)) +
            geom_point() +
            theme_bw(base_size = bs) +
            scale_size(range=c(4,10))+
            guides(colour="none")+
            theme(panel.border = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_text(size=bs),
                  text = element_text(size = bs),
                  axis.text.y = element_text(size = bs),
                  legend.position = "right", legend.margin = margin(10,0,0,0), legend.box.margin = margin(-5,-5,-5,5),
                  legend.text = element_text(margin = margin(l = -7, unit = "pt")), legend.key.size = unit(1.2,"line"),
                  axis.text.x = element_text(angle = 45, hjust=1))+
            scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "_"),
                                                           width = 20))

          if(data_type == "seurat"){
            data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][["global"]]<-dotplot_ct
          } else {
            data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][["global"]]<-dotplot_ct

          }



        } else if ("cell" %in% resolution){
          group<-unique(top_celltypes$cell)
        }

        for (gr in group){
          #lolliplot with top N cell types per cluster

           if("cluster" %in% resolution){
              top_celltypes_cl<-top_celltypes[group == gr]
              colnames(top_celltypes_cl)[colnames(top_celltypes_cl) == "celltype_impact_score"] <- "impact_score"
              colnames(top_celltypes_cl)[colnames(top_celltypes_cl) == eval(info_to_plot_per_cluster)] <- "cell_type"
              top_celltypes_cl<-top_celltypes_cl[order(-impact_score)]
              name<-paste0(gr,"_", unique(top_celltypes_cl$cell_type[1]))

            top_celltypes_cl<-top_celltypes_cl[order(impact_score)]
            top_celltypes_cl[,cell_type:=factor(cell_type,levels = unique(cell_type))]

            if(info_to_plot == "accordion" & length(top_celltypes_cl$cell_ID) > 1){
              onto_plot<-onto_plot2(cell_onto, top_celltypes_cl$cell_ID)
              onto_plot@nodes<-gsub("(.{10,}?)\\s", "\\1\n", onto_plot@nodes, perl = TRUE)

              onto_igraph<-graph_from_graphnel(onto_plot, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
              V(onto_igraph)$CL<- str_split_i(onto_plot@nodes, "CL:", i= -1)
              V(onto_igraph)[V(onto_igraph)$CL %in% str_split_i(top_celltypes_cl$cell_ID, "CL:", i= -1)]$color <- "#8B1A1A"
              V(onto_igraph)[!(V(onto_igraph)$CL %in% str_split_i(top_celltypes_cl$cell_ID, "CL:", i= -1))]$color <- "gray50"

                pl <- ggplot(top_celltypes_cl, aes(impact_score, cell_type)) +
                  geom_vline(xintercept = 0, linetype = 2) +
                  geom_segment(aes(x = 0, xend = impact_score, y = cell_type, yend = cell_type),color = "#8B1A1A",linewidth = bs/10, show.legend = F) +
                  geom_point(aes(size=perc_celltype_cluster), color="#8B1A1A" , alpha = 1, shape = 16) +
                  geom_text(aes(label = paste0(perc_celltype_cluster, "%")), size= bs*0.3, color="#8B1A1A", position = position_dodge(width = .9), vjust = 0.5,hjust= -0.3,show.legend = FALSE)+
                  theme_bw(base_size = bs) +
                  scale_size(range=c(4,10))+
                  theme(panel.border = element_blank()) +
                  labs(x = "Cell type impact score", y = "") +
                  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                                                 width = 20))+
                  theme(panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.ticks.y = element_blank(),
                        strip.background = element_blank(),
                        strip.text = element_text(size=bs),
                        text = element_text(size = bs),
                        plot.title = element_blank(), axis.text.y = element_text(size = bs),legend.position = "none",
                        plot.margin = unit(c(1,1,1,1), "cm"))+
                  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
                  scale_x_continuous(limits = c(0, round_any(max(top_celltypes_cl$impact_score)*1.25+2, 1, f = ceiling)))


                tree_plot <- ggraph(onto_igraph,layout = 'tree') +
                  geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name,padding = margin(10, 10, 10, 10, "mm"))),
                                 arrow = arrow(type = "closed", length = unit(3, 'mm')))+
                  geom_label(aes(x = x, y = y, label = name), nudge_y = 0.1, label.size = NA,size = bs*0.28, colour=V(onto_igraph)$color, label.padding = unit(0.1, "lines")) +
                  theme_graph()

                como_plot<-plot_grid(pl, tree_plot, rel_widths = c(1, 2), nrow=1, scale = c(1, 1))
                if(data_type == "seurat"){
                  data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][[name]] <- como_plot
                } else{
                  data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][[name]]<-como_plot

                }
            } else {
              pl <- ggplot(top_celltypes_cl, aes(impact_score, cell_type)) +
                geom_vline(xintercept = 0, linetype = 2) +
                geom_segment(aes(x = 0, xend = impact_score, y = cell_type, yend = cell_type),color = "#8B1A1A",linewidth = bs/10, show.legend = F) +
                geom_point(aes(size=perc_celltype_cluster), color="#8B1A1A" , alpha = 1, shape = 16) +
                geom_text(aes(label = paste0(perc_celltype_cluster, "%")), size= bs*0.3, color="#8B1A1A", position = position_dodge(width = .9), vjust = 0.5,hjust= -0.3,show.legend = FALSE)+
                theme_bw(base_size = bs) +
                scale_size(range=c(4,10))+
                theme(panel.border = element_blank()) +
                labs(x = "Cell type impact score", y = "") +
                scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                                               width = 20))+
                theme(panel.grid.major.y = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks.y = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(size=bs),
                      text = element_text(size = bs),
                      plot.title = element_blank(), axis.text.y = element_text(size = bs),legend.position = "none",
                      plot.margin = unit(c(1,1,1,1), "cm"))+
                theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
                scale_x_continuous(limits = c(0, round_any(max(top_celltypes_cl$impact_score)*1.25+2, 1, f = ceiling)))
              if(data_type == "seurat"){
                data@misc[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][[name]] <- pl
              } else{
                data[[info_to_plot]][[resolution_slot]][["detailed_annotation_info"]][[celltype_slot_plot]][[name]]<-pl

              }
            }

          }
        }



  return(data)

}
