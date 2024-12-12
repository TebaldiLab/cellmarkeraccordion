#' Automatically annotating and interpreting aberrant single-cell populations
#' with the built-in Cell Marker Accordion disease database
#'
#' This function identified aberrant populations exploiting the built-in
#' Accordion gene marker disease database. It takes in input either a Seurat
#' object or a raw or normalized count matrix and return in output the cell
#' types assignment and the detailed informations of the annotation results
#' (added to the Seurat object or as a list).
#'
#' @param data Either a  Seurat object (version 4 or 5) or a raw or normalized
#'   count matrix with genes on rows and cells on columns. If raw counts are
#'   provided, data are log-normalized exploiting the NormalizeData() function
#'   from the Seurat package.
#' @param disease Character string or character string vector specifying
#'   diseases to consider. Run the function "list_diseases()" to obtain the
#'   available diseases. If NULL, all diseases are considered. Default is NULL.
#' @param tissue Character string or character string vector specifying the tissue.
#'  Run the function "list_tissues()" to obtain the available tissues. If multiple
#'  tissues are selected cell types and markers from the selected tissues
#'  are aggregated. If NULL, all tissues are considered. Default is NULL.
#' @param include_descendants  Logical value indicating whether include all the
#'  tissues that are descendants of the selected tissue(s) according to the uberon
#'  ontology. If TRUE,cell types and markers from the selected tissues and their
#'  descendants are aggregated. Default is FALSE
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
#' @param NCIT_celltypes Character string or character string vector specifying the
#'   cell types to annotate. Run the function "list_aberrant_celltypes()" to obtain the
#'   available aberrant cell types. If this parameter is not specified, all aberrant
#'   cell types are used for the annotation. Default is NULL.
#' @param species Character string or character string vector specifying the
#'   species. Currently, either “Human” and/or “Mouse” are supported. If
#'   multiple species are selected, marker genes are merged together. Default is
#'   “Human”.
#' @param evidence_consistency_score_threshold Integer value
#'   specifying the minimum evidence consistency (EC) score for each
#'   marker. Only markers >= this threshold are kept. If NULL, no filter is
#'   applied. Default is NULL.
#' @param specificity_score_threshold numeric value in (0,1) specifying the
#'   minimum specificity score for each marker. Only markers <= this threshold
#'   are kept. If  NULL, no filter is applied. Default is NULL.
#' @param log2FC_threshold numeric value specifying the
#'   minimum log2FC threshold for each marker reporting this information.
#'   Only markers <= this threshold or without any log2FC
#'   are kept. If  NULL, no filter is applied. Default is NULL.
#' @param malignant_quantile_threshold numeric value in (0,1) specifying the malignant
#'  quantile threshold. If an aberant cell type is provide in input,
#'  a malignant score specific for each cell is computed.
#'  The \code{malignant_quantile_threshold} is
#'  computed across cells, and only cells with a score greater than
#'  the \code{top_cell_score_quantile_threshold} are labeled as malignant.
#'  If \code{annotation_resolution} is set to "cell", a cell
#' @param min_n_marker Integer value specifying the minimum number of markers to
#'   keep for each cell type. Only cell types with a number of markers >= this
#'   threshold are kept.  Default is 5.
#' @param max_n_marker Integer value specifying the maximum number of markers to
#'   keep for each cell type. For the selection, markers are ranked according to
#'   their combined score, obtained by multiplying evidence consistency score
#'   and specificity score. If  NULL, no filter is applied. Default is NULL.
#' @param combined_score_quantile_threshold numeric value in (0,1) specifying
#'   the combined score quantile threshold. For the selection, marers are
#'   ranked according to their combined score,  obtained by multiplying evidence
#'   consistency score and specificity score. Only markers >  the
#'   quantile_threshold are kept. If  NULL, no filter is applied. Default is
#'   NULL.
#' @param disease_vs_healthy Logical value indicating whether to compare the
#'   markers associated with disease cell types with respect to the markers
#'   associated with the corresponding healthy cell types. If TRUE the
#'   specificity score is calculated considering if a gene is also a marker for
#'   the corresponding healthy cell type. For each cell, a specific score for
#'   disease cell types and another score for the healthy types (named as the
#'   cell types label) are returned. Default is TRUE.
#' @param annotation_resolution Character string or character string vector
#'   specifying the resolution of the annotation. Either “cluster” and/or “cell”
#'   are supported. Default is “cluster”.
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
#' @param annotation_name Character string specifying the name of the column
#'   in either the metadata of the input Seurat object or in the input
#'   \code{cluster_info} where the annotation will be stored. Per cluster and
#'   per cell annotation results will be stored in the
#'   \code{annotation_name}_per_cluster and \code{annotation_name}_per_cell
#'   columns respectively.
#'   If \code{include_detailed_annotation_info} parameter is set to TRUE, the
#'   detailed information the stored in a list named \code{annotation_name}.
#'   Default is “accordion_disease”.
#' @param include_detailed_annotation_info Logical value indicating whether to
#'   store information on the top cell types and markers in the output. If TRUE,
#'   a nested list named \code{annotation_name} is created. If
#'   \code{resolution_annotation} is set to “cluster” and/or “cell, sublists
#'   named “cluster_resolution” and/or “cell_resolution” are then added. Inside
#'   the sublist “detailed_annotation_info” the \code{n_top_markers} markers,
#'   group by \code{group_markers_by} and the \code{n_top_celltypes} cell types
#'   are then included. If a Seurat object is provided as input the list is
#'   stored in the misc slot of the object (object@misc@\code{annotation_name}). If
#'   the input is a count matrix, the list is returned in the final output.
#'   Default is TRUE.
#' @param condition_group_info in case \code{data} is a Seurat object,
#'  \code{condition_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains condition ids for each cell;
#'  if \code{data} is a count matrix, \code{condition_group_info} should be need to be a
#'   data frame or data table containing condition identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “condition”, which specifies
#'   the condition id’s for each cell.  Default is NULL.
#' @param NCIT_celltype_group_info in case \code{data} is a Seurat object,
#'  \code{NCIT_celltype_group_info} should be need to be a character string specifying the
#'  name of the column in the metadata that contains cell types ids for each cell;
#'  if \code{data} is a count matrix, \code{NCIT_celltype_group_info} should be need to be a
#'   data frame or data table containing cell types identity for each cell. The
#'   data frame or data table should contain at least two columns, one  named
#'   “cell”, which specifies cell id’s, and one named “NCIT_celltype”, which specifies
#'   the cell types for each cell.  Default is NULL.
#' @param group_markers_by Character string or character string vector
#'   specifying the classification of marker genes. It possible to retrieve
#'   \code{n_top_markers} marker genes for each cell type identified with
#'   cluster ("celltype_cluster") or cell ("celltype_cell") resolution;
#'   \code{n_top_markers} marker genes per cluster ("cluster") or per cell
#'   ("cell") can be also obtained. Additionally, by setting \code{group_markers_by}
#'  to "score_cell", the \code{n_top_markers} marker genes only for
#'  cells with a score greater than \code{top_cell_score_quantile_threshold} are
#'  retrieved. Either "celltype_cluster", "celltype_cell",
#'  "cluster", "cell" or "score_cell". Default is "celltype_cluster".
#' @param top_cell_score_quantile_threshold numeric value in (0,1) specifying
#'  the cell score quantile threshold. For each cell type a score specific for
#'  each cell is computed. The \code{top_cell_score_quantile_threshold} is
#'  computed across cells belonging to the same cell type, and only
#'  cells with a score greater than the \code{top_cell_score_quantile_threshold}
#'  are kept. This parameter is necessary only when \code{group_markers_by}
#'  is set to "score_cell". Default is 0.90.
#' @param n_top_celltypes Integer value specifying the number of the top cell
#'   types to be included in the output for each cluster and cell depending on
#'   the selected \code{annotation_resolution} parameter Default is 5.
#' @param n_top_markers Integer value specifying the number of the top markers
#'   to be included in the output for each cell type, cluster or cell depending
#'   on the selected \code{annotation_resolution} and \code{group_markers_by}
#'   parameters. Default is 5.
#' @param top_marker_score_quantile_threshold numeric value in (0,1) specifying
#'   the marker score quantile threshold. For each marker a score specific for
#'   each cell is computed. To identify the \code{n_top_markers} for a cluster
#'   cl or a cell type ct, the \code{top_marker_score_quantile_threshold} is
#'   computed across cells belonging to that cluster or labeled as ct, and the
#'   \code{n_top_markers} with the maximum score are reported. Default is 0.75.
#' @param plot Logical value indicating whether to store plots displaying
#'   detailed annotation information.  This parameter can be set to TRUE only
#'   when \code{include_detailed_annotation_info} is set to TRUE. If TRUE,
#'   lollipop plots displaying the top \code{n_top_markers} group by
#'   \code{group_markers_by} and  op \code{n_top_celltypes} for each
#'   \code{annotation_resolution} together with the cell types hierarchies based
#'   on the cell ontology structure are stored in the "accordion" list. Default
#'   is TRUE.
#' @param color_by Character string specifying if the plot reporting the top
#' cell types for each cluster/cell is colored based on the assigned cell type
#' ("CL_celltype") or on cluster id ("cluster"). Default is "CL_celltype.
#'
#'
#' @return A Seurat object or a list
#' @details
#' If a Seurat object was provided in input, the function returns the Seurat
#' object with markers-based scaled data in the scale.data slot and cell types
#' annotation results in the metadata. If
#' \code{include_detailed_annotation_info} and \code{plot} were set to TRUE, a
#' list containing cell types and markers information, together with ggplot
#' objects, is stored in the “misc@\code{annotation_name}” slot. If a count matrix
#' was provided in input, the function returns a list containing the following
#' elements:
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
#' @importFrom stats aggregate
#' @importFrom methods as
#' @importFrom stats quantile
#' @export
accordion_disease<-function(data,
                            disease = NULL,
                            tissue = NULL,
                            cluster_info = "seurat_clusters",
                            assay = "RNA",
                            NCIT_celltypes = NULL,
                            species = "Human",
                            include_descendants = FALSE,
                            evidence_consistency_score_threshold = NULL,
                            specificity_score_threshold = NULL,
                            log2FC_threshold = NULL,
                            malignant_quantile_threshold = 0.95,
                            min_n_marker = 5,
                            max_n_marker = NULL,
                            combined_score_quantile_threshold = NULL,
                            disease_vs_healthy = TRUE,
                            annotation_resolution = "cluster",
                            cluster_score_quantile_threshold = 0.75,
                            allow_unknown = TRUE,
                            annotation_name = "accordion_disease",
                            include_detailed_annotation_info = TRUE,
                            condition_group_info = NULL,
                            NCIT_celltype_group_info = NULL,
                            group_markers_by = "celltype_cluster",
                            top_cell_score_quantile_threshold = 0.90,
                            n_top_celltypes = 5,
                            n_top_markers = 5,
                            top_marker_score_quantile_threshold = 0.75,
                            plot = TRUE,
                            color_by = "cell_type"
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

  data("disease_accordion_marker", package = "cellmarkeraccordion",envir = environment())
  disease_accordion_marker<-disease_accordion_marker[marker %in% rownames(data)]

  #load the Cell Marker Accordion database based on the condition selected
  #for those markers with log2FC keep only the genes with log2FC above the threshold selected
  if(!is.null(log2FC_threshold)){
    if(!is.numeric(log2FC_threshold)){
      warning("Invalid log2FC_threshold type. Parameter log2FC_threshold must be a numeric value. No filter is applied")
    } else{
      disease_accordion_marker<-disease_accordion_marker[is.na(log2FC)| log2FC >= log2FC_threshold]
    }
  }


  #select species
  #change name to the input species
  if(length(species) ==1){
    if(!(species %in% c("Human","Mouse"))){
      warning("Invalid species type")
      if(all(grepl("^[[:upper:]]+$", rownames(data)[1:10]))){
        disease_accordion_marker<-disease_accordion_marker[species %in% "Human"]
        warning("The dataset might be human. Human markers are indeed used.")
      } else{
        disease_accordion_marker<-disease_accordion_marker[species %in% "Mouse"]
        warning("The dataset might be mouse. Mouse markers are indeed used.")
      }
    } else if (species %in% c("Human","Mouse")){
      input_species<-species
      disease_accordion_marker<-disease_accordion_marker[species %in% input_species]
    }
    #if more than one species is selected aggregate genes and in case of common genes between the species the relative EC score are summed
  } else if(length(species) >=2){
    if(all(grepl("^[[:upper:]]+$", rownames(data)[1:10]))){ #convert to human
      disease_accordion_marker[,marker:= toupper(marker)] # convert lower case in upper case (human symbol)
    } else{
      disease_accordion_marker[,marker:= str_to_title(marker)] # convert upper case in lower case (mouse symbol)
    }
    input_species<-species
    ec_score<-unique(disease_accordion_marker[,c("NCIT_celltype","marker","marker_type","EC_score_NCIT_global")])
    ec_score[,EC_score_sum:= sum(EC_score_NCIT_global), by=c("NCIT_celltype","marker","marker_type")]
    disease_accordion_marker[,species:=paste(input_species,collapse=", ")]
    disease_accordion_marker<-merge(disease_accordion_marker,ec_score, by=c("NCIT_celltype","marker","marker_type"))
    disease_accordion_marker<-unique(disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score_sum","resource")])
    colnames(disease_accordion_marker)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score","resource")
  }

  #check disease selected
  if(is.null(disease)){
    warning("no specific disease selected. The entire disease Accordion database is considered")
  } else{
      disease_accordion_marker<-disease_accordion_marker[toupper(DO_diseasetype) %in% toupper(disease)]
      if(nrow(disease_accordion_marker) == 0){
        warning("disease not found. The entire disease Accordion database is considered")
        disease_accordion_marker<-disease_accordion_marker
      }
      if(length(disease)>1){
        disease_accordion_marker[, DO_ID:="ALL"]
        disease_accordion_marker[, DO_diseasetype:="ALL"]
      }
  }
  disease_accordion_marker<-disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score_NCIT_global", "resource")]
  colnames(disease_accordion_marker)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score","resource")

  #keep only tissue gives in input and in case all its descendants
  if(!is.null(tissue)){
    #check tissue name in input
    tissue<-tolower(tissue)
    tissue_not_in_accordion<-tissue[!(tissue %in% unique(disease_accordion_marker$Uberon_tissue))]
    accordion_tissue<-disease_accordion_marker[Uberon_tissue %in% tissue]
    if(length(tissue_not_in_accordion) == 1){
      if(uniqueN(accordion_tissue$Uberon_tissue) == 0){
        ct_not_present<-knitr::combine_words(tissue_not_in_accordion[1:length(tissue_not_in_accordion)])
        warning(eval(ct_not_present), " tissue is not present. Annotation performed considering all tissues in the database")
        tissue<-unique(disease_accordion_marker$Uberon_tissue)

      } else{
        ct_not_present<-knitr::combine_words(tissue_not_in_accordion[1:length(tissue_not_in_accordion)])
        ct_present<-knitr::combine_words(unique(accordion_tissue$Uberon_tissue))
        warning(eval(ct_not_present), " tissue is not present. Annotation performed considering only ", eval(ct_present), " tissues")
        disease_accordion_marker<-disease_accordion_marker[Uberon_tissue %in% ct_present]
        tissue<-ct_present

      }
    } else if(length(tissue_not_in_accordion) > 1) {
      if(uniqueN(accordion_tissue$Uberon_tissue) == 0){
        ct_not_present<-knitr::combine_words(tissue_not_in_accordion[1:length(tissue_not_in_accordion)])
        warning(eval(ct_not_present), " tissues are not present. Annotation performed considering all tissues in the database")
        tissue<-unique(disease_accordion_marker$Uberon_tissue)

      } else{
        ct_not_present<-knitr::combine_words(tissue_not_in_accordion[1:length(tissue_not_in_accordion)])
        ct_present<-knitr::combine_words(unique(accordion_tissue$Uberon_tissue))
        warning(eval(ct_not_present), " tissues are not present. Annotation performed considering only ", eval(ct_present), " tissues")
        disease_accordion_marker<-disease_accordion_marker[Uberon_tissue %in% ct_present]
        tissue<-ct_present

      }
    }

    if(include_descendants == TRUE){
      data("uberon_onto", package = "cellmarkeraccordion",envir = environment())
      root_id<-unique(disease_accordion_marker[Uberon_tissue %in% tissue]$Uberon_ID)
      desc<-as.data.table(get_descendants(uberon_onto, roots=eval(root_id)))
      disease_accordion_marker<-disease_accordion_marker[Uberon_ID %in% desc$V1]

      if(length(root_id) > 1){
        disease_accordion_marker[, Uberon_ID:=paste(root_id,collapse=", ")]
        disease_accordion_marker[, Uberon_tissue:=paste(tissue,collapse=", ")]

      } else{
        disease_accordion_marker[,Uberon_ID:=eval(root_id)]
        disease_accordion_marker[,Uberon_tissue:=eval(tissue)]
      }

    } else{
      disease_accordion_marker<-disease_accordion_marker[Uberon_tissue %in% tissue]
      if(length(tissue) > 1){
        ub_tissue<-unique(disease_accordion_marker$Uberon_tissue)
        ub_id<-unique(disease_accordion_marker$Uberon_ID)

        disease_accordion_marker[, Uberon_ID:=paste(ub_tissue,collapse=", ")]
        disease_accordion_marker[, Uberon_tissue:=paste(ub_id,collapse=", ")]
      }
    }

    #calculate EC score based on filtered tissue
    disease_accordion_marker<-unique(disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
    ec_score<-ddply(disease_accordion_marker,.(species,DO_diseasetype,DO_ID, Uberon_tissue,Uberon_ID,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
    colnames(ec_score)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","EC_score")
    disease_accordion_marker<-merge(disease_accordion_marker,ec_score,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type"), all.x = TRUE)

  } else{
    disease_accordion_marker[, Uberon_ID:="ALL"]
    disease_accordion_marker[, Uberon_tissue:="ALL"]
    disease_accordion_marker<-disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score", "resource")]
    colnames(disease_accordion_marker)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score","resource")

  }

  # keep only cell types gives in input and markers found in data
  # assigned the parameter NCIT_celltype to the more specific "input_NCIT_celltype"
  input_NCIT_celltype <- NCIT_celltypes

  if(is.null(input_NCIT_celltype)){
    disease_accordion_marker<-disease_accordion_marker[marker %in% rownames(data)]
  } else {
    input_NCIT_celltype_not_in_accordion<-input_NCIT_celltype[!(input_NCIT_celltype %in% unique(disease_accordion_marker$NCIT_celltype))]
    disease_present<-disease_accordion_marker[NCIT_celltype %in% input_NCIT_celltype][marker %in% rownames(data)]
    if(length(input_NCIT_celltype_not_in_accordion) == 1){
      if(uniqueN(disease_present$NCIT_celltype) == 0){
        ct_not_present<-knitr::combine_words(input_NCIT_celltype_not_in_accordion[1:length(input_NCIT_celltype_not_in_accordion)])
        warning(eval(ct_not_present), " cell type is not present. Annotation performed considering all aberrant cell types in the database")
        input_NCIT_celltype<-unique(disease_accordion_marker$NCIT_celltype)
      } else{
        ct_not_present<-knitr::combine_words(input_NCIT_celltype_not_in_accordion[1:length(input_NCIT_celltype_not_in_accordion)])
        ct_present<-knitr::combine_words(unique(disease_accordion_marker$NCIT_celltype))
        warning(eval(ct_not_present), " cell type are not present. Annotation performed considering only ", eval(ct_present), " cell types")
        disease_accordion_marker<-disease_accordion_marker[NCIT_celltype %in% ct_present]
        input_NCIT_celltype<-ct_present

      }
    } else if(length(input_NCIT_celltype_not_in_accordion) > 1) {
      if(uniqueN(disease_present$NCIT_celltype) == 0){
        ct_not_present<-knitr::combine_words(input_NCIT_celltype_not_in_accordion[1:length(input_NCIT_celltype_not_in_accordion)])
        warning(eval(ct_not_present), " cell types are not present. Annotation performed considering all aberrant cell types in the database")
        input_NCIT_celltype<-unique(disease_accordion_marker$NCIT_celltype)
      } else{
        ct_not_present<-knitr::combine_words(input_NCIT_celltype_not_in_accordion[1:length(input_NCIT_celltype_not_in_accordion)])
        ct_present<-knitr::combine_words(unique(disease_accordion_marker$NCIT_celltype))
        warning(eval(ct_not_present), " cell types are not present. Annotation performed considering only ", eval(ct_present), "  cell types")
        disease_accordion_marker<-disease_accordion_marker[NCIT_celltype %in% ct_present]
        input_NCIT_celltype<-ct_present

      }
    } else if(length(input_NCIT_celltype_not_in_accordion) ==0){
      disease_accordion_marker<-disease_accordion_marker[NCIT_celltype %in% input_NCIT_celltype]

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


  #penalized discordant markers (ex. CD38 positive and CD38 negative for B cell)
  disease_accordion_marker[,group:=paste0(NCIT_celltype, "_", marker)]
  dup<-unique(disease_accordion_marker[,-"resource"])
  discordant_marker<-dup[duplicated(dup[,c("group")])]

  if(nrow(discordant_marker) > 0){
    disease_accordion_marker[,EC_score:=ifelse(group %in% discordant_marker$group & marker_type == "negative", -EC_score, EC_score)]
    ec_score <- aggregate(EC_score~group,disease_accordion_marker[group %in% discordant_marker$group], sum)
    colnames(ec_score)<-c("group","EC_score_new")
    ec_score<-as.data.table(ec_score)[,marker_type_new:=ifelse(EC_score_new > 0, "positive","negative")]
    disease_accordion_marker<-merge(disease_accordion_marker, ec_score, by="group", all.x = T)
    disease_accordion_marker[,EC_score_new:=abs(EC_score_new)]
    disease_accordion_marker[,EC_score:=ifelse(group %in% discordant_marker$group, EC_score_new, EC_score)]
    disease_accordion_marker[,marker_type:=ifelse(group %in% discordant_marker$group, marker_type_new, marker_type)]
    disease_accordion_marker<-unique(disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type", "EC_score")])
  }

  #remove marker with EC = 0
  disease_accordion_marker<-disease_accordion_marker[EC_score > 0]
  disease_accordion_marker<-unique(disease_accordion_marker)

  disease_accordion_marker<-unique(disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type", "EC_score")])
  #keep only marker genes with EC score above the selected threshold
  if(!is.null(evidence_consistency_score_threshold)){
    if(!is.numeric(evidence_consistency_score_threshold) | !(evidence_consistency_score_threshold %in% 1 == 0)){
      warning("Invalid evidence_consistency_score_threshold type. Parameter evidence_consistency_score_threshold must be an integer value (currently in [1,16]). No filter is applied")
    } else{
      disease_accordion_marker<-disease_accordion_marker[EC_score >= evidence_consistency_score_threshold]
    }
  }


  #evidence consistency score log-transformed
  disease_accordion_marker[,EC_score_scaled := log10(EC_score)+1]

  #compute specificity for positive and negative markers
  mark_spec<-ddply(disease_accordion_marker,.(marker,marker_type),nrow)
  colnames(mark_spec)<-c("marker","marker_type","specificity")
  disease_accordion_marker<-merge(disease_accordion_marker,mark_spec,by=c("marker","marker_type"),all.x = TRUE)

  length_ct_pos<-uniqueN(disease_accordion_marker[marker_type=="positive"]$NCIT_celltype)
  length_ct_neg<-uniqueN(disease_accordion_marker[marker_type=="negative"]$NCIT_celltype)

  #scale and log transforme specificity
  disease_accordion_marker<-disease_accordion_marker[marker_type=="positive",specificity_scaled := scales::rescale(as.numeric(specificity), to = c(1,length_ct_pos),from = c(length_ct_pos,1))
  ][marker_type=="negative",specificity_scaled := scales::rescale(as.numeric(specificity), to = c(1,length_ct_neg),from = c(length_ct_neg,1))
  ][,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","marker","marker_type","specificity","specificity_scaled","EC_score_scaled","EC_score")]

  disease_accordion_marker<-disease_accordion_marker[marker_type=="positive",specificity_scaled := scales::rescale(as.numeric(specificity_scaled), to = c(min(disease_accordion_marker[marker_type=="positive"]$EC_score),max(disease_accordion_marker[marker_type=="positive"]$EC_score)))
  ][marker_type=="negative",specificity_scaled := scales::rescale(as.numeric(specificity_scaled), to = c(min(disease_accordion_marker[marker_type=="negative"]$EC_score),max(disease_accordion_marker[marker_type=="negative"]$EC_score)))
  ][,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","marker","marker_type","specificity","specificity_scaled","EC_score_scaled","EC_score")]

  disease_accordion_marker[,specificity:=as.numeric(format(round(1/specificity,2), nsmall=2))]
  disease_accordion_marker[,specificity_scaled:=log10(specificity_scaled)+1]


  #keep only marker genes with EC score above the selected threshold
  if(!is.null(evidence_consistency_score_threshold)){
    if(!is.numeric(evidence_consistency_score_threshold) | !is.integer(evidence_consistency_score_threshold)){
      warning("Invalid evidence_consistency_score_threshold type. Parameter evidence_consistency_score_threshold must be an integer value (currently in [1,7]). No filter is applied")
    } else{
      disease_accordion_marker<-disease_accordion_marker[EC_score >= evidence_consistency_score_threshold]
    }
  }

    if(disease_vs_healthy == T){   # compare the healthy and the disease cell types if compare is set to TRUE
      #disease_accordion_marker[,cellID_healthy:= tstrsplit(NCIT_ID, "-", keep=2)]
      data("accordion_marker", package = "cellmarkeraccordion",envir = environment())
      accordion_healthy<-accordion_marker
      accordion_healthy<-accordion_healthy[species %in% input_species & Uberon_tissue %in% tissue & CL_ID %in% unique(disease_accordion_marker$CL_ID) & marker %in% rownames(data)]
      accordion_healthy<-accordion_healthy[!is.na(CL_ID)]
      accordion_healthy[,DO_diseasetype:=NA][,DO_ID:=NA][,NCIT_celltype:=CL_celltype][,NCIT_ID:=CL_ID]
      accordion_healthy<-unique(accordion_healthy)
      if(nrow(accordion_healthy)>0){
        disease_accordion_marker<-rbind(disease_accordion_marker[,c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score")],accordion_healthy[,c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score")])
      }
    } else if (disease_vs_healthy == F){
      disease_accordion_marker<- disease_accordion_marker
    }


  setkey(disease_accordion_marker,marker,NCIT_celltype)

  # merge Z_scaled_dt and accordion table
  disease_accordion_marker[,combined_score := specificity_scaled * EC_score_scaled]

  # filter markers according to the quantile threshold set
  if(!is.null(combined_score_quantile_threshold)){
    if(!is.numeric(combined_score_quantile_threshold) | combined_score_quantile_threshold > 1 | combined_score_quantile_threshold <= 0){
      warning("Invalid combined_score_quantile_threshold type. Parameter combined_score_quantile_threshold must be a numeric value in (0,1]. No filter is applied")
    } else{
      disease_accordion_marker[,quantile_combined_score:=quantile(combined_score, probs= combined_score_quantile_threshold), by="NCIT_celltype"]
      disease_accordion_marker<-disease_accordion_marker[combined_score>quantile_combined_score]
    }
  }

  # keep only the max_n_marker genes for each cell type
  if(!is.null(max_n_marker)){
    if(!is.numeric(max_n_marker) | !(max_n_marker %in% 1 == 0)){
      warning("Invalid max_n_marker type. Parameter max_n_marker must be an integer value. No filter is applied")
    } else {
      disease_accordion_marker<-disease_accordion_marker[order(-combined_score)][,head(.SD, max_n_marker), by="NCIT_celltype"]
    }
  }

  # number of markers for each cell type
  disease_accordion_marker[,length:= .N, by="NCIT_celltype"]
  if(!is.null(min_n_marker)){
    if(!is.numeric(min_n_marker) | !(min_n_marker %in% 1 == 0)){
      warning("Invalid min_n_marker type. Parameter min_n_marker must be an integer value. No filter is applied")
    } else{
      disease_accordion_marker<-disease_accordion_marker[length >= min_n_marker]
    }
  }


  if(!is.numeric(cluster_score_quantile_threshold) | cluster_score_quantile_threshold >1 | cluster_score_quantile_threshold <= 0){
    warning("Invalid cluster_score_quantile_threshold type. Parameter cluster_score_quantile_threshold must be a numeric value in (0,1]. Default 0.75 will be used.")
    cluster_score_quantile_threshold <- 0.75
  }

  if(!is.numeric(top_marker_score_quantile_threshold) | top_marker_score_quantile_threshold >1 | top_marker_score_quantile_threshold <= 0){
    warning("Invalid top_marker_score_quantile_threshold type. Parameter top_marker_score_quantile_threshold must be a numeric value in (0,1]. Default 0.75 will be used.")
    top_marker_score_quantile_threshold <- 0.75
  }

  # store original scale.data slot if present
  if(sum(dim(GetAssayData(data, assay=assay, slot='scale.data')))!=0){
    orig.scale_data<-GetAssayData(data, assay=assay, slot='scale.data')
  }

  # scale data based on markers used for the annotation
  data<-ScaleData(data, features = unique(disease_accordion_marker$marker))
  Zscaled_data<-GetAssayData(data, assay=assay, slot='scale.data')
  Zscaled_data<-as.data.table(as.data.frame(Zscaled_data),keep.rownames = "marker")
  setkey(Zscaled_data, marker)
  Zscaled_m_data<-melt.data.table(Zscaled_data,id.vars = c("marker"))
  colnames(Zscaled_m_data)<-c("marker","cell","expr_scaled")

  # compute the score for each cell
  dt_score<-merge.data.table(Zscaled_m_data,disease_accordion_marker, by="marker",allow.cartesian = TRUE)
  dt_score[,score := expr_scaled * combined_score]
  dt_score_ct <- unique(dt_score[, c("NCIT_celltype", "cell")])
  setkey(dt_score, NCIT_celltype, cell, marker_type)
  sum_dt <- dt_score[data.table("NCIT_celltype" = rep(dt_score_ct$NCIT_celltype, each = 2),
                                "cell" = rep(dt_score_ct$cell, each = 2),
                                "marker_type" = c("positive", "negative")),
                     .(score= (sum(score)/(sqrt((sum(EC_score_scaled * specificity_scaled)))))), by = .EACHI]


  sum_dt<-unique(sum_dt)
  sum_dt[is.na(score), score := 0]
  final_dt <- sum_dt[marker_type == "positive"
  ][, diff_score := score - sum_dt[marker_type == "negative", score]
  ][, marker_type := NULL][,score := NULL]

 # annotation per cluster
    if ("cluster" %in% annotation_resolution){
      setkey(cluster_table, cell)
      final_dt_cluster<-merge.data.table(final_dt, cluster_table, by="cell")
        final_dt_cluster[, quantile_score_cluster:= quantile(diff_score,probs = cluster_score_quantile_threshold, na.rm=TRUE), by=c("seurat_clusters","NCIT_celltype")]
          if (allow_unknown == T){
            final_dt_cluster[quantile_score_cluster < 0, annotation_per_cell:= "unknown"][quantile_score_cluster > 0, annotation_per_cell := NCIT_celltype]
          } else {
            final_dt_cluster[, annotation_per_cell := NCIT_celltype]
          }
       # add the annotation results in the metadata of the Seurat data
        anno_dt_cell<-final_dt_cluster[order(-diff_score)][,head(.SD, 1),"cell"]
        anno_dt_cell_ptc<-anno_dt_cell[,ncell_celltype_cluster:= .N,by=c("seurat_clusters","annotation_per_cell")]
        anno_dt_cell_ptc[,ncell_tot_cluster:= .N, by="seurat_clusters"]
        anno_dt_cell_ptc[,perc_celltype_cluster:= round((ncell_celltype_cluster/ncell_tot_cluster)*100, digits = 2)]
        anno_dt_cl<-anno_dt_cell_ptc[order(-quantile_score_cluster)][,head(.SD, 1),"seurat_clusters"][,-c("cell","NCIT_celltype","diff_score","ncell_tot_cluster","ncell_celltype_cluster")]
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

          accordion_output<-list(GetAssayData(data, assay=assay, slot='scale.data'), cluster_table)
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

    if(uniqueN(anno_dt_cell$NCIT_celltype) ==1){
      thr_assig<-quantile(anno_dt_cell$diff_score, probs = malignant_quantile_threshold)
      if (allow_unknown == T){
        anno_dt_cell[diff_score < thr_assig, annotation_per_cell:= "unknown"][diff_score > thr_assig, annotation_per_cell := NCIT_celltype]

      } else {
        anno_dt_cell[, annotation_per_cell := NCIT_celltype]

      }
    } else{
      if (allow_unknown == T){
        anno_dt_cell[diff_score < 0, annotation_per_cell:= "unknown"][diff_score > 0, annotation_per_cell := NCIT_celltype]

      } else {
        anno_dt_cell[, annotation_per_cell := NCIT_celltype]

      }
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
                                                      cluster_info,
                                                      final_dt,
                                                      anno_dt_cell,
                                                      n_top_celltypes,
                                                      n_top_markers,
                                                      top_marker_score_quantile_threshold,
                                                      top_cell_score_quantile_threshold,
                                                      condition_group_info,
                                                      NCIT_celltype_group_info)
    } else{
      accordion_output<-include_detailed_annotation_info_helper(accordion_output,
                                                                data_type,
                                                                annotation_resolution,
                                                                final_dt_cluster,
                                                                anno_dt_cl,
                                                                dt_score,
                                                                annotation_name,
                                                                group_markers_by,
                                                                cluster_info,
                                                                final_dt,
                                                                anno_dt_cell,
                                                                n_top_celltypes,
                                                                n_top_markers,
                                                                top_marker_score_quantile_threshold,
                                                                top_cell_score_quantile_threshold,
                                                                condition_group_info,
                                                                NCIT_celltype_group_info)
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

