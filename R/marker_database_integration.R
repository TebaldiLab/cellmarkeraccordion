#' Create an augmented database by integrated custom set of markers genes with
#' the Accordion database
#'
#' This function performs markers integration provided by the users and the
#' accordion database, either physiological ("healthy") or disease ("disease). It
#' takes in input a table of marker genes associated to cell types and return
#' the integrated database of custom and accordion markers. This new integrated
#' database can be used as input for the "accordion" or "disease_accordion" function
#' to provide the annotation.
#'@param marker_table Data table or data frame containing cell type
#'  markers. The table needs to have at least two columns, the
#'  \code{celltype_column},  which specifies cell types, and the
#'  \code{marker_column}, which specifies the corresponding markers on each row.
#'  Columns indicating the marker type (either positive or negative), species, condition,
#'  tissue and resource can be optionally included.
#'@param database String characters specifying the accordion database to use for
#'  integration. Either “healthy” or “disease” are supported, where "healthy" refers
#'  to the physiological accordion database, while "disease" refers to the disease
#'  accordion database. Default is "healthy"
#'@param species_column String characters specifying the name of the
#'  \code{marker_table} column containing species.
#'  Default is “species”.
#'@param disease_column String characters specifying the name of the
#'  \code{marker_table} column containing diseases.
#'  Default is “disease”.
#'@param tissue_column String characters specifying the name of the
#'  \code{marker_table} column containing tissues.
#'  Default is “tissue”.
#'@param celltype_column String characters specifying the name of the
#'  \code{marker_table} column containing cell types.
#'  Default is “cell_type”.
#'@param marker_column String characters specifying the name of the
#'  \code{marker_table} column containing markers.
#'  Default is “marker”.
#'@param marker_type_column String characters specifying the name of the
#'  \code{marker_table} column containing marker types.
#'  Default is “marker_type”.
#'@param marker_type_column String characters specifying the name of the
#'  \code{marker_table} column containing resources.
#'  Default is “resource”.
#' @importFrom plyr ddply
#' @import data.table
#' @importFrom stringr str_remove_all
#' @export
#' @return A data table
marker_database_integration<-function(marker_table,
                           database = "healthy",
                           species_column = "species",
                           disease_column = "disease",
                           tissue_column = "tissue",
                           celltype_column = "cell_type",
                           marker_column = "marker",
                           marker_type_column = "marker_type",
                           resource_column = "resource"

){
  if(database == "healthy"){
    data("accordion_marker", package = "cellmarkeraccordion",envir = environment())
  } else if(database == "disease"){
    data("disease_accordion_marker", package = "cellmarkeraccordion",envir = environment())
  } else {
    stop("The selected database does not exist. Please select one database between: \"healthy\" or \"disease\"")
  }

  marker_table<-as.data.table(marker_table)
  #healthy integration
  if(database == "healthy"){
    if(!celltype_column %in% colnames(marker_table)){
      stop("\"celltype_column\" not found. Please specificy the \"celltype_column\" contaning the list of cell types to integrate." )
    } else{
      setnames(marker_table, eval(celltype_column), "CL_celltype")

    }
    if(!marker_column %in% colnames(marker_table)){
      stop("\"marker_column\" not found.  Please specificy the \"marker_column\" contaning the list of markers to integrate.")
    }else{
      setnames(marker_table, eval(marker_column), "marker")

    }
    #check input table
    if(!species_column %in% colnames(marker_table)){
      warning("\"species_column\" not found. By default Human will be considered.")
      marker_table[,species:="Human"]
    } else{
      setnames(marker_table, eval(species_column), "species")
    }

    if(!tissue_column %in% colnames(marker_table)){
      warning("\"tissue_column\" not found. The integration will not consider the tissues.")
      marker_table[,Uberon_tissue:=NA][,original_tissue:=NA]
      accordion_marker[,Uberon_tissue:=NA][,Uberon_ID:=NA][,original_tissue:=NA]
      accordion_marker<-unique(accordion_marker)
    } else{
      setnames(marker_table, eval(tissue_column), "Uberon_tissue")
      marker_table[,original_tissue:=Uberon_tissue]

    }
    if(!marker_type_column %in% colnames(marker_table)){
      warning("\"marker_type_column\" not found. By default genes will be considered as positive markers.")
      marker_table[,marker_type:="positive"]
    } else{
      setnames(marker_table, eval(marker_type_column), "marker_type")

    }
    if(!resource_column %in% colnames(marker_table)){
      warning("\"resource_column\" not found. By default a unique resource will be associated to possible recurrence markers.")
        marker_table[,resource:="custom_set"]
    } else{
      setnames(marker_table, eval(resource_column), "resource")

    }

    marker_table[,original_celltype:=CL_celltype]
    marker_table[,log2FC:=NA]
    marker_table[,p.value:=NA]
    marker_table[,adjusted_p.value:=NA]
    marker_table[,pct1:=NA]


    accordion_marker$SPs_tissue_specific<-NULL
    accordion_marker$SPs_global<-NULL
    accordion_marker$ECs_tissue_specific<-NULL
    accordion_marker$ECs_global<-NULL

    marker_table<-unique(marker_table)

    #add cell_ID
    data("cell_onto", package = "cellmarkeraccordion",envir = environment())

    ontology_celltype<-as.data.frame(cell_onto[["name"]])
    colnames(ontology_celltype)<-"CL_celltype"
    ontology_celltype$CL_ID<-rownames(ontology_celltype)
    ontology_celltype$cell_definition<-cell_onto[["def"]]
    ontology_celltype<-as.data.table(ontology_celltype)

    marker_table<-merge(marker_table, ontology_celltype, by="CL_celltype", all.x=TRUE)

    #add Uberon_ID
    data("uberon_onto", package = "cellmarkeraccordion",envir = environment())
    ontology_tissue<-as.data.frame(uberon_onto[["name"]])
    colnames(ontology_tissue)<-"Uberon_tissue"
    ontology_tissue$Uberon_ID<-rownames(ontology_tissue)
    ontology_tissue$tissue_definition<-uberon_onto[["def"]]
    ontology_tissue<-as.data.table(ontology_tissue)
    marker_table<-merge(marker_table, ontology_tissue, by="Uberon_tissue", all.x=TRUE)



    #add gene description
    load("C:/Users/emmab/Desktop/PhD/CellMarkerAccordion_Rpackage/github/R2/data/gene_description.rda")
    marker_table<-merge(marker_table, gene_description, by.x="marker",by.y="gene_symbol")


    marker_table<-marker_table[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1")]
    marker_table<-unique(marker_table)
    marker_table[,cell_definition:=str_remove_all(cell_definition, '"')]
    marker_table[,cell_definition:=str_remove_all(cell_definition, "\\$")]
    marker_table[,cell_definition:=str_remove_all(cell_definition, r"(\\)")]
    marker_table[,cell_definition:=tstrsplit(cell_definition, "[", fixed = TRUE, keep = 1)]
    marker_table[,gene_description:=tstrsplit(gene_description,"[",fixed=TRUE,keep=1)]
    marker_table[,tissue_definition:=str_remove_all(tissue_definition, '"')]
    marker_table[,tissue_definition:=str_remove_all(tissue_definition, "\\$")]
    marker_table[,tissue_definition:=str_remove_all(tissue_definition, r"(\\)")]
    marker_table[,tissue_definition:=tstrsplit(tissue_definition, "[", fixed = TRUE, keep = 1)]
    marker_database<-rbind(accordion_marker, marker_table)

    #Compute score

      #calculate EC score considering tissues
      table<-unique(marker_database[,c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker","marker_type","resource")])
      ECs<-ddply(table,.(species,Uberon_tissue,Uberon_ID,CL_celltype,CL_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker", "marker_type","ECs_tissue_specific")
      marker_database<-merge(marker_database,ECs,by=c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker","marker_type"), all.x = TRUE)

      #calculate EC score not tissues
      table<-unique(marker_database[,c("species","CL_celltype","CL_ID","marker","marker_type","resource")])
      ECs<-ddply(table,.(species,CL_celltype,CL_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","CL_celltype","CL_ID","marker", "marker_type","ECs_global")
      marker_database<-merge(marker_database,ECs,by=c("species","CL_celltype","CL_ID","marker","marker_type"), all.x = TRUE)


      #Compute specificity
      table<-unique(marker_database[,c("species","CL_celltype","CL_ID","marker","marker_type")])
      mark_spec<-ddply(table,.(species, marker, marker_type),nrow)
      colnames(mark_spec)<-c("species","marker","marker_type","SPs_global")
      marker_database<-merge(marker_database,mark_spec,by=c("species","marker","marker_type"),all.x = TRUE)
      marker_database[,SPs_global:=format(round(1/SPs_global,2), nsmall=2)]

      #Compute specificity: Tissue specific
      table<-unique(marker_database[,c("Uberon_tissue","Uberon_ID","species","CL_celltype","CL_ID","marker","marker_type","resource")])
      mark_spec<-ddply(table,.(species, Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
      colnames(mark_spec)<-c("species","Uberon_tissue","Uberon_ID","marker","marker_type", "SPs_tissue_specific")
      marker_database<-merge(marker_database,mark_spec,by=c("species","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
      marker_database[,SPs_tissue_specific:=format(round(1/SPs_tissue_specific,2), nsmall=2)]

      marker_database<-marker_database[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1","ECs_global","ECs_tissue_specific","SPs_global","SPs_tissue_specific")]


  } else if(database == "disease"){ #disease integration
      if(!celltype_column %in% colnames(marker_table)){
        stop("\"celltype_column\" not found. Please specificy the \"celltype_column\" contaning the list of cell types to integrate." )
      } else{
        setnames(marker_table, eval(celltype_column), "NCIT_celltype")
      }
      if(!marker_column %in% colnames(marker_table)){
        stop("\"marker_column\" not found.  Please specificy the \"marker_column\" contaning the list of markers to integrate.")
      }else{
        setnames(marker_table, eval(marker_column), "marker")

      }
      if(!disease_column %in% colnames(marker_table)){
        warning("\"disease_column\" not found. The integration will not consider the disease.")
        marker_table[,DO_diseasetype:=NA][,DO_ID:=NA][,original_diseasetype:=NA][,DO_definition:=NA]
        disease_accordion_marker[,DO_diseasetype:=NA][,DO_ID:=NA][,original_diseasetype:=NA][,DO_definition:=NA]
        disease_accordion_marker<-unique(disease_accordion_marker)
      } else{
        setnames(marker_table, eval(disease_column), "DO_diseasetype")
        marker_table[,original_diseasetype:=DO_diseasetype][,DO_ID:=NA][,DO_definition:=NA]
      }
      if(!species_column %in% colnames(marker_table)){
        warning("\"species_column\" not found. By default Human will be considered.")
        marker_table[,species:="Human"]
      } else{
        setnames(marker_table, eval(species_column), "species")
      }

      if(!tissue_column %in% colnames(marker_table)){
        warning("\"tissue_column\" not found. The integration will not consider the tissues.")
        marker_table[,Uberon_tissue:=NA][,original_tissue:=NA]
        disease_accordion_marker[,Uberon_tissue:=NA][,Uberon_ID:=NA][,original_tissue:=NA]
        disease_accordion_marker<-unique(disease_accordion_marker)
      } else{
        setnames(marker_table, eval(tissue_column), "Uberon_tissue")
        marker_table[,original_tissue:=Uberon_tissue]

      }
      if(!marker_type_column %in% colnames(marker_table)){
        warning("\"marker_type_column\" not found. By default genes will be considered as positive markers.")
        marker_table[,marker_type:="positive"]
      } else{
        setnames(marker_table, eval(marker_type_column), "marker_type")

      }
      if(!resource_column %in% colnames(marker_table)){
        warning("\"resource_column\" not found. By default a unique resource will be associated to possible recurrence markers.")
        marker_table[,resource:="custom_set"]
      } else{
        setnames(marker_table, eval(resource_column), "resource")

      }


      marker_table[,NCIT_cell_definition:=NA]
      marker_table[,NCIT_ID:=NA]
      marker_table[,original_celltype:=NCIT_celltype]
      marker_table[,log2FC:=NA]
      marker_table[,p.value:=NA]
      marker_table[,adjusted_p.value:=NA]
      marker_table[,pct1:=NA]

      disease_accordion_marker$SPs_disease_tissue_specific<-NULL
      disease_accordion_marker$SPs_disease_specific<-NULL
      disease_accordion_marker$SPs_global<-NULL
      disease_accordion_marker$ECs_disease_tissue_specific<-NULL
      disease_accordion_marker$ECs_disease_specific<-NULL
      disease_accordion_marker$ECs_global<-NULL
      marker_table<-unique(marker_table)

      #Add Uberon_ID
      data("uberon_onto", package = "cellmarkeraccordion",envir = environment())

      ontology_tissue<-as.data.frame(uberon_onto[["name"]])
      colnames(ontology_tissue)<-"Uberon_tissue"
      ontology_tissue$Uberon_ID<-rownames(ontology_tissue)
      ontology_tissue$tissue_definition<-uberon_onto[["def"]]
      ontology_tissue<-as.data.table(ontology_tissue)
      marker_table<-merge(marker_table, ontology_tissue, by="Uberon_tissue", all.x=TRUE)

      marker_table<-marker_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1")]
      marker_table<-unique(marker_table)

      marker_table[,cell_definition:=str_remove_all(cell_definition, '"')]
      marker_table[,cell_definition:=str_remove_all(cell_definition, "\\$")]
      marker_table[,cell_definition:=str_remove_all(cell_definition, r"(\\)")]
      marker_table[,cell_definition:=tstrsplit(cell_definition, "[", fixed = TRUE, keep = 1)]
      marker_table[,gene_description:=tstrsplit(gene_description,"[",fixed=TRUE,keep=1)]
      marker_table[,tissue_definition:=str_remove_all(tissue_definition, '"')]
      marker_table[,tissue_definition:=str_remove_all(tissue_definition, "\\$")]
      marker_table[,tissue_definition:=str_remove_all(tissue_definition, r"(\\)")]
      marker_table[,tissue_definition:=tstrsplit(tissue_definition, "[", fixed = TRUE, keep = 1)]

      marker_database<-rbind(disease_accordion_marker, marker_table)


      #Compute score

      #calculate EC score considering disease and tissues
      table<-unique(marker_database[,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
      ECs<-ddply(table,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs_disease_tissue_specific")
      marker_database<-merge(marker_database,ECs,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type"), all.x = TRUE)

      #calculate EC score considering disease
      table<-unique(marker_database[,c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
      ECs<-ddply(table,.(species,DO_diseasetype,DO_ID,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs_disease_specific")
      marker_database<-merge(marker_database,ECs,by=c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker","marker_type"), all.x = TRUE)

      #calculate EC score global
      table<-unique(marker_database[,c("species","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
      ECs<-ddply(table,.(species,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs_global")
      marker_database<-merge(marker_database,ECs,by=c("species","NCIT_celltype","NCIT_ID","marker","marker_type"), all.x = TRUE)


      #Compute specificity
      table<-unique(marker_database[,c("species","NCIT_celltype","NCIT_ID","marker","marker_type")])
      mark_spec<-ddply(table,.(species, marker, marker_type),nrow)
      colnames(mark_spec)<-c("species","marker","marker_type","SPs_global")
      marker_database<-merge(marker_database,mark_spec,by=c("species","marker","marker_type"),all.x = TRUE)
      marker_database[,SPs_global:=format(round(1/SPs_global,2), nsmall=2)]

      #Compute specificity: disease specific
      table<-unique(marker_database[,c("DO_diseasetype","DO_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
      mark_spec<-ddply(table,.(species,DO_diseasetype,DO_ID, marker,marker_type),nrow)
      colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","marker","marker_type", "SPs_disease_specific")
      marker_database<-merge(marker_database,mark_spec,by=c("species","DO_diseasetype","DO_ID","marker","marker_type"),all.x = TRUE)
      marker_database[,SPs_disease_specific:=format(round(1/SPs_disease_specific,2), nsmall=2)]

      #Compute specificity: disease and Tissue specific
      table<-unique(marker_database[,c("DO_diseasetype","DO_ID","Uberon_tissue", "Uberon_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
      mark_spec<-ddply(table,.(species,DO_diseasetype,DO_ID, Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
      colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type", "SPs_disease_tissue_specific")
      marker_database<-merge(marker_database,mark_spec,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
      marker_database[,SPs_disease_tissue_specific:=format(round(1/SPs_disease_tissue_specific,2), nsmall=2)]


      marker_database<-marker_database[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1","ECs_global","ECs_disease_specific","ECs_disease_tissue_specific","SPs_global","SPs_disease_specific","SPs_disease_specific")]

  }

   return(marker_database)
  }



