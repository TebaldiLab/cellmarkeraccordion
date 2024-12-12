
[![R-CMD-check](https://github.com/TebaldiLab/cellmarkeraccordion/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TebaldiLab/cellmarkeraccordion/actions/workflows/R-CMD-check.yaml)

![Logo](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/f71d49b1-72c9-4c45-99d8-e682248154ab)
# cellmarkeraccordion
### R package for automatically annotating and interpreting single-cell populations.
## Installation 
To install the cellmarkeraccordion package directly from GitHub the devtools package is required. If not already installed on your system, run
```bash
install.packages("devtools")
```
Otherwise, load devtools and install the cellmarkeraccordion by
```bash
library(devtools)
install_github("TebaldiLab/cellmarkeraccordion", dependencies = TRUE)
```
## Loading
To load the cellmarkeraccordion run
```bash
library(cellmarkeraccordion)
library(Seurat)
library(data.table)
```
## Input data
All the functions of the <strong>cellmarkeraccordion</strong> accept as input either a Seurat object or a raw or normalized count matrix. 
As an example we used a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. 
Load the raw counts and create a Seurat object
```bash
load(system.file("extdata", "counts.rda", package = "cellmarkeraccordion")) #counts data
# Create Seurat Object
data <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
```
Process and cluster the data
```bash
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.8)
data <- RunUMAP(data, dims = 1:10)
```
## Annotate and interprete single-cell populations with the built-in Cell Marker Accordion database
<strong>cellmarkeraccordion</strong> allows to automatically identifies cell populations in multiple tissue in single-cell dataset by running function ``` accordion ```. 
It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns. The cell types annotation is performed by exploiting the built-in Cell Marker Accordion database of marker genes. In addition, this function provides an easy interpretation of the results by reporting the for each group of cells the top marker genes which mostly impacted the annotation, together with the top cell types and their relationship based on the cell ontology tree (thanks to the *include_detailed_annotation_info* and *plot* parameters). 

Run list_tissues() function to explore which tissues are available in the Cell Marker Accordion. 
```bash  
available_tissue<-list_tissues(species = "Human")
available_tissue[1:20]
```

To perform cell types identification by cluster and obtain detailed annotation information simply run:
```bash  
# Input: Seurat object
# Output: Seurat object with annotation results 
data <- accordion(data, assay ="RNA", species ="Human", tissue="blood", annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

```bash
DimPlot(data, group.by = "accordion_per_cluster")
```
![Annotation_pbmc](https://github.com/user-attachments/assets/9fea9931-a0b7-47bf-8748-b6ef614acfc5)

Or you can use raw counts matrix and specify cluster's id for each cell.
```bash
# Input: raw counts and clusters id  
raw_counts <- GetAssayData(data, assay="RNA", slot='counts')
clusters<- data.table(cell = rownames(data@meta.data), cluster = data@meta.data$seurat_clusters)
# Output: list with annotation results 
output <- accordion(counts, assay ="RNA", species ="Human", tissue="blood", cluster_info = clusters, annotation_resolution= "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

## Improve the biological interpretation of the results
The Cell Marker Accordion has been developed to improve the biological interpretation of the results, by returning a dot plot listing the top N cell types for each cluster. The dot size is proportional to the impact score, and the winning annotation is highlighted.

```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["global"]]
```

![Top_ct_global](https://github.com/user-attachments/assets/4ecdcf91-de64-4edc-8e04-20e61de99379)


Further insight is provided by looking at the percentage of cells assigned to the top scoring cell types, and their similarity based on the Cell Ontology hierarchy
```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["3_B cell"]]
```
![Top_ct_Bcell](https://github.com/user-attachments/assets/61301129-3d0c-43b3-9fc2-bc81c570b5a5)


Finally the top N marker genes contributing to the annotation of each specific cell type can be explored 
```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cluster_plot"]][["global"]]
```
![Top_marker_global](https://github.com/user-attachments/assets/aedc38ec-9cdc-4e03-a9c1-c682d0ecd582)


## Cell type or pathways identification with custom genes sets
<strong>cellmarkeraccordion</strong> performs automatic identification of cell populations based on a custom input set of marker genes by running function ```accordion_custom ```. It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns and a table of marker genes associated to cell types or  to pathways. The marker table should contains at least two columns, the *category_column*,  which specifies cell types or categories, and the *marker_column*, which specifies the corresponding markers on each row. Columns indicating the marker type (either positive or negative), and the marker weight can be optionally included. We used a published human retinal dataset (Lu et al., Dev Cell. 2020) and we included a table of well-known markers associated to retinal cell types.

```bash
load(system.file("extdata", "retina_markers.rda", package = "cellmarkeraccordion"))
head(retina_markers)
```

| cell_type  | marker |                                         
| ------------- | ------------- | 
| Cones | ARR3 | 
| Retinal ganglion cells	 | ATOH7 | 
| Retinal ganglion cells	 | POU4F1 |
| Rods | C11orf96 | 
| Bipolar cells	 | CA10 |
| Bipolar cells	 | CADPS |


```bash
load(system.file("extdata", "retinal_data.rda", package = "cellmarkeraccordion")) 
```

To perform the annotation with the custom marker genes run:
```bash
retinal_data <-accordion_custom(retinal_data, annotation_resolution = "cluster", marker_table  = retina_markers, category_column = "cell_type", marker_column = "marker", min_n_marker = 2, plot=T, annotation_name = "cell_type_retina")

DimPlot(retinal_data, group.by = "cell_type_retina_per_cluster", reduction = "umap.integrated", label=T) + NoLegend()
```

![Retina_ct_accordion_custom](https://github.com/user-attachments/assets/688652ba-7066-45ab-bc25-1e697d3c6bee) 

You can also exploit the ```accordion_custom``` function to explore the expression of group of genes associated to a specific pathway. As an example: 
```bash
load(system.file("extdata", "marker_table_pathway.rda", package = "cellmarkeraccordion")) 
```
| pathway  | genes |                                         
| ------------- | ------------- | 
| apoptosis | AKT3 | 
| apoptosis | IRAK3 | 
| apoptosis | CHP1 |
| apoptosis | CHUK | 
| apoptosis | CSF2RB |
| apoptosis | DFFA |
| apoptosis | DFFB |
| apoptosis | ENDOG |
| apoptosis | AKT1 |
| apoptosis | AKT2 |

And simply run: 
```bash
retinal_data<-accordion_custom(retinal_data, marker_table_pathway, category_column= "pathway", marker_column ="genes", annotation_resolution = "cell",annotation_name = "apoptosis_signature")

FeaturePlot(retinal_data, features = "apoptosis_signature_per_cell_score",  max.cutoff = "q90")
```

![Retina_apoptosis](https://github.com/user-attachments/assets/41a887ce-e98e-48bf-98c5-2b71956e87bd)

## Automatically identify and interpreting cell cycle state of single-cell populations
<strong>cellmarkeraccordion</strong> provides the ```accordion_cellcycle``` function to automatically assing cell cycle state to cell populations. This function exploits the built-in collection of
marker genes associated to each cell cycle phase (G0, G1, G2M, S). It takes in input either a Seurat object or a raw or normalized count matrix. 
To perform cell cycle identification run: 
```bash
data<-accordion_cellcycle(data)
```

![CellCycle](https://github.com/user-attachments/assets/9a1f7e1d-5a48-4fbc-ade1-a29a6d7c6b2c)

## Annotate and interprete aberrant single-cell populations with the built-in Cell Marker Accordion disease database
<strong>cellmarkeraccordion</strong> includes the ```accordion_disease``` function which allows the identification of aberrant populations exploiting the built-in Accordion gene marker disease database. 
This function requires in input either a Seurat object or a raw or normalized count matrix. It is possible to specific both disease and critical cells to identify, thanks to *disease* and *cell_types* parameters. We analyzed a published scRNA-seq dataset of CD34+ bone marrow cells from 5 healthy controls and 14 acute myeloid leukemia patients.

```bash
load(system.file("extdata", "bone_marrow_data.rda", package = "cellmarkeraccordion"))
```


To identify for example "leukemia stem cell" in "acute myeloid leukemia" samples run: 
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, assay = "RNA", species="Human",tissue="bone marrow", disease= "acute myeloid leukemia", NCIT_celltypes = "Leukemic Hematopoietic Stem Cell",annotation_resolution = "cell", max_n_marker = 30, log2FC_threshold = 1, plot=F, annotation_name = "LHSC")

FeaturePlot(bone_marrow_data, features = "LHSC_per_cell_score", min.cutoff = "q15",max.cutoff = "q95", split.by="condition", cols = c("gray","red"), order = T) 
```

![LSCH_score](https://github.com/user-attachments/assets/1b11885b-5345-4c17-9027-9559665c9952)

To identify instead neoplastic monocyte run:
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, assay = "RNA", species="Human",tissue="bone marrow", disease= "acute myeloid leukemia", NCIT_celltypes = "Neoplastic Monocyte",annotation_resolution = "cell",max_n_marker = 30, log2FC_threshold = 1, plot=F, annotation_name = "Neoplastic_monocyte")

FeaturePlot(bone_marrow_data, features = "Neoplastic_monocyte_per_cell_score", min.cutoff = "q10",max.cutoff = "q90", split.by="condition", cols = c("gray","blue"), order = T) 
```
![Neo_mono_score](https://github.com/user-attachments/assets/f6b54ca9-f7a3-45c7-ad3f-00d19ce6b77d)





