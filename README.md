![Logo](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/f71d49b1-72c9-4c45-99d8-e682248154ab)
# cellmarkeraccordion
 R package for automatically annotating and interpreting single-cell populations.
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
# Input data
All the functions of the <strong>cellmarkeraccordion</strong> accept as input either a Seurat object or a raw or normalized count matrix. 
As an example we used a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. 
Load the raw counts and create a Seurat object
```bash
data(counts) #raw counts
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
# Annotate and interprete single-cell populations with the built-in Cell Marker Accordion database
<strong>cellmarkeraccordion</strong> allows to automatically identifies hematopoietic populations in single-cell dataset by running function ``` accordion ```. 
It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns. The cell types annotation is performed by exploiting the built-in Cell Marker Accordion database of marker genes. In addition, this function provides an easy interpretation of the results by reporting the for each group of cells the top marker genes which mostly impacted the annotation, together with the top cell types and their relationship based on the cell ontology tree (thanks to the *include_detailed_annotation_info* and *plot* parameters). 
To perform cell types identification by cluster and obtain detailed annotation information simply run:
```bash  
# Input: Seurat object
# Output: Seurat object with annotation results 
data <- accordion(data, annotation_resolution = "cluster", include_detailed_annotation_info = TRUE, plot = TRUE)
```
Or 
```bash
# Input: raw counts and clusters id  
raw_counts <- data@assays[["RNA"]]@counts
clusters<- data.table(cell = rownames(data@meta.data), cluster = data@meta.data$seurat_clusters)
# Output: list with annotation results 
output <- accordion(counts, cluster_info = clusters, annotation_resolution= "cluster", include_detailed_annotation_info = TRUE, plot = TRUE)
```

# Cell type or pathways identification with custom genes sets
<strong>cellmarkeraccordion</strong> performs automatic identification of cell populations based on a custom input set of marker genes by running function ```accordion_custom ```.It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns and a table of marker genes associated to cell types or  to pathways. The marker table should contains at least two columns, the *category_column*,  which specifies cell types or categories, and the *marker_column*, which specifies the corresponding markers on each row. Columns indicating the marker type (either positive or negative), and the marker weight can be optionally included. An example of marker table:

| cell_type  | marker |
| ------------- | ------------- |
| Naive CD4+ T | IL7R |
| Naive CD4+ T | CCR7 |
| CD14+ Mono | CD14 |
| CD14+ Mono | LYZ |
| Memory CD4+ | IL7R |
| Memory CD4+ | S100A4 |
| B | MS4A1 |
| CD8+ T | CD8A |
| FCGR3A+ Mono | MS4A7 |
| FCGR3A+ Mono | FCGR3A |
| NK | NKG7 |
| NK | GNLY |
| DC | FCER1A |
| DC | CST3 |
| Platelet | PPBP |

To perform custom annotation run:
```bash
data<-accordion_custom(data, marker_table, category_column= "cell_type", marker_column ="marker", marker_type_column = "marker_type", weight_column = "weight")
```

# Automatically identify and interpreting cell cycle state of single-cell populations
<strong>cellmarkeraccordion</strong> provides the ```accordion_cellcycle``` function to automatically assing cell cycle state to cell populations. This function exploits the built-in collection of
marker genes associated to each cell cycle phase (G0, G1, G2M, S). It takes in input either a Seurat object or a raw or normalized count matrix. 
To perform cell cycle identification run: 
```bash
data<-accordion_cellcycle(data)
```
# Annotate and interprete aberrant single-cell populations with the built-in Cell Marker Accordion disease database
<strong>cellmarkeraccordion</strong> includes the ```accordion_disease``` function which allows the identification of aberrant populations exploiting the built-in Accordion gene marker disease database. 
This function requires in input either a Seurat object or a raw or normalized count matrix. It is possible to specific both disease and critical cells to identify, thanks to *disease* and *cell_types* parameters.
To identify for example "leukemia stem cell" in "acute myeloid leukemia" samples run: 
```bash
data<-accordion_disease(data, disease="acute myeloid leukemia", cell_types ="leukemia stem cell")
```






