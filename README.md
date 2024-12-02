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
data <- accordion(data,tissue="blood", annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

```bash
DimPlot(data, group.by = "accordion_per_cluster")
```
![Annotation_example](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/673e5368-0014-444d-916c-873d0b522b7e)

Or you can use raw counts matrix and specify cluster's id for each cell.
```bash
# Input: raw counts and clusters id  
raw_counts <- GetAssayData(data, assay="RNA", slot='counts')
clusters<- data.table(cell = rownames(data@meta.data), cluster = data@meta.data$seurat_clusters)
# Output: list with annotation results 
output <- accordion(counts, tissue="blood", cluster_info = clusters, annotation_resolution= "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

## Improve the biological interpretation of the results
The Cell Marker Accordion has been developed to improve the biological interpretation of the results, by returning a dot plot listing the top N cell types for each cluster. The dot size is proportional to the impact score, and the winning annotation is highlighted.

```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["global"]]
```

![Pbmc_top_cell_type](https://github.com/user-attachments/assets/29c20803-e134-44ab-ab0f-083256c266c1)

Further insight is provided by looking at the percentage of cells assigned to the top scoring cell types, and their similarity based on the Cell Ontology hierarchy
```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["2_naive B cell"]]
```
![Pbmc_rank_naive](https://github.com/user-attachments/assets/1372c184-f3e1-43eb-908b-1419e7c21698)

Finally the top N marker genes contributing to the annotation of each specific cell type can be explored 
```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cluster_plot"]][["global"]]
```
![Pbmc_top_markers](https://github.com/user-attachments/assets/282337b9-b897-4880-a422-ebfec7ecc7a5)

## Cell type or pathways identification with custom genes sets
<strong>cellmarkeraccordion</strong> performs automatic identification of cell populations based on a custom input set of marker genes by running function ```accordion_custom ```. It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns and a table of marker genes associated to cell types or  to pathways. The marker table should contains at least two columns, the *category_column*,  which specifies cell types or categories, and the *marker_column*, which specifies the corresponding markers on each row. Columns indicating the marker type (either positive or negative), and the marker weight can be optionally included. We used a published human retinal dataset (Lu et al., Dev Cell. 2020) and we included a table of well-known markers associated to retinal cell types.

```bash
read_excel("retina_markers.xlsx")
head(retina_marker)
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
data(retinal_data)
```

To perform the annotation with the custom marker genes run:
```bash
retinal_data <-accordion_custom(retinal_data, annotation_resolution = "cluster", marker_table  = retina_markers, category_column = "cell_type", marker_column = "marker", min_n_marker = 2, plot=T, annotation_name = "cell_type_retina")

DimPlot(retinal_data, group.by = "cell_type_retina_per_cluster", reduction = "umap.integrated", label=T) + NoLegend()
```

![Retina_ct_accordion_custom](https://github.com/user-attachments/assets/688652ba-7066-45ab-bc25-1e697d3c6bee) 

You can also exploit the ```accordion_custom``` function to explore the expression of group of genes associated to a specific pathway. As an example: 
```bash
data(marker_table_pathway)
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
load(bone_marrow_data)
```


To identify for example "leukemia stem cell" in "acute myeloid leukemia" samples run: 
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, disease= "acute myeloid leukemia", NCIT_celltype = "Leukemia Hematopoietic Stem Cell",combined_score_quantile_threshold = 0.75, annotation_resolution = "cell", plot=F, annotation_name = "LSC")

FeaturePlot(bone_marrow_data, features = "LSC_per_cell_score", min.cutoff = "q10",max.cutoff = "q75", cols = c("gray","red"), order = T) 
```

![LSC](https://github.com/user-attachments/assets/11ed18ad-fdbb-4e6e-b7d0-ce68346ddad6)


