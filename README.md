![Logo](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/f71d49b1-72c9-4c45-99d8-e682248154ab)
# cellmarkeraccordion
### R package for the automatic annotation of and interpretation of single-cell populations.
## Installation 
To install the <strong>cellmarkeraccordion</strong> package directly from GitHub, the devtools package is required. If not already installed on your system, run:
```bash
install.packages("devtools")
```
Otherwise, load devtools and install the <strong>cellmarkeraccordion</strong> by:
```bash
library(devtools)
install_github("TebaldiLab/cellmarkeraccordion", dependencies = TRUE)
```
## Loading
To load the <strong>cellmarkeraccordion</strong>, run:
```bash
library(cellmarkeraccordion)
library(Seurat)
library(data.table)
```
## Input data
All the functions of the <strong>cellmarkeraccordion</strong> accept as input either a Seurat object (version 4 or 5 are supported) or a matrix with raw or normalized counts. 
As an introductory example, we used a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. 
First, load the raw counts and create a Seurat object:
```bash
data(counts) #raw counts
# Create Seurat Object
data <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
```
Then, process and cluster the data with a standard Seurat workflow:
```bash
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.8)
data <- RunUMAP(data, dims = 1:10)
```
## Annotate and interpret single-cell populations with the built-in Cell Marker Accordion database
The <strong>cellmarkeraccordion</strong> allows the automatic identification of hematopoietic populations in single-cell datasets by running function ``` accordion ```. 
This function requires in input only a Seurat object (v4 or v5), or a matrix with raw or normalized counts (genes/features as rows and cells as columns). By default, the annotation of cell types is performed by exploiting the built-in Cell Marker Accordion database of marker genes. In addition, this function provides an extensive and easy interpretation of the results by reporting, for each cluster or annotated cell type, a) the top n marker genes which mostly impacted the annotation, b) the top n cell types competing for the final annotation and c) their relationship based on the cell ontology tree (these plots are activated by the *include_detailed_annotation_info* and *plot* parameters). 
To identify the best matching cell types for each cluster and obtain detailed annotation information, simply run:
```bash  
# Input: Seurat object
# Output: Seurat object with annotation results 
data <- accordion(data, annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

```bash
DimPlot(data, group.by = "accordion_per_cluster")
```
![Annotation_example](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/673e5368-0014-444d-916c-873d0b522b7e)

Alternatively, the raw count matrix can be used as input, together with a vector assigning each cell to a cluster.
```bash
# Input: raw counts and clusters id  
raw_counts <- GetAssayData(data, assay="RNA", slot='counts')
clusters<- data.table(cell = rownames(data@meta.data), cluster = data@meta.data$seurat_clusters)
# Output: list with annotation results 
output <- accordion(counts, cluster_info = clusters, annotation_resolution= "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

## Enhanced biological interpretation of results
The <strong>cellmarkeraccordion</strong> has been developed to improve the biological interpretation of results. On key output is a dot plot displaying the top N scoring cell types associated with each cluster. The dot size is proportional to the cell type impact score, and the winning annotation is highlighted.

```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["global"]]
```

![Pbmc_top_cell_type](https://github.com/user-attachments/assets/29c20803-e134-44ab-ab0f-083256c266c1)

Further insight is provided by displaying, for each cluster, the percentages of cells directly assigned to the top-scoring cell types and their similarity based on the Cell Ontology hierarchy
```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["2_naive B cell"]]
```
![Pbmc_rank_naive](https://github.com/user-attachments/assets/1372c184-f3e1-43eb-908b-1419e7c21698)

Finally, the top N marker genes contributing to the annotation of each specific cell type can be inspected 
```bash
seurat_obj@misc[["accordion_pbmc"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cluster_plot"]][["global"]]
```
![Pbmc_top_markers](https://github.com/user-attachments/assets/282337b9-b897-4880-a422-ebfec7ecc7a5)

## Cell type identification or pathway activity quantification using custom gene sets
The <strong>cellmarkeraccordion</strong> enables the automatic identification of cell populations based on custom sets of input marker genes (bring your own markers). The function ```accordion_custom ```. requires as input: a) a Seurat object or a raw or normalized count matrix (with genes on rows and cells on columns) and b) a table of marker genes associated with cell types to pathways. This marker table should contain at least two columns: the *category_column*,  which specifies cell types or categories, and the *marker_column*, which specifies the corresponding markers. Optional columns can be used to specify the marker type (either positive or negative) and the marker weight. An example of this procedure is shown based on a published human retinal dataset (Lu et al., Dev Cell. 2020). The annotation was performed with a table of well-known markers associated with retinal cell types.

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

To perform the annotation with custom marker genes, run:
```bash
retinal_data <-accordion_custom(retinal_data, annotation_resolution = "cluster", marker_table  = retina_markers, category_column = "cell_type", marker_column = "marker", min_n_marker = 2, plot=T, annotation_name = "cell_type_retina")

DimPlot(retinal_data, group.by = "cell_type_retina_per_cluster", reduction = "umap.integrated", label=T) + NoLegend()
```

![Retina_ct_accordion_custom](https://github.com/user-attachments/assets/688652ba-7066-45ab-bc25-1e697d3c6bee) 

The ```accordion_custom``` function can be used to quantify the activity of specific pathways based on the expression of associated gene signatures or modules. As an example: 
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

Cell-specific pathway activity scores can be displayed by running: 
```bash
retinal_data<-accordion_custom(retinal_data, marker_table_pathway, category_column= "pathway", marker_column ="genes", annotation_resolution = "cell",annotation_name = "apoptosis_signature")

FeaturePlot(retinal_data, features = "apoptosis_signature_per_cell_score",  max.cutoff = "q90")
```

![Retina_apoptosis](https://github.com/user-attachments/assets/41a887ce-e98e-48bf-98c5-2b71956e87bd)

## Automatic identification and interpretation of single-cell cycle phases
The <strong>cellmarkeraccordion</strong> provides the ```accordion_cellcycle``` function to automatically assign cell cycle phases to cell populations. This function exploits the built-in collection of marker genes associated with each cell cycle phase (G0, G1, G2M, S). This function takes as input either a Seurat object or a raw or normalized count matrix. A published scRNA-seq dataset of bone marrow of Mettl3 conditional knockout mice (Cheng at al., Cell Rep, 2019).

To perform cell cycle identification, run: 
```bash
data<-accordion_cellcycle(data, species = "Mouse")
```

![CellCycle](https://github.com/user-attachments/assets/9a1f7e1d-5a48-4fbc-ade1-a29a6d7c6b2c)

## Identification of disease-critical single-cell populations with the built-in Cell Marker Accordion disease database
The <strong>cellmarkeraccordion</strong> includes the ```accordion_disease``` function, allowing the identification of aberrant or disease-related cell populations. This function exploits the built-in Accordion gene marker disease database. 
This function requires as input either a Seurat object (v4 or v5) or a matrix with raw or normalized counts. It is possible to specify both the disease and the critical cell type to identify, thanks to the *disease* and *cell_types* parameters. A published scRNA-seq dataset of CD34+ bone marrow cells from 5 healthy controls and 14 acute myeloid leukemia patients is used as an example (Van Galen et al., Cell, 2019).

```bash
load(bone_marrow_data)
```

To identify the cell type "leukemia stem cell" associated with the disease "acute myeloid leukemia", run: 
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, disease= "acute myeloid leukemia", cell_types = "leukemia stem cell",combined_score_quantile_threshold = 0.75, annotation_resolution = "cell", plot=F, annotation_name = "LSC")

FeaturePlot(bone_marrow_data, features = "LSC_per_cell_score", min.cutoff = "q10",max.cutoff = "q75", cols = c("gray","red"), order = T) 
```

![LSC](https://github.com/user-attachments/assets/11ed18ad-fdbb-4e6e-b7d0-ce68346ddad6)


