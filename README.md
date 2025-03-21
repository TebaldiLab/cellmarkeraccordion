
[![R-CMD-check](https://github.com/TebaldiLab/cellmarkeraccordion/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TebaldiLab/cellmarkeraccordion/actions/workflows/R-CMD-check.yaml)

![Logo](https://github.com/TebaldiLab/cellmarkeraccordion/assets/68125242/f71d49b1-72c9-4c45-99d8-e682248154ab)
# cellmarkeraccordion

### R package for automated annotation and interpretation of single-cell and spatial omics data.
## Overview
Single-cell technologies offer a unique opportunity to explore cellular heterogeneity in health and
disease. However, reliable identification of cell types and states represents a bottleneck. Available
databases and analysis tools employ dissimilar markers, leading to inconsistent annotations and
poor interpretability. Furthermore, current tools focus mostly on physiological cell types, limiting their
applicability to disease. <br>
We developed Cell Marker Accordion, a user-friendly platform that includes both an R package and a [Shiny app](https://rdds.it/CellMarkerAccordion/). This tool provides automated annotation and biological interpretation of single-cell populations using consistency-weighted markers.
Cell Marker Accordion enhances annotation accuracy in single-cell and spatial omics datasets across various human and murine tissues. Additionally, it can identify disease-critical cells and pathological processes, helping to extract potential biomarkers across diverse disease contexts. <br>

To explore, filter and download the Cell Marker database we recommend using our [Shiny app](https://rdds.it/CellMarkerAccordion/).

## Citing the cellmarkeraccordion package
Please cite the following article when using the cellmarkeraccordion package:

<strong>Cell Marker Accordion: interpretable single-cell and spatial omics annotation in health and disease</strong>

Emma Busarello, Giulia Biancon, Ilaria Cimignolo, Fabio Lauria, Zuhairia Ibnat, Christian Ramirez, Gabriele Tomè, Marianna Ciuffreda, Giorgia Bucciarelli, Alessandro Pilli, Stefano Maria Marino, Vittorio Bontempi, Kristin R. Aass, Jennifer VanOudenhove, Maria Caterina Mione, Therese Standal, Paolo Macchi, Gabriella Viero, Stephanie Halene, Toma Tebaldi

bioRxiv 2024.03.08.584053; doi: https://doi.org/10.1101/2024.03.08.584053 

## Installation 
To install the `cellmarkeraccordion` package directly from GitHub the `devtools` package is required. If not already installed on your system, run:
```bash
install.packages("devtools")
```
Next, load `devtools` and install `cellmarkeraccordion` using:
```bash
library(devtools)
install_github("TebaldiLab/cellmarkeraccordion", dependencies = TRUE)
```
## Loading the package
Once installed, load `cellmarkeraccordion` along with `Seurat` and `data.table` packages required for this tutorial:
```bash
library(cellmarkeraccordion)
library(Seurat)
library(data.table)
library(ggplot2)
```
## Access and download the Accordion database
To access the *healthy* Accordion database run:
```bash
data(accordion_marker)
```
To access the *disease* Accordion database run:
```bash
data(disease_accordion_marker)
```
To download the Accordion database as an Excel file click the Download button in the Cell Marker Accordion Shiny app available at: https://rdds.it/CellMarkerAccordion/.<br>

<img src= https://github.com/user-attachments/assets/ea25a808-68ed-406f-a655-16f1ebbe00ac style="width:30%; height:30%;"> <br />


Alternatively, download the "AccordionDB.xlsb" file from the Shiny app’s GitHub repository: [Download AccordionDB from Shiny app repo](https://github.com/TebaldiLab/shiny_cellmarkeraccordion/blob/main/AccordionDB.xlsb).

## Input data
All the functions of the <strong>cellmarkeraccordion</strong> accept as input either a Seurat object or a raw or normalized count matrix. 
As an example we used a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. 
Load the raw counts and create a Seurat object:
```bash
load(system.file("extdata", "counts.rda", package = "cellmarkeraccordion")) #counts data
# Create Seurat Object
data <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
```
Process and cluster the data:
```bash
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.4)
data <- RunUMAP(data, dims = 1:10)
```
## Annotate and interpret single-cell and spatial omics data with the built-in Cell Marker Accordion database
<strong>cellmarkeraccordion</strong> allows to automatically identifies cell populations in multiple tissues in single-cell dataset by running function ``` accordion ```. 
It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns. The cell types annotation is performed by exploiting the built-in Cell Marker Accordion database of marker genes. In addition, this function provides an easy interpretation of the results by reporting for each group of cells the top marker genes which mostly impacted the annotation, together with the top cell types and their relationship based on the cell ontology tree (thanks to the *include_detailed_annotation_info* and *plot* parameters). 

Run ``` list_tissues()```   function to explore which tissues are available in the Cell Marker Accordion:
```bash  
available_tissue<-list_tissues(species = "Human")
available_tissue[1:20]
```

To perform cell types identification by cluster (*annotation_resolution = "cluster"* by default) and obtain detailed annotation information simply run:
```bash  
# Input: Seurat object
# Output: Seurat object with annotation results 
data <- accordion(data, assay ="RNA", species ="Human", tissue="blood", annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

```bash
DimPlot(data, group.by = "accordion_per_cluster")
```
![Annotation_pbmc](https://github.com/user-attachments/assets/9fea9931-a0b7-47bf-8748-b6ef614acfc5)

Or you can use raw counts matrix and specify cluster's id for each cell:
```bash
# Input: raw counts and clusters id  
raw_counts <- GetAssayData(data, assay="RNA", slot='counts')
clusters<- data.table(cell = rownames(data@meta.data), cluster = data@meta.data$seurat_clusters)
# Output: list with annotation results 
output <- accordion(counts, assay ="RNA", species ="Human", tissue="blood", cluster_info = clusters, annotation_resolution= "cluster", max_n_marker = 30, include_detailed_annotation_info = TRUE, plot = TRUE)
```

## Improve the biological interpretation of the results
The Cell Marker Accordion has been developed to improve the biological interpretation of the results, by returning a dot plot listing the top N cell types for each cluster. The dot size is proportional to the impact score, and the winning annotation is highlighted. Detailed annotation information is stored in the *misc* slot of the Seurat object:
```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["global"]]
```

![Top_ct_global](https://github.com/user-attachments/assets/4ecdcf91-de64-4edc-8e04-20e61de99379)


Further insight is provided by looking at the percentage of cells assigned to the top scoring cell types, and their similarity based on the Cell Ontology hierarchy:
```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_celltypes_plot"]][["3_B cell"]]
```
![Top_ct_Bcell](https://github.com/user-attachments/assets/61301129-3d0c-43b3-9fc2-bc81c570b5a5)


Finally the top N marker genes contributing to the annotation of each specific cell type can be explored:
```bash
data@misc[["accordion"]][["cluster_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cluster_plot"]][["global"]]
```
![Top_marker_global](https://github.com/user-attachments/assets/aedc38ec-9cdc-4e03-a9c1-c682d0ecd582)


## Cell type or pathways identification with custom gene sets
<strong>cellmarkeraccordion</strong> performs automatic identification of cell populations based on a custom input set of marker genes by running the ```accordion_custom ``` function. It requires in input only a Seurat object or a raw or normalized count matrix with genes on rows and cells on columns and a table of marker genes associated to cell types or  to pathways. The marker table should contains at least two columns, the *category_column*,  which specifies cell types or categories, and the *marker_column*, which specifies the corresponding markers on each row. Columns indicating the marker type (either positive or negative), and the marker weight can be optionally included. We used a published human retinal dataset (Lu et al., Dev Cell. 2020) and we included a table of well-known markers associated to retinal cell types.

Load custom markers:
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

Load the already processed Seurat object:
```bash
load(system.file("extdata", "retinal_data.rda", package = "cellmarkeraccordion")) 
```

To perform the annotation by cluster with the custom marker genes run:
```bash
retinal_data <-accordion_custom(retinal_data, annotation_resolution = "cluster", marker_table  = retina_markers, category_column = "cell_type", marker_column = "marker", min_n_marker = 2, plot=T, annotation_name = "cell_type_retina")

DimPlot(retinal_data, group.by = "cell_type_retina_per_cluster", reduction = "umap.integrated", label=T) + NoLegend()
```

![Retina_ct_accordion_custom](https://github.com/user-attachments/assets/688652ba-7066-45ab-bc25-1e697d3c6bee) 

You can also exploit the ```accordion_custom``` function to explore the expression of group of genes associated to a specific pathway. As an example: 
```bash
load(system.file("extdata", "marker_table_pathway.rda", package = "cellmarkeraccordion"))
head(marker_table_pathway, 10)
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

And simply run the *accordion_custom* function by setting *annotation_resolution = cell* : 
```bash
retinal_data<-accordion_custom(retinal_data, marker_table_pathway, category_column= "pathway", marker_column ="genes", annotation_resolution = "cell",annotation_name = "apoptosis_signature")

FeaturePlot(retinal_data, features = "apoptosis_signature_per_cell_score", reduction="umap.integrated",order = T, max.cutoff = "q90")
```
![Retina_fp](https://github.com/user-attachments/assets/0f763a48-22b9-46de-9064-7053cdbf4859)


## Identification of pathway-specific top genes across multiple cell types and conditions
Recent studies turned the spotlight on aberrant activation of innate immune pathways as a consequence of response to the pharmacological inhibition of the m6A methyltransferase Mettl3. To explore the impact of the inhibition of Mettl3 on immunity in single-cell datasets, the Cell Marker Accordion can be exploit to compute an “innate immune response” score based on the activation of genes associated with this signature. 

As an example dataset we used a published bone marrow dataset from mice upon pharmacological inhibition of Mettl3 with STM245 (Sturgess et al., Leukemia, 2023). We included a table of genes associated to innate immune response signature.

Load the Seurat object already processed:
```bash
load(system.file("extdata", "mouse_bm_data.rda", package = "cellmarkeraccordion"))
table(mouse_bm_data$condition)
```
| Vehicle | STM2457 |
|-----------|-------|
| 480   | 395   |
    
Load the innate immune response signature table:
```bash
load(system.file("extdata", "in_im_resp_sig.rda", package = "cellmarkeraccordion"))
head(in_im_resp_sig)
```
| Symbol  | terms                   |
|---------|-------------------------|
| Acod1   | innate_immune_response  |
| Actg1   | innate_immune_response  |
| Actr2   | innate_immune_response  |
| Actr3   | innate_immune_response  |
| Adam8   | innate_immune_response  |
| Adam15  | innate_immune_response  |


First, cell types annotation can be performed by running the ```accordion``` function, specyfing *species ="Mouse"* and *tissue="bone marrow"*:
```bash
mouse_bm_data <- accordion(mouse_bm_data, assay ="RNA", species ="Mouse", tissue="bone marrow", annotation_resolution = "cluster", max_n_marker = 30, allow_unknown=F, include_detailed_annotation_info = F, plot = F)
DimPlot(mouse_bm_data, group.by = "accordion_per_cluster")
```
![Mouse_anno](https://github.com/user-attachments/assets/7a862144-7bcc-4129-acf1-b4a7f558d393)

Next, the ```accordion_custom``` function can be used to explore the expression of innate immune response genes following Mettl3 inhibition. To identify the most impactful condition-specific genes in vehicle- and STM245-treated mice respectively, we can specify in the *condition_group_info* parameter the column name in the metadata of the Seurat object that contains the cell condition information:
```bash
mouse_bm_data <-accordion_custom(mouse_bm_data, marker_table = in_im_resp_sig,  category_column= "terms", marker_column ="Symbol",  annotation_resolution = "cell", 
                                     condition_group_info = "condition", annotation_name = "innate_immune_response_condition")

#visualize the top markers associated to the innate immune response, for vehicle- and STM245-treated mice respectively:
mouse_bm_data@misc[["innate_immune_response_condition"]][["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell_plot"]][["innate_immune_response"]]
```
![Top_markers_cond](https://github.com/user-attachments/assets/bf6d0c95-447a-437d-b857-3a465a1bae17)


Moreover, the <strong>cellmarkeraccordion</strong> allows to further identify the top N (5 by default) cell type-condition-specific genes, by specifying in the *condition_group_info* and *celltype_group_info* parameters both the condition and the cell type annotation columns of the metadata:
```bash
mouse_bm_data <-accordion_custom(mouse_bm_data, marker_table = in_im_resp_sig,  category_column= "terms", marker_column ="Symbol",  annotation_resolution = "cell", 
                                     condition_group_info = "condition", celltype_group_info = "accordion_per_cluster", annotation_name = "innate_immune_response_celltype_condition")
head(mouse_bm_data@misc[["innate_immune_response_celltype_condition"]][["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell"]], n = 10)
```
| innate_immune_response | condition | celltype                     | marker  | marker_type | gene_impact_score | weight | SPs |
|------------------------|-----------|-----------------------------|---------|-------------|-------------------|--------|-----|
| innate_immune_response | Vehicle   | dendritic cell              | H2-Ab1  | positive    | 8.443083          | 1      | 1   |
| innate_immune_response | Vehicle   | dendritic cell              | Rnase6  | positive    | 6.746420          | 1      | 1   |
| innate_immune_response | Vehicle   | dendritic cell              | Irf8    | positive    | 5.249773          | 1      | 1   |
| innate_immune_response | Vehicle   | dendritic cell              | Unc93b1 | positive    | 4.886153          | 1      | 1   |
| innate_immune_response | Vehicle   | dendritic cell              | Grn     | positive    | 4.766681          | 1      | 1   |
| innate_immune_response | STM2457   | plasmacytoid dendritic cell | Stat1   | positive    | 6.115051          | 1      | 1   |
| innate_immune_response | STM2457   | plasmacytoid dendritic cell | Isg15   | positive    | 5.217695          | 1      | 1   |
| innate_immune_response | STM2457   | plasmacytoid dendritic cell | Gbp7    | positive    | 5.029485          | 1      | 1   |
| innate_immune_response | STM2457   | plasmacytoid dendritic cell | Isg20   | positive    | 4.258571          | 1      | 1   |
| innate_immune_response | STM2457   | plasmacytoid dendritic cell | Irgm1   | positive    | 4.122629          | 1      | 1   |


We can extract the annotation results from the misc slot and visualize the top 5 genes for common lymphoid progenitor and megakaryocyte populations for vehicle- and STM245-treated mice respectively:
```bash
#extract annotation results for common lymphoid progenitor and megakaryocyte populations
dt <- mouse_bm_data@misc[["innate_immune_response_celltype_condition"]][["cell_resolution"]][["detailed_annotation_info"]][["top_markers_per_celltype_cell"]]
dt_filt<- dt[accordion_per_cluster %in% c("mast cell", "megakaryocyte")]
dt_filt<- dt_filt[order(gene_impact_score_per_celltype_cell)][,marker:=factor(marker, levels = unique(marker))]

#customize the plot with ggplot
bs<-20
ggplot(dt_filt, aes(gene_impact_score_per_celltype_cell, marker)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(x = 0, xend = gene_impact_score_per_celltype_cell, y = marker, yend = marker, color=condition),linewidth = bs/10, show.legend = F) +
  geom_point(aes(color=condition),size=6, alpha= 1, shape = 16) +
  theme_bw(base_size = bs) +
  facet_grid(accordion_per_cluster ~ condition, scale="free_y") + 
  scale_color_manual(values = c("gray50","#8B1A1A"), guide="none") + 
  theme(panel.border = element_blank()) +
  labs(x = "Gene impact score", y = "") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=bs, face = "bold"),
        text = element_text(size = bs),
        axis.text.y = element_text(size = bs)) 
```
![Top_markers_cond_celltype](https://github.com/user-attachments/assets/cb723c32-1b43-4855-925d-1879e2d049ec)



## Automatically identify and interpret cell cycle state of single-cell populations
<strong>cellmarkeraccordion</strong> provides the ```accordion_cellcycle``` function to automatically assign cell cycle state to cell populations. This function exploits the built-in collection of
marker genes associated to each cell cycle phase (G0, G1, G2M, S). It takes in input either a Seurat object or a raw or normalized count matrix. 
To perform cell cycle identification run: 
```bash
mouse_bm_data<-accordion_cellcycle(mouse_bm_data)
DimPlot(mouse_bm_data, group.by="accordion_cell_cycle_per_cell")
```
![Cellcycle](https://github.com/user-attachments/assets/c603f14f-00d1-4bc9-948d-5ae283b561a5)

## Annotate and interpret aberrant single-cell populations with the built-in Cell Marker Accordion disease database
<strong>cellmarkeraccordion</strong> includes the ```accordion_disease``` function which allows the identification of aberrant populations exploiting the built-in Accordion gene marker disease database. 
This function requires in input either a Seurat object or a raw or normalized count matrix. It is possible to specify both disease and critical cells to identify, thanks to *disease* and *cell_types* parameters. We analyzed a published scRNA-seq dataset of CD34+ bone marrow cells from 5 healthy controls and 14 acute myeloid leukemia patients.

Load the already processed Seurat object:
```bash
load(system.file("extdata", "bone_marrow_data.rda", package = "cellmarkeraccordion"))
```


To identify for example "Leukemic Hematopoietic Stem Cell" in "acute myeloid leukemia" patients from "bone marrow" samples run: 
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, assay = "RNA", species="Human",tissue="bone marrow", disease= "acute myeloid leukemia", NCIT_celltypes = "Leukemic Hematopoietic Stem Cell",annotation_resolution = "cell", max_n_marker = 30, log2FC_threshold = 1, plot=F, annotation_name = "LHSC")

FeaturePlot(bone_marrow_data, features = "LHSC_per_cell_score", min.cutoff = "q15",max.cutoff = "q95", split.by="condition", cols = c("gray","red"), order = T) 
```

![LSCH_score](https://github.com/user-attachments/assets/1b11885b-5345-4c17-9027-9559665c9952)

To identify instead "Neoplastic Monocyte" run:
```bash
bone_marrow_data = accordion_disease(bone_marrow_data, assay = "RNA", species="Human",tissue="bone marrow", disease= "acute myeloid leukemia", NCIT_celltypes = "Neoplastic Monocyte",annotation_resolution = "cell",max_n_marker = 30, log2FC_threshold = 1, plot=F, annotation_name = "Neoplastic_monocyte")

FeaturePlot(bone_marrow_data, features = "Neoplastic_monocyte_per_cell_score", min.cutoff = "q10",max.cutoff = "q90", split.by="condition", cols = c("gray","blue"), order = T) 
```
![Neo_mono_score](https://github.com/user-attachments/assets/f6b54ca9-f7a3-45c7-ad3f-00d19ce6b77d)


## Integrate custom set of markers with the Accordion database
The <strong>cellmarkeraccordion</strong> package includes the ```marker_database_integration``` function, which allows users to integrate a custom set of marker genes into the Accordion database—either for healthy or disease conditions.

<strong>Usage</strong>

Set the *database* parameter to either:
- "healthy" → Integrate with the healthy Accordion database
- "disease" → Integrate with the disease Accordion database

<strong>Input Requirements</strong>

The function requires a marker gene table with at least two columns:
- "cell_type" – Specifies the cell type
- "marker" – Lists the marker genes

To ensure proper integration, cell types nomenclature should be standardized:
- Healthy database → Use Cell Ontology
- Disease database → Use NCI Thesaurus

If non-standardized cell types are provided, they will be added as "new" cell types in the database.

<strong>Optional Columns</strong>

Additional columns can be included:
- "species": Specifies the species (default: "Human").
- "tissue": Specifies the related tissue. Standardization with Uberon Ontology is recommended for effective integration. Non-standardized tissues will be added as "new" tissues. If omitted, integration will ignore tissue specificity.
- "marker_type": Defines marker type ("positive" or "negative"; default: "positive").
- "resource": Indicates the data source. If omitted, markers are labeled as "custom_set".
- "disease": Required if database = "disease". Standardization with Disease Ontology is recommended. Non-standardized diseases will be added as "new" diseases. If omitted, disease specificity is ignored.

<strong>Running the Integration</strong>

Load custom set of marker genes:
```bash
load(system.file("extdata", "custom_markers_to_integrate.rda", package = "cellmarkeraccordion"))
head(custom_markers_to_integrate, 10)
```

| species | Uberon_tissue | CL_celltype         | marker  | resource     |
|---------|--------------|---------------------|---------|-------------|
| Mouse   | brain        | glutamatergic neuron | Satb2   | custom_set_1 |
| Mouse   | brain        | glutamatergic neuron | Satb2   | custom_set_2 |
| Mouse   | brain        | glutamatergic neuron | Slc17a6 | custom_set_1 |
| Mouse   | brain        | glutamatergic neuron | Slc17a7 | custom_set_1 |
| Mouse   | brain        | glutamatergic neuron | Slc17a7 | custom_set_2 |
| Mouse   | brain        | pyramidal neuron     | Sv2b    | custom_set_1 |
| Mouse   | brain        | pyramidal neuron     | Calb1   | custom_set_1 |
| Mouse   | brain        | pyramidal neuron     | Pde1a   | custom_set_1 |
| Mouse   | brain        | pyramidal neuron     | Pde1a   | custom_set_2 |
| Mouse   | brain        | pyramidal neuron     | Pde1a   | custom_set_3 |

To integrate the custom table with the healthy Accordion database, use:
```bash
database_integrated<-marker_database_integration(marker_table = custom_markers_to_integrate,
                           database = "healthy",
                           species_column = "species",
                           disease_column = "disease",
                           tissue_column = "Uberon_tissue",
                           celltype_column = "CL_celltype",
                           marker_column = "marker",
                           marker_type_column = "marker_type",
                           resource_column = "resource")
```

To perform automatic cell type annotation using the previously integrated marker database, pass the output table from  ```marker_database_integration```  to the *database* parameter of the ```accordion``` function (or ```accordion_disease``` if the integration has been performed with the disease database of the Cell Marker Accordion).

## Annotate and interpret single-cell and spatial omics data with the integrated marker database
As an example we used an adult mouse brain MERFISH dataset (Zhuang et al., 2024) with a panel of 1122 genes. 

Load already processed Seurat object:
```bash
load(system.file("extdata", "brain_data.rda", package = "cellmarkeraccordion"))
```

First, perform cell type annotation with the Cell Marker Accordion database only and visualize the result:
```bash
brain_data <- accordion(brain_data, assay ="SCT", species ="Mouse", tissue="brain", annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = F, plot = F, allow_unknown = F)
DimPlot(brain_data, group.by="accordion_per_cluster")
```
![Merfish_anno_accordion](https://github.com/user-attachments/assets/d2a9e34d-d63a-43e2-83a5-f8ae0b9cfdb5)

Then, perform cell with the integrated database by setting *database = database_integrated* and compare the result:
```bash
brain_data <- accordion(brain_data, assay ="SCT",database=database_integrated, species ="Mouse", tissue="brain", annotation_resolution = "cluster", max_n_marker = 30, include_detailed_annotation_info = F, plot = F, allow_unknown = F, annotation_name = "integrated_database")
```
We can notice that glutamatergic neuron are now identified:
```bash
DimPlot(brain_data, group.by="integrated_database_per_cluster")
```
![Merfish_anno_integratedDB](https://github.com/user-attachments/assets/902c2a4d-6e14-4db4-885b-58cfb9db9e4d)

In the metadata of the Seurat object are stored the x and y coordinates of each cell. We can then visualize the annotation results on the brain tissue:
```bash
ggplot(brain_data@meta.data, aes(x=x, y=-y, color=integrated_database_per_cluster))+
  geom_point(size=1)+
  theme_classic()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    axis.ticks = element_blank(),  
    axis.line = element_blank()    
  )
```

![Merfish](https://github.com/user-attachments/assets/23d43cc0-aabd-4810-8bed-92e69009f10b)


