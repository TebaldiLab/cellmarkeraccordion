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
```

# Annotate and interprete single-cell populations with the built-in Cell Marker Accordion database
*cellmarkeraccordion* allows to automatically identifies hematopoietic populations in single-cell dataset by exploiting the built-in collections of marker genes. 
In addition, the package provides an easy interpretation of the results by reporting the for each group of cells the top marker genes which mostly impacted the annotation, together with the top cell types and their relationship based on the cell ontology tree (thanks to the include_detailed_annotation_info and plot parameters). 
To perform cell types identification and obtain detailed annotation information simply run:
```bash
data<-accordion(data)
```



