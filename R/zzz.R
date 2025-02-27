utils::globalVariables(c("include_detailed_annotation_info_helper","include_detailed_annotation_info_helper_custom","data",
                         "cell"       ,         "..col"       ,        "log2FC" ,             "marker"            ,  "ECs_sum" ,       "ECs_global"    , "Uberon_tissue"  ,     "uberon_onto"        ,
                         "Uberon_ID"   ,        "CL_celltype",         "CL_ID"   ,            "marker_type"      ,   "group"         ,      "ECs"          ,  "marker_type_new",     "ECs_new"       ,
                         "ECs_scaled",     "SPs_scaled",  "SPs",         "combined_score"  ,    "score"          ,     "expr_scaled"      ,   "diff_score"      ,    "annotation_per_cell",
                         "ncell_tot_cluster",   "percentage"   ,       "seurat_clusters",     "weight"         ,     "cell_type"       ,    "..col_vec"       ,    "DO_ID"            ,   "NCIT_ID",
                         "cell_onto"         ,  "impact_score",        "color_combo"     ,    "ECs_range",      "weight_range"     ,   "win_ct"         ,     "win_ct_border"     ,  "freq",
                         "condition"          , "node1.name",          "node2.name"       ,   "x",                   "y"                 ,  "NCIT_celltype" ,      "accordion_marker" ,
                         "ncell_celltype_cluster","perc_celltype_cluster","quantile_score_cluster","annotation_per_cluster","quantile_score_marker","quantile_diff_score",
                         "quantile_combined_score","perc_celltype_cluster","annotation_per_cluster",
                         "weight_scaled","dt_top_marker","ECs_NCIT_global","DO_diseasetype","SPs_ratio","SPs_range","SPs_positive","SPs_negative",
                         "tot_cell_ct_cluster","disease_accordion_marker", "cell_onto","uberon_onto", "species", "original_tissue", "resource","original_celltype", "p.value","adjusted_p.value", "pct1","gene_description",
                         "cell_definition","tissue_definition", "SPs_global","SPs_tissue_specific","original_diseasetype","DO_definition","NCIT_cell_definition","SPs_disease_specific","SPs_disease_tissue_specific","ECs_reg","SPs_reg",""




                         ))
.onAttach <- function(libname, pkgname) {
  art <- "

                          ...............
    .....                 ...............
  ........                ...............
   ....                   .             .
    ..                    .             .                  ......
    ...                  .              .                 .......
     ...                 .             .                 ...   ..
      ...                .             .                ...
    ......               .        ......         ....  ...
  .........        .......       .......       ..........
 ...........      .......         ....         .........
  .........                                    ........
   .......                                       ....

"
  packageStartupMessage(art)
  packageStartupMessage("")
}
