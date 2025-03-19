library(shinydashboard)
library(shiny)
library(stringr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(igraph)
library(shinyWidgets)
library(ontologyIndex)
library(plyr)
library(ontologyPlot)
library(ontoProc)
library(writexl)
library(DT)
library(dplyr)
library(shinyBS)
library(shinyhelper)
library(shinydashboardPlus)
library(readxl)
library(plyr)
library(rstatix)
library(tidyverse)
library(data.table)
library(shinyjs)
library(rstudioapi)
source('helper_function.R')

# REMEMBER TO CHANGE WHEN  -----
setwd(dirname(getActiveDocumentContext()$path)) # to run the app locally
#setwd("/home/rdds/www/apps/CellMarkerAccordion/") # to run the online version of the app on

#load data
load("data/accordion_complete.rda")
load("data/cell_onto.rda")


onto_plot<-onto_plot2(cell_onto, unique(accordion_complete$celltype_ID))
nodes<-as.data.table(onto_plot@nodes)
nodes<-nodes[,V1:=tstrsplit(nodes$V1,"CL", fixed = TRUE, keep = 1)]
nodes<-nodes[,V1:=tstrsplit(nodes$V1,"DOID", fixed = TRUE, keep = 1)]
nodes<-nodes[, V1:=str_replace_all(V1,"\n"," ")]
nodes<-nodes[, V1:=str_sub(V1,1,nchar(V1)-1)]
onto_plot@nodes<-nodes$V1


ontology_def<-unique(accordion_complete[,c("celltype_ID","cell_definition","celltype")])
ontology_celltype<-unique(accordion_complete[,c("celltype","celltype_ID")])


nmarker_values<-c("ALL")
disease_list<-unique(accordion_complete$DO_diseasetype)
names(disease_list)<-unique(accordion_complete$DO_diseasetype)

tissue_list<-unique(accordion_complete$Uberon_tissue)
names(tissue_list)<-unique(accordion_complete$Uberon_tissue)


celltype_list<-unique(accordion_complete$celltype)
names(celltype_list)<-unique(accordion_complete$celltype)
css <- HTML(
  ".dataTables_scrollBody {
    transform:rotateX(180deg);
}
.dataTables_scrollBody table {
    transform:rotateX(180deg);
}"
)

ui <- dashboardPage(
  dashboardHeader(titleWidth  = 500,title = tagList(
    tags$span(
      class = "logo-mini", "CMA"
    ),
    tags$span(
      class = "logo-lg", "The Cell Marker Accordion"
    )
  )
  ),

    ## Sidebar content
  dashboardSidebar(width=300,
                   sidebarMenu(style = "position: relative; overflow: visible;",
                               menuItem("Homepage", tabName = "dashboard", icon = icon("home")),
                               menuItem("Search by tissue and cell types", tabName = "celltype_h", icon = icon("circle-notch")),
                               menuItem("Search by marker genes", tabName = "marker_h", icon = icon("dna")),
                               menuItem("Custom markers integration", tabName = "integration", icon = icon("gear")),
                               menuItem("Marker enrichment analysis", tabName = "anno", icon = icon("stack-overflow")),

                               # Centered Download Button
                               div(style = "display: flex; justify-content: center; margin-top: 20px;",
                                   downloadButton("downloadAccordionDB", HTML("Download the <br>Cell Marker Accordion database"), style = "background-color: transparent;
                            color: white; border: 1px solid white; padding: 10px 15px;
                            text-align: center; width: 80%; font-size: 16px;")
                               )
                   )
  ),
  dashboardBody(tags$head(
    tags$style("#shiny-modal img { max-width: 100%; }"),
    tags$script(HTML("$('body').addClass('sidebar-mini');")),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style("@import url(https://use.fontawesome.com/releases/v6.1.1/css/all.css);")),
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                HTML("<h3>The Cell Marker Accordion</h3>"),
                HTML("<h2>A Web Tool for Single-Cell and Spatial Omics Annotation</h2>"),
                titlePanel(
                  div(
                    style = "display: flex; align-items: center; gap: 5px; text-align: left;",  # Align left with small gap
                    icon("github", class = "fa-lg", lib = "font-awesome", style = "color: black; font-size: 20px;"),
                    shiny::span(
                      HTML("<p style='margin: 0; font-size: 18px;'><h>Encountering issues? Check out our <a href='https://github.com/TebaldiLab/shiny_cellmarkeraccordion' target='_blank'>GitHub page</a>!</p></h>")
                    )
                  )
                ),
                div(img(src = "Logo.png"), style = "text-align: center;"),
                br(),
                br(),
                tags$style("@import url('https://use.fontawesome.com/releases/v6.1.1/css/all.css');"),
                HTML("<p style='text-align: justify;'>
            <h><strong>The Cell Marker Accordion</strong> is a powerful and user-friendly platform designed to enhance the accuracy and interpretation of normal and aberrant cell populations in single-cell and spatial omics data.
            Our framework includes both a <strong>Shiny app</strong> and an
            <a href='https://github.com/TebaldiLab/cellmarkeraccordion' target='_blank'> <strong>R package</strong></a>.
            <br><br>
            The Cell Marker Accordion web interface allows users to easily explore the integrated built-in database of consistency-weighted markers.
            Specifically, it enables:
        </p></h>")
              ),

              titlePanel(shiny::span((icon("circle-notch",class = "about-icon fa-pull-left", lib = "font-awesome")), p(style="text-align: justify;", HTML("<h>  Search and download lists of marker genes associate with input cell types across different tissues in health and disease. </h>")))),
              titlePanel(shiny::span((icon("dna",class = "about-icon fa-pull-left", lib = "font-awesome", style= "margin-top:-30px;")), p(style="text-align: justify;margin-top:-30px;", HTML("<h>  Search and download lists of cell types associated with input marker genes across different tissues in health and disease. </h>")))),
              titlePanel(shiny::span((icon("gear",class ="about-icon fa-pull-left", lib = "font-awesome", style= "margin-top:-30px;")), p(style="text-align: justify;margin-top:-30px;", HTML("<h> Integrate custom set of marker genes with the Cell Marker Accordion database. </h>")))),
              titlePanel(shiny::span((icon("stack-overflow",class ="about-icon fa-pull-left", lib = "font-awesome",style= "margin-top:-30px;")), p(style="text-align: justify;margin-top:-30px;", HTML("<h> Perform cell type marker enrichment analysis across tissues in health and disease. </h>")))),
              HTML("<p style='text-align: left; margin-left:-30px;'> <h>Additionally, in all sections users can easily: </p></h>"),
              

              titlePanel(shiny::span((icon("sitemap",class = "about-icon fa-pull-left", lib = "font-awesome")), p(style="text-align: justify;", HTML("<h>  Browse hierarchies of cell types following the Cell Ontology structure in order to obtain the desired level of specificity in the markers in both search options. </h>")))),
              titlePanel(shiny::span((icon("arrow-down-short-wide",class ="about-icon fa-pull-left", lib = "font-awesome", style= "margin-top:-30px;")), p(style="text-align: justify; margin-top:-30px;", HTML("<h> Rank and select marker genes by their evidence consistency and specificity scores. </h>"))))),




      # Second tab content: 
      # Search markers according to cell types ----
      tabItem(tabName = "celltype_h",
              tags$h2(
                "Search and download lists of marker genes by cell types across different tissues in health and disease",
                style = "font-weight: bold; text-align: center;"
              ),
              tags$style(css),
              div(fluidRow(column(width=12,wellPanel(id="sidebar",                                                                
                                                    div(
                                                      style = "display: flex; justify-content: center; align-items: center;",
                                                      actionButton(
                                                        'celltypeInfo', 'Info',
                                                        icon = icon("info"),
                                                        style = 'margin-top: 5px; margin-bottom: 5px;'
                                                      )
                                                    ),
                                                    
                                                     checkboxGroupInput("species", "Select species:",
                                                                       choiceNames =
                                                                         list(tags$img(src = "human.png"),tags$img(src = "mouse.png")),
                                                                       choiceValues =
                                                                         list("Human","Mouse"),selected = c("Human"),inline=TRUE),
                                                    br(),
                                                    #pickerInput('disease', 'Hematologic disorder', choices= NULL ,multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE)),
                                                    pickerInput('disease', 'Condition',  choices= NULL, selected = "healthy", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    pickerInput('tissue', 'Tissue',  choices= NULL,selected="blood", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    #checkboxInput("collapse_tissue","",value=FALSE))),
                                                    checkboxInput("tissue_aware","Tissue aware",value=TRUE),
                                                    pickerInput('celltype', 'Cell type', choices= NULL ,multiple=TRUE,selected=c("hematopoietic stem cell (Hs, Mm)", "hematopoietic multipotent progenitor cell (Hs, Mm)","hematopoietic lineage restricted progenitor cell (Hs, Mm)"), options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                                                                        br(),
                                                    pickerInput('descendantsof', 'See subtypes of:', choices= NULL,multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes"),choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), 141))),
                                                    checkboxInput("cellid","Plot celltype_ID",value=FALSE)))),
                           
                           #change style sliderinput
                           tags$style(HTML("
  .irs-single, .irs-bar-edge, .irs-bar {
    background: #990000 !important;
  }
")),

                  #HTML(icon("hand-back-point-up"),"<h4> <strong> Click </strong> on a node to look at cell type description</h4>"),
                  fluidRow(column(12,   downloadButton("save_plot", "Save Ontology plot",style = "padding: 10px 20px; font-size: 18px;")), column(12, radioButtons("file_format", "", choices = c("PNG", "PDF"), inline=TRUE),
                                                                                                   
                 tags$style(HTML("
                #file_format .shiny-options-group label {
             font-size: 12px !important;
           }
           #file_format {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))),  # Format selection),
                 column(12, titlePanel(
                   shiny::span(
                     style = "display: flex; align-items: center; justify-content: center; margin-top: -50px;",
                     icon("hand-pointer", class = "about-icon", lib = "font-awesome"),
                     p(style = "text-align: center; ", HTML("<h><strong>Click</strong> on a node to look at cell type description</h>"))
                   )
                 )),  column(12, sliderInput("zoom", "Zoom:", min = 0.5, max = 10, value = 1, step = 0.2),style = "display: flex; align-items: center;justify-content: center;margin-top: -40px;"), # Zoom slider)               
                    column(12,align="center",div(style='width:100%;overflow-x: scroll;height:100%;overflow-y: scroll;', uiOutput('scalableplot')))),
                  #grVizOutput("plot1"),
                  tags$style(
                    '#test {
    cursor: grap;
    color: black;
    }'),
                  conditionalPanel(id="test",condition= "input.plot1_click",wellPanel(textOutput("celltype_def"))),
                  br(),
                  br(),

                  conditionalPanel(condition= "input.descendantsof != ''",titlePanel(shiny::span(p(style="text-align: justify;", HTML("<h>Table options</h>"),actionButton('help', 'Info',icon= icon("info"), align="left"))))),
                  conditionalPanel(condition= "input.descendantsof == ''",titlePanel(shiny::span(p(style="text-align: justify;", HTML("<h>Table options</h>"),actionButton('help_empty', 'Info',icon= icon("info"), align="left"))))),

                  fluidRow(column(width=12,wellPanel(id="sidebar2",
                                                              
                                                              fluidRow(column(3,sliderInput('EC_score','ECs',min = 1, max = 17, value = 1, step=1)),
                                                                       column(3,sliderInput('specificity','SPs', min = 0, max = 1, value = 0, step=0.1)),
                                                                       column(3,radioButtons("tabletype","Table type",c("Simple","Complete"),selected="Simple")),          
                                                              
                                                              conditionalPanel(condition= "input.descendantsof != ''", column(2,radioButtons("mergeDescendant","Merge subtypes", c("Yes","No"),selected="No"))))))),
                  
                  fluidRow(
                    column(12, 
                           downloadButton("downloadData", "Download", 
                                          style = "padding: 10px 20px; font-size: 18px;")), # Bigger button
                    column(12, 
                           radioButtons("downloadType", "",  # No title
                                        choices = c("CSV" = ".csv",
                                                    "XLSX" = ".xlsx",
                                                    "TSV" = ".tsv"), 
                                        inline = TRUE),
                           tags$style(HTML("
                #downloadType .shiny-options-group label {
             font-size: 12px !important;
           }
           #downloadType {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))) # CSS to make choices smaller
                  ),
                  
                  

                  tags$head(tags$style(HTML('
         #sidebar2 {
            background-color: #ad000019;
            font-size: 16px;
        }'))),
                  # bsTooltip("EC_score", "Evidence consensus score (number of databases)", placement = "top", trigger = "hover",
                  #           options = NULL),
                  br(),
                  br(),

                  DTOutput("table1", width = "100%"),
                  br()
              )),

      # search by markers -----
      tabItem(tabName="marker_h",
              tags$style(css),
              tags$h2(
                "Search and download lists of cell types by marker genes across different tissues in health and disease",
                style = "font-weight: bold; text-align: center;"
              ),
              div(fluidRow(column(width=12,wellPanel(id="sidebar",
                                                    
                                                    div(
                                                      style = "display: flex; justify-content: center; align-items: center;",
                                                      actionButton(
                                                        'markerInfo', 'Info',
                                                        icon = icon("info"),
                                                        style = 'margin-top: 5px; margin-bottom: 5px;'
                                                      )
                                                    ),
                                                    
                                                    checkboxGroupInput("speciesM", "Select species:",
                                                                                    choiceNames =
                                                                                      list(tags$img(src = "human.png"),tags$img(src = "mouse.png")),
                                                                                    choiceValues =
                                                                                      list("Human","Mouse"),selected = c("Human"),inline=TRUE),
                                                    br(),
                                                    textInput("marker", "Insert marker genes", value = "OLIG2", width = NULL, placeholder = NULL),
                                                    fileInput("markerfile", "Upload file with marker genes ",buttonLabel=list(icon("upload")),
                                                              multiple = FALSE),
                                                    pickerInput('diseaseM', 'Condition', selected = "healthy", choices= NULL, multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    pickerInput('tissueM', 'Tissue',  choices= NULL,selected="brain", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    checkboxInput("tissue_awareM","Tissue aware",value=TRUE),

                                                    checkboxInput("cellidM","Plot celltype_ID",value=FALSE)))),
                  
                  
                  

                  
                
                  tags$head(tags$style(HTML('
         #sidebar {
            background-color: #ad000019;
            font-size: 16px;
        }'))),

                  fluidRow(column(12,   downloadButton("save_plotM", "Save Ontology plot",style = "padding: 10px 20px; font-size: 18px;")), column(12, radioButtons("file_formatM", "", choices = c("PNG", "PDF"), inline=TRUE),
                                                                                                                                                  
                                                                                                                                                  tags$style(HTML("
                #file_formatM .shiny-options-group label {
             font-size: 12px !important;
           }
           #file_formatM {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))),
                           # Format selection),
                           column(12, titlePanel(
                             shiny::span(
                               style = "display: flex; align-items: center; justify-content: center; margin-top: -50px;",
                               icon("hand-pointer", class = "about-icon", lib = "font-awesome"),
                               p(style = "text-align: center; ", HTML("<h><strong>Click</strong> on a node to look at cell type description</h>"))
                             )
                           )),  column(12, sliderInput("zoomM", "Zoom:", min = 0.5, max = 10, value = 1, step = 0.2),style = "display: flex; align-items: center;justify-content: center;margin-top: -40px;"), # Zoom slider)               
                           column(12,align="center",div(style='width:100%;overflow-x: scroll;height:100%;overflow-y: scroll;', uiOutput('scalableplotM')))),
                  #grVizOutput("plot1"),
                  tags$style(
                    '#testM {
    cursor: grap;
    color: black;
    }'),
                  conditionalPanel(id="testM",condition= "input.plot1M_click",wellPanel(textOutput("celltype_defM"))),
                  br(),
                  br(),
                  titlePanel(shiny::span(p(style="text-align: justify;", HTML("<h>Table options</h>"),actionButton('helpM', 'Info',icon= icon("info"), align="left")))),
                  fluidRow(column(width=12,wellPanel(id="sidebar2M",
                                                     # fluidRow(column(3,radioButtons('EC_scoreM','ECs', c(">=1",">=2",">=3",">=4",">=5",">=6",">=7"), selected = ">=1")),
                                                     #          column(3,radioButtons('specificityM','SPs', c(">=0",">=0.25",">=0.5","=1"),selected = ">=0")),
                                                     #          column(3,radioButtons("tabletypeM","Table type",c("Simple","Complete"),selected="Simple")))))),
                  
                                                     
                                                      fluidRow(column(3,sliderInput('EC_scoreM','ECs',min = 1, max = 17, value = 1, step=1)),
                                                               column(3,sliderInput('specificityM','SPs', min = 0, max = 1, value = 0, step=0.1)),
                                                               column(3,radioButtons("tabletypeM","Table type",c("Simple","Complete"),selected="Simple")))))),                                  
                                                                                        
                  
                  fluidRow(
                    column(12, 
                           downloadButton("downloadDataM", "Download", 
                                          style = "padding: 10px 20px; font-size: 18px;")), # Bigger button
                    column(12, 
                           radioButtons("downloadTypeM", "",  # No title
                                        choices = c("CSV" = ".csv",
                                                    "XLSX" = ".xlsx",
                                                    "TSV" = ".tsv"), 
                                        inline = TRUE),
                           tags$style(HTML("
                #downloadTypeM .shiny-options-group label {
             font-size: 12px !important;
           }
           #downloadTypeM {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))) # CSS to make choices smaller
                  ),
                  
                  tags$head(tags$style(HTML('
         #sidebar2M {
            background-color: #ad000019;
            font-size: 16px;
        }'))),
                  # bsTooltip("EC_score", "Evidence consensus score (number of databases)", placement = "top", trigger = "hover",
                  #           options = NULL),
                  br(),
                  br(),
                  dataTableOutput('table1M', width = "100%"),
                  br()
              )),



      # Custom markers integration ----
      tabItem(tabName="integration",
              tags$style(css),
              tags$h2(
                "Integrate custom set of marker genes with the Cell Marker Accordion database",
                style = "font-weight: bold; text-align: center;"
              ),

              # div(fluidRow(wellPanel(id="sidebar",fluidRow(column(8,splitLayout(cellWidths = c("75%", "25%"),fileInput("usermarker", "Load your custom set to be integrated with the Accordion database",buttonLabel=list(icon("upload")),multiple = FALSE),
              #                                                   actionButton('usermarkerinfo', 'InputFile',icon= icon("file-circle-question"), align="left",style='margin-top:30px'))),
              #                                              column(3,radioButtons("databaseInt", "Select the Accordion database:",
              #                                                                    choices = c("Healthy" = "Healthy", "Disease" = "Disease"), inline=TRUE))),
              #                        fluidRow(column(width = 2, actionButton("demo_ex", HTML("Load demo example"),style="font-size: 14px;",icon= icon("cloud-arrow-up")),uiOutput("success_icon_int"), align="left"),
              #                                 column(width = 2, actionButton("start_int", HTML("Start the integration!"),style="font-size: 14px;",icon= icon("gear")), align="left"))))),
              div(
                fluidRow(
                  wellPanel(
                    id = "sidebar",
                    div(
                      style = "display: flex; justify-content: center; align-items: center;",
                      actionButton(
                        'intInfo', 'Info',
                        icon = icon("info"),
                        style = 'margin-top: 5px; margin-bottom: 5px;'
                      )
                    ),
                    # File input and database selection
                    fluidRow(
                      column(
                        8,
                        splitLayout(
                          cellWidths = c("75%", "25%"),
                          fileInput(
                            "usermarker",
                            HTML("Upload custom marker genes"),
                            buttonLabel = list(icon("upload")),
                            multiple = FALSE
                          ),
                          actionButton(
                            'usermarkerinfo',
                            'InputFile',
                            icon = icon("file-circle-question"),
                            align = "left",
                            style = 'margin-top:30px'
                          )
                        )
                      ),
                      column(
                        3,
                        radioButtons(
                          "databaseInt",
                          "Select the Accordion database:",
                          choices = c("Healthy" = "Healthy", "Disease" = "Disease"),
                          inline = TRUE
                        )
                      )
                    ),

                    fluidRow(
                      column(
                        width = 4,
                        div(
                          style = "display: flex; align-items: center; gap: 10px;",
                          actionButton(
                            "demo_ex",
                            HTML("Load demo example"),
                            style = "font-size: 14px;",
                            icon = icon("cloud-arrow-up")
                          ),
                          uiOutput("success_icon_int")  # Checkmark icon appears here
                        )
                      ),
                      
                      column(
                        width = 8,
                        div(
                          style = "display: flex; align-items: center; height: 100%;",  # Ensures vertical alignment
                          conditionalPanel(
                            condition = "input.usermarker > 0 | input.demo_ex > 0",
                            div(
                              style = "display: flex; align-items: center; margin-top: -5px;",  # Adjust this value if needed
                              actionButton(
                                "start_int",
                                HTML("Start the integration!"),
                                style = "font-size: 18px;",
                                icon = icon("gear")
                              )
                            )
                          )
                        )
                      )
                    )
                
                  )
                )
              ),

              conditionalPanel(condition = "input.start_int > 0",
              
                        fluidRow(column(12,radioButtons("tabletypeInt0","Table type",c("Simple","Complete"),selected="Simple",inline = TRUE))),
              
                       
                       fluidRow(
                         column(12, 
                                downloadButton("downloadDataInt", "Download", 
                                               style = "padding: 10px 20px; font-size: 18px;")), # Bigger button
                         column(12, 
                                radioButtons("downloadTypeInt", "",  # No title
                                             choices = c("CSV" = ".csv",
                                                         "XLSX" = ".xlsx",
                                                         "TSV" = ".tsv"), 
                                             inline = TRUE),
                                tags$style(HTML("
                #downloadTypeInt .shiny-options-group label {
             font-size: 12px !important;
           }
           #downloadTypeInt {
             margin-top: -20px; /* Reduces space between elements */
           }
         ")))), # CSS to make choices smaller
                                                                      br(),
                                                      dataTableOutput('table1Int', width = "100%"),
              br(),
              br(),


              useShinyjs(),
              h2("Filter the integrated database"),
              div(fluidRow(column(width=8,wellPanel(id="sidebar",
                                                    div(
                                                      style = "display: flex; justify-content: center; align-items: center;",
                                                      actionButton(
                                                        'celltypeInfo', 'Info',
                                                        icon = icon("info"),
                                                        style = 'margin-top: 5px; margin-bottom: 5px;'
                                                      )
                                                    ),
                                                    checkboxGroupInput("speciesInt", "Select species:",
                                                                       choiceNames =
                                                                         list(tags$img(src = "human.png"),tags$img(src = "mouse.png")),
                                                                       choiceValues =
                                                                         list("Human","Mouse"),selected = c("Human"),inline=TRUE),
                                                    br(),
                                                    #pickerInput('disease', 'Hematologic disorder', choices= NULL ,multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE)),
                                                    pickerInput('diseaseInt', 'Condition',  choices= NULL, selected = "healthy", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    conditionalPanel(condition = "output.showTissuePanel == true",

                                                    pickerInput('tissueInt', 'Tissue',  choices= NULL,selected=tissue_list, multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    #checkboxInput("collapse_tissue","",value=FALSE))),
                                                    checkboxInput("tissue_awareInt","Tissue aware",value=TRUE)),


                                                    pickerInput('celltypeInt', 'Cell type', choices= NULL ,multiple=TRUE,selected=celltype_list, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                                                    br(),
                                                    pickerInput('descendantsofInt', 'See subtypes of:', choices= NULL,multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes"),choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), 141))),
                                                    checkboxInput("cellidInt","Plot celltype_ID",value=FALSE)))),
         
                          
                        
                  tags$head(tags$style(HTML('
         #sidebar {
            background-color: #ad000019;
            font-size: 16px;
        }'))),

                  fluidRow(column(12,   downloadButton("save_plotInt", "Save Ontology plot",style = "padding: 10px 20px; font-size: 18px;")), column(12, radioButtons("file_formatInt", "", choices = c("PNG", "PDF"), inline=TRUE),
                                                                                                                                                   
                                                                                                                                                   tags$style(HTML("
                #file_formatInt .shiny-options-group label {
             font-size: 12px !important;
           }
           #file_formatInt {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))),   
                           column(12, titlePanel(
                             shiny::span(
                               style = "display: flex; align-items: center; justify-content: center; margin-top: -50px;",
                               icon("hand-pointer", class = "about-icon", lib = "font-awesome"),
                               p(style = "text-align: center; ", HTML("<h><strong>Click</strong> on a node to look at cell type description</h>"))
                             )
                           )),  column(12, sliderInput("zoomInt", "Zoom:", min = 0.5, max = 10, value = 1, step = 0.2),style = "display: flex; align-items: center;justify-content: center;margin-top: -40px;"), # Zoom slider)               
                           column(12,align="center",div(style='width:100%;overflow-x: scroll;height:100%;overflow-y: scroll;', uiOutput('scalableplotInt')))),
                  
                  #grVizOutput("plot1"),
                  tags$style(
                    '#testInt {
    cursor: grap;
    color: black;
    }'),
                  conditionalPanel(id="testInt",condition= "input.plot1Int_click",wellPanel(textOutput("celltype_defInt"))),
                  br(),
                  br(),

                  conditionalPanel(condition= "input.descendantsofInt != ''",titlePanel(shiny::span(p(style="text-align: justify;", HTML("<h>Table options</h>"),actionButton('helpInt', 'Info',icon= icon("info"), align="left"))))),
                  conditionalPanel(condition= "input.descendantsofInt == ''",titlePanel(shiny::span(p(style="text-align: justify;", HTML("<h>Table options</h>"),actionButton('help_emptyInt', 'Info',icon= icon("info"), align="left"))))),

                  fluidRow(column(width=12,wellPanel(id="sidebar2",
                                                    
                                                              fluidRow(column(3,sliderInput('EC_scoreInt','ECs',min = 1, max = 17, value = 1, step=1)),
                                                                       column(3,sliderInput('specificityInt','SPs', min = 0, max = 1, value = 0, step=0.1)),
                                                                       column(3,radioButtons("tabletypeInt","Table type",c("Simple","Complete"),selected="Simple")),          
                                                                       
                                                              conditionalPanel(condition= "input.descendantsofInt != ''", column(2,radioButtons("mergeDescendantInt","Merge subtypes", c("Yes","No"),selected="No"))))))),


                  fluidRow(
                    column(12, 
                           downloadButton("downloadDataInt2", "Download", 
                                          style = "padding: 10px 20px; font-size: 18px;")), # Bigger button
                    column(12, 
                           radioButtons("downloadTypeInt2", "",  # No title
                                        choices = c("CSV" = ".csv",
                                                    "XLSX" = ".xlsx",
                                                    "TSV" = ".tsv"), 
                                        inline = TRUE),
                           tags$style(HTML("
                #downloadTypeInt2 .shiny-options-group label {
             font-size: 12px !important;
           }
           #downloadTypeInt2 {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))) # CSS to make choices smaller
                  ),

                  tags$head(tags$style(HTML('
         #sidebar2 {
            background-color: #ad000019;
            font-size: 16px;
        }'))),
                  # bsTooltip("EC_score", "Evidence consensus score (number of databases)", placement = "top", trigger = "hover",
                  #           options = NULL),
                  br(),
                  br(),
                  dataTableOutput('filteredInt', width = "100%")),
                  br()
              )),

      # enrichment with Fisher test-----

      tabItem(tabName="anno",
              tags$style(css),
              tags$h2(
                "Perform cell type marker enrichment analysis across tissues in health and disease",
                style = "font-weight: bold; text-align: center;"
              ),
              div(fluidRow(column(width=6,wellPanel(id="sidebar",
                                                    
                                                    div(
                                                      style = "display: flex; justify-content: center; align-items: center;",
                                                      actionButton(
                                                        'annoInfo', 'Info',
                                                        icon = icon("info"),
                                                        style = 'margin-top: 5px; margin-bottom: 5px;'
                                                      )
                                                    ),
                                                    splitLayout(cellWidths = c("75%", "25%"),fileInput("clusterfile", "Load file",buttonLabel=list(icon("upload")),multiple = FALSE),
                                                                actionButton('userclusterfileinfo', 'InputFile',icon= icon("file-circle-question"), align="left", style='margin-top:30px; margin-bottom:0px')),
                                                    
                                                    fluidRow(column(width=6,radioButtons("nmarkerpos", HTML("Specificy the number of <em>positive</em> genes to keep for each cluster"), choices = nmarker_values)),
                                                             column(width=6, radioButtons("nmarkerneg", HTML("Specificy the number of <em>negative</em> genes to keep for each cluster"), choices = nmarker_values))),
                                                    fluidRow(column(width=6,textInput("nmarkerotherpos", HTML("Type in number of <em>positive</em> genes"))),
                                                             column(width=6,textInput("nmarkerotherneg", HTML("Type in number of <em>negative</em> genes")))),
                                                    fluidRow(column(width=6,actionButton("addpos", HTML("Add value"))),
                                                             column(width=6,actionButton("addneg", HTML("Add value")))),
                                                    #checkboxInput("collapse_tissue","",value=FALSE))),
                                                    

                                                    style = "padding-bottom: 5px;")),
                           column(
                             6, offset = 0,  # Centers the column in a 12-grid layout
                             div(
                               style = "display: flex; align-items: center; justify-content: center; gap: 10px; height: 100%;",
                               actionButton(
                                 "demo_ex_anno",
                                 HTML("Load demo example"),
                                 style = "font-size: 14px;",
                                 icon = icon("cloud-arrow-up")
                               ),
                               div(style = "display: flex; align-items: center;", uiOutput("success_icon_anno"))  # Keeps icon aligned
                             )
                           ),
                           
                           
                           
                            br(),
                            column(6, HTML("Output table generated by the FindAllMarkers function in the Seurat package"), align="center"),
                           br(),
                           column(6,  dataTableOutput('exampletable', width = "100%"))),
                  br(),

                  fluidRow(
                    column(
                      width = 6,
                      wellPanel(
                        id = "sidebar",

                        # Centered and larger "Filters" title
                        div(style = "text-align: center; font-size: 24px; font-weight: bold; margin-bottom: 10px;",
                            HTML("Filters"),
                            actionButton('filterhelpA', 'Info', icon = icon("info"))
                        ),

                        checkboxGroupInput("speciesA", "Select species:",
                                           choiceNames =
                                             list(tags$img(src = "human.png"),tags$img(src = "mouse.png")),
                                           choiceValues =
                                             list("Human","Mouse"),selected = c("Human"),inline=TRUE),
                        pickerInput('diseaseA', 'Condition',  choices= disease_list, selected = "healthy", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                        pickerInput('tissueA', 'Tissue',  choices= unique(accordion_complete[DO_diseasetype == "healthy"]$Uberon_tissue),selected="blood", multiple=TRUE, options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                        pickerInput('celltypeA', 'Cell type', choices= unique(accordion_complete[DO_diseasetype == "healthy" & Uberon_tissue=="blood"]$celltype) ,multiple=TRUE,selected=unique(accordion_complete[DO_diseasetype == "healthy" & Uberon_tissue=="blood"]$celltype), options = list(`actions-box` = TRUE,`live-search`=TRUE, style="box-celltypes")),
                        
                          br(),
                          fluidRow(column(4,sliderInput('EC_scoreA','ECs',min = 1, max = 17, value = 1, step=1)),
                                   column(4,sliderInput('specificityA','SPs', min = 0, max = 1, value = 0, step=0.1)),
                        
                          
                          
                          column(width = 4,
                                 textInput("max_n_marker",
                                           value = 50,
                                           placeholder = NULL,
                                           label = HTML("Maximum number of markers to keep for each cell type"))
                          )
                        )
                      )
                    ),




                  column(6,div(img(src = "Logo.png",height="40%", width="40%"),style="text-align: center;"),
                         br(),
                         actionButton("button", HTML("<strong>Start enrichment!</strong>"),style="font-size: 20px;",icon= icon("rocket")), align="center")),

                  conditionalPanel(
                    condition= "input.button > 0",           
  
                  fluidRow(
                    column(12, 
                           downloadButton("downloadDataA", "Download", 
                                          style = "padding: 10px 20px; font-size: 18px;")), # Bigger button
                    column(12, 
                           radioButtons("downloadTypeA", "",  # No title
                                        choices = c("CSV" = ".csv",
                                                    "XLSX" = ".xlsx",
                                                    "TSV" = ".tsv"), 
                                        inline = TRUE),
                           tags$style(HTML("
                #downloadTypeA .shiny-options-group label {
             font-size: 12px !important;
           }
           #downloadTypeA {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))) # CSS to make choices smaller
                  ),
                  br(),
                  dataTableOutput('table1A', width = "100%"),
                  br(),
          

                    checkboxInput("cellidA","Plot celltype_ID",value=FALSE),
                    #Ontology plot

                  fluidRow(column(12,   downloadButton("save_plotA", "Save Ontology plot",style = "padding: 10px 20px; font-size: 18px;")), column(12, radioButtons("file_formatA", "", choices = c("PNG", "PDF"), inline=TRUE),
                                                                                                                                                     
                                                                                                                                                     tags$style(HTML("
                #file_formatA .shiny-options-group label {
             font-size: 12px !important;
           }
           #file_formatA {
             margin-top: -20px; /* Reduces space between elements */
           }
         "))),
                           column(12, titlePanel(
                             shiny::span(
                               style = "display: flex; align-items: center; justify-content: center; margin-top: -50px;",
                               icon("hand-pointer", class = "about-icon", lib = "font-awesome"),
                               p(style = "text-align: center; ", HTML("<h><strong>Click</strong> on a node to look at cell type description</h>"))
                             )
                           )),  column(12, sliderInput("zoomA", "Zoom:", min = 0.5, max = 10, value = 1, step = 0.2),style = "display: flex; align-items: center;justify-content: center;margin-top: -40px;"), # Zoom slider)               
                           column(12,align="center",div(style='width:100%;overflow-x: scroll;height:100%;overflow-y: scroll;', uiOutput('scalableplotA')))),
                  
                  #grVizOutput("plot1"),
                  tags$style(
                    '#testA {
                              cursor: grap;
                              color: black;
                            }'),
                  conditionalPanel(id="testA",condition= "input.plot1A_click",wellPanel(textOutput("celltype_defA"))),
                  br(),
                  tags$head(tags$style(HTML('
         #sidebar {
            background-color: #ad000019;
            font-size: 16px;
        }'))),



                  tags$head(tags$style(HTML('
         #sidebar2A {
            background-color: #ad000019;
            font-size: 16px;
        }'))),
                  # bsTooltip("EC_score", "Evidence consensus score (number of databases)", placement = "top", trigger = "hover",
                  #           options = NULL),
                  br(),
                  br()
              ))
    )
  )
  )
)
server <- function(input, output, session) {

  # server for cell type search ----
  # Download handler
  output$downloadAccordionDB <- downloadHandler(
    filename = function() {
      "TheCellMarkerAccordion_database_v0.9.5.xlsx"  # Name of the downloaded file
    },

    content = function(file) {
      withProgress(message = "Preparing file...", value = 0, {
      # Simulate preparation (e.g., loading, processing)
      incProgress(0.1, detail = "Downloading The Cell Marker Accordion database...")  # Initial step (10%)
      Sys.sleep(1)  # Simulate processing time

      # Load data
      load("./data/AccordionDB_list.rda")
      incProgress(0.2, detail = "Downloading The Cell Marker Accordion database...")  # Update progress (20%)
      Sys.sleep(1)  # Simulate some more processing time

      # Write the file
      write_xlsx(data_list, file)
      incProgress(0.5, detail = "Downloading The Cell Marker Accordion database...")  # Another update (50%)
      Sys.sleep(1)  # Simulate the writing process

      # Finish the progress (total progress = 1)
      incProgress(0.2, detail = "Downloading The Cell Marker Accordion database...")  # Final step (20%)
      })
    }
  )


  observeEvent(input$celltypeInfo,{
    showModal(modalDialog(
      title = "Inputs information",
      HTML(     "
<ul>
  <li><strong>Select species:</strong> currently Human and/or Mouse.</li>
  <li><strong>Condition:</strong> healthy or multiple diseases.</li>
  <li><strong>Tissue:</strong> select one or multiple tissues from the list. When the <strong>'Tissue aware'</strong> button is enabled, tissue specificity is maintained (i.e., tissues remain separate). Otherwise, selected tissues will be combined and analyzed together.</li>
  <li><strong>Cell type:</strong> select one or multiple cell types from the list.</li>
  <li><strong>See subtypes of:</strong> displays the list of cell type descendants of the previously selected cell types. Users can select one or more subtypes to visualize in the Ontology tree and the output table.</li>
</ul>
"
        )
        ))
  })
  
  observeEvent(input$help,{
    showModal(modalDialog(
      title = "Table Option Information",
      HTML("Marker genes can be ranked and selected by: <br>
      <ul><li> <strong> ECs</strong>: evidence consistency score, measuring the agreement of different annotation sources </li>
       <li> <strong> SPs</strong>: specificity score, indicating whether a gene is a marker for different cell types present in all the accordion database </li>

      In addition:  <br>
      <strong> Merged subtypes </strong>:  if subtypes of at least one input cell type is displayed
      <br> <ul><li> yes: merge together the subtypes gene markers and assign them to selected cell type </li><li> no: merge of subtypes is not performed </li></ul>
      <strong> Table type </strong> <br> <ul><li><strong>simple</strong>: retrieves a compact and easy table format</li><li> <strong>complete</strong>: additional information are added</li></ul>")
      #textInput('text2', 'You can also put UI elements here')
    ))
  })


  observeEvent(input$help_empty,{
    showModal(modalDialog(
      title = "Table Option Information",
      HTML("Marker genes can be ranked and selected by: <br>
             <ul><li> <strong> ECs</strong>: evidence consistency score, measuring the agreement of different annotation sources </li>
             <li> <strong> SPs</strong>: specificity score, indicating whether a gene is a marker for different cell types present in all the accordion database </li>
             In addition: <br>
      <strong> Table type </strong> <br> <ul><li><strong>simple</strong>: retrieves a compact and easy table format</li><li> <strong>complete</strong>: additional information are added</li></ul>")
      #textInput('text2', 'You can also put UI elements here')
    ))
  })

  
  toListen_disease <- reactive({
    list(input$species)
  })
  
  observeEvent(toListen_disease(),{
    updatePickerInput(session,'disease', selected="healthy",
                      choices=unique(accordion_complete[species %in% input$species]$DO_diseasetype),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(accordion_complete[species %in% input$species]$DO_diseasetype))))
  })
  
  
  markerTableComplete <- reactive({
    validate(need(length(input$disease) >=1, "Please select at least one condition"))
    accordion_complete<-accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species]
    accordion_complete

  })

  toListen_tissue <- reactive({
    list(input$disease,input$species)
  })

  observeEvent(toListen_tissue(),{
    updatePickerInput(session,'tissue', selected="blood",
                      choices=unique(accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species]$Uberon_tissue),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species]$Uberon_tissue))))
  })

  toListen <- reactive({
    list(input$species, input$disease, input$tissue)
  })

  observeEvent(toListen(),{
    updatePickerInput(session,'celltype', selected=c("hematopoietic stem cell (Hs, Mm)", "hematopoietic multipotent progenitor cell (Hs, Mm)","hematopoietic lineage restricted progenitor cell (Hs, Mm)"),
                      choices=unique(accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species & Uberon_tissue %in% input$tissue]$celltype_species),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species & Uberon_tissue %in% input$tissue]$celltype))))
  })


  controlDescendant<-reactive({
    if("healthy" %in% input$disease){
      onto_igraph<-graph_from_graphnel(onto_plot, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
      node<-as.data.table(unique(markerTableComplete()[celltype_species %in% input$celltype & Uberon_tissue %in% input$tissue & DO_diseasetype=="healthy"]$celltype))
      node_present<-V(onto_igraph)$name
      node_present <- gsub("CL:\\d+", "", gsub("\n", " ", node_present))
      # Trim extra spaces
      node_present <- trimws(node_present)      
      dt<- data.table(celltype=character(), distance=numeric())
      for (i in 1:nrow(node)){
        if (node[i]$V1 %in% node_present) {
        distan<-max(eccentricity(onto_igraph, vids = node[i]$V1, mode = c("out")))
        data<-as.data.table(t(c(node[i]$V1,distan)))
        colnames(data)<-c("celltype","distance")
        dt<-rbind(dt,data)
        }
      }  
      as.data.table(dt)
    } else{
      data.table(celltype="", distance =0)
    }
  })


  observeEvent(input$celltype,{
    updatePickerInput(session,'descendantsof',
                      choices=unique(accordion_complete[DO_diseasetype %in% input$disease & species %in% input$species & Uberon_tissue %in% input$tissue & celltype %in% controlDescendant()[distance!=0]$celltype]$celltype_species),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), 105)))
  })

  descendantTable<-reactive ({
    if(length(input$descendantsof)>0){
      select_descendant(onto_plot, markerTableComplete(), input$descendantsof)
    }

  })

  mergeDescendantTable<-reactive ({
    if(length(input$descendantsof)>0 & input$mergeDescendant=="Yes"){
      table_merge_descendant(onto_plot, cell_onto, markerTableComplete(), input$descendantsof, input$species, input$tissue, colnames(accordion_complete))
    }
  })

  markerTablePreOutput <- reactive ({

    #table with descandants merged
    if(length(input$descendantsof)!=0 & input$mergeDescendant=="Yes"){
      table_merged_desc_other(markerTableComplete(),input$celltype,mergeDescendantTable(),input$species,input$tissue,colnames(accordion_complete))
    }
    #table with descendants NOT merged
    else if (length(input$descendantsof)!=0 & input$mergeDescendant=="No"){
      table_desc_other(markerTableComplete(),input$celltype,descendantTable(),input$species,input$tissue,colnames(accordion_complete))
    } #table without descendant
    else if(length(input$descendantsof)==0){
      table_input_celltypes(markerTableComplete(),input$celltype,input$species,input$tissue, input$disease)
    }

  })

  markerTableOutput <- reactive ({
    if(nrow(markerTablePreOutput() > 0)){

      #keep tissue separated
      if(input$tissue_aware == TRUE){
        #compute EC and specificity condition&tissue specific

        #EC score
        st_table_tissue_specific<-unique(markerTablePreOutput()[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID,celltype,celltype_ID,marker,marker_type),nrow)
        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(markerTablePreOutput(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]

        mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type", "specificity_score")

        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition", "original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      } else{ #not consider tissue
        #EC_score
        st_table_tissue_specific<-unique(markerTablePreOutput()[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]

        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,celltype,celltype_ID,marker,marker_type),nrow)

        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(markerTablePreOutput(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID, marker,marker_type),nrow)

        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","marker","marker_type", "specificity_score")

        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]
        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      }
    }
  })

  markerTablePlot <- reactive ({
    if(nrow(markerTableComplete()>0)){
      #plot with descendants
      if (length(input$descendantsof)!=0){
        table_desc_other(markerTableComplete(),input$celltype,descendantTable(),input$species,input$tissue, colnames(accordion_complete))
      }
      else if(length(input$descendantsof)==0){
        table_input_celltypes(markerTableComplete(),input$celltype,input$species,input$tissue, input$disease)
      }
    }
  })

  Plot <- reactive({
    
    if(!is.null(markerTablePlot())){
      healthy_ct<-unique(markerTablePlot()[DO_diseasetype=="healthy"]$celltype_ID)
    
      if("healthy" %in% input$disease){
        if (length(healthy_ct) > 1 | length(input$descendantsof)>=1){
          ontosubplot<-onto_plot2(cell_onto,healthy_ct,cex=0.8)
          nodes<-as.data.table(ontosubplot@nodes)
          nodes<-nodes[,V1:=tstrsplit(nodes$V1,"CL", fixed = TRUE, keep = 1)]
          nodes<-nodes[,V1:=tstrsplit(nodes$V1,"DOID", fixed = TRUE, keep = 1)]

          nodes<-nodes[, V1:=str_replace_all(V1,"\n"," ")]
          nodes<-nodes[, V1:=str_sub(V1,1,nchar(V1)-1)]
          ontosubplot@nodes<-nodes$V1
          ontosubplot
        } else if(length(healthy_ct) == 1){
          ontosubplot<-onto_plot(cell_onto,healthy_ct,cex=0.8)
          ontosubplot2<-make_graphNEL_from_ontology_plot(ontosubplot)
          dt_onto2<-as.data.table(ontosubplot2@nodes)
          label<-merge(dt_onto2,ontology_celltype,by.x="V1",by.y="celltype_ID")
          ontosubplot[["node_attributes"]][["label"]][[label$V1]]<-label$celltype
          ontosubplot2@nodes<-ontosubplot[["node_attributes"]][["label"]]

          ontosubplot2

        }
      }
    }
  })


  # Reactive expression to store the Graphviz plot
  graphPlot <- reactive({
    if(!is.null(markerTablePlot())) {
      healthy_ct <- markerTablePlot()[DO_diseasetype == "healthy"]$celltype_ID
      if (length(healthy_ct) >= 1) {
        if (length(input$descendantsof) >= 1) {
          return(hierac_plot1_desc(markerTableComplete(), Plot(), healthy_ct, descendantTable(), input$cellid, input$disease))
        } else {
          return(hierac_plot1(markerTableComplete(), Plot(), healthy_ct, input$cellid, input$disease))
        }
      }
    } else {
    return(NULL)  #Return NULL if no plot is available
    }#
  })
  
  output$plot1 <- renderGrViz({

      
      graphPlot()  # Render the reactive graphPlot()

    })
  

  output$scalableplot <- renderUI({
    tagList( tags$div(
      #div(grVizOutput('plot1',height = input$height, width = input$width)))
     grVizOutput('plot1', height = "700px", width = "700px"),
     style = paste("transform: scale(", input$zoom, "); transform-origin: center;")
    )
    )
  })

  
  output$save_plot <- downloadHandler(
    filename = function() {
      if (input$file_format == "PNG") {
      "Ontology_plot.png"
      } else if(input$file_format == "PDF") {
        "Ontology_plot.pdf"
      }
    },
    content = function(file) {
      grViz_obj <- graphPlot()  # Use the stored reactive Graphviz plot
      
      if (is.null(grViz_obj)) {
        stop("No plot available to save.")
      }
      
      svg_file <- tempfile(fileext = ".svg")
      
      # Convert Graphviz to SVG and save
      svg_content <- export_svg(grViz_obj)  # Export Graphviz as SVG
      writeLines(svg_content, svg_file)  # Write to temporary SVG file
      
      # Get actual width from UI
      default <-700

      # Scale dimensions based on zoom input
      scaled_zoom <- default * input$zoom

      
      # Convert SVG to PNG with specified width & height
      if (input$file_format == "PNG") {
      suppressWarnings(
        rsvg_png(svg_file, file, width = scaled_zoom, height = scaled_zoom)
      
      )} else if(input$file_format == "PDF"){
        suppressWarnings(
        rsvg_pdf(svg_file, file, width = scaled_zoom / 100, height = scaled_zoom / 100) 
        )}
    }
  )

  

  click_plot<-reactive ({
    click_node(accordion_complete,ontology_celltype,Plot(),markerTableComplete(),input$cellid, ontology_def, input$disease)
  })


  txt <- reactive({
    req(input$plot1_click)
    nodeval <- input$plot1_click$nodeValues[1:length(input$plot1_click$nodeValues)]
    nodeval<-as.data.table(paste(nodeval,sep="", collapse=" "))
    return(click_plot()[label %in% nodeval$V1]$cell_definition)

  })
  output$celltype_def <- renderText({
    req(txt())
    txt()
  })


  outputTable<-reactive({
    if(!is.null(markerTableOutput())){
      if("celltype" %in% colnames(markerTableOutput())){
        markerTableOutput()[DO_diseasetype=="healthy",CL_celltype:=celltype][DO_diseasetype=="healthy",CL_ID:=celltype_ID][DO_diseasetype=="healthy",CL_cell_definition:=cell_definition]
        markerTableOutput()[DO_diseasetype!="healthy",NCIT_celltype:=celltype][DO_diseasetype!="healthy",NCIT_ID:=celltype_ID][DO_diseasetype!="healthy",NCIT_cell_definition:=cell_definition]
        markerTableOutput()[,c("celltype","celltype_ID","cell_definition"):=NULL]
      }
      if("EC_score" %in% colnames(markerTableOutput())){
        setnames(markerTableOutput(), "EC_score", "ECs")
        setnames(markerTableOutput(), "specificity_score", "SPs")
      }

      if(length(input$celltype)!=0 & input$tabletype=="Complete" & input$tissue_aware==TRUE){
        if(input$mergeDescendant=="Yes"){
          markerTableOutput()[ECs>= input$EC_score& SPs >= input$specificity][,
                                                                                                                                     c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","celltype_ancestor","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }else{
          markerTableOutput()[ECs>= input$EC_score & SPs >= input$specificity][,
                                                                                                                                     c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }
      } else if(length(input$celltype)!=0 & input$tabletype=="Complete" & input$tissue_aware==FALSE){
        if(input$mergeDescendant=="Yes"){
          markerTableOutput()[ECs>= input$EC_score & SPs >= input$specificity][,
                                                                                                                                     c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","celltype_ancestor","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }else{
          markerTableOutput()[ECs>= input$EC_score & SPs >= input$specificity][,
                                                                                                                                     c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }
      }  else if(length(input$celltype)!=0 & input$tabletype=="Simple" & input$tissue_aware==TRUE){
        if(input$mergeDescendant=="Yes"){
          unique(markerTableOutput()[ECs>= input$EC_score & SPs >= input$specificity][,
                                                                                                                                            c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype_ancestor","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }else{
          unique(markerTableOutput()[ECs>=input$EC_score & SPs >= input$specificity][,
                                                                                                                                            c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }
      } else if(length(input$celltype)!=0 & input$tabletype=="Simple" & input$tissue_aware==FALSE){
        if(input$mergeDescendant=="Yes"){
          unique(markerTableOutput()[ECs>= input$EC_score & SPs >=input$specificity][,
                                                                                                                                            c("species","DO_diseasetype","DO_ID","celltype_ancestor","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }else{
          unique(markerTableOutput()[ECs>= input$EC_score & SPs >= input$specificity][,
                                                                                                                                            c("species","DO_diseasetype","DO_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }
      }
    }

  })

  outputTable2<-reactive({
    if(!is.null(outputTable())){
      if(all(is.na(outputTable()$DO_ID))){
        suppressWarnings({
          outputTable()[,c("original_diseasetype","DO_diseasetype","DO_ID","DO_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition"):=NULL]
        })
      }
      if(all(is.na(outputTable()$CL_ID))){
        suppressWarnings({
          outputTable()[,c("CL_celltype","CL_ID","CL_cell_definition"):=NULL]
        })
      }

      outputTable()
    }

  })

  output$table1 <- renderDataTable({
    if(!is.null(outputTable2())){
      datatable(outputTable2(),rownames = FALSE,options = list(
        pageLength=5, scrollX='400px',autoWidth = TRUE,
        columnDefs = list(
          list(width = '130px', targets = "_all"))), filter = 'top')
  }
    })


  # Downloadable csv of selected dataset

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("TheCellMarkerAccordion_database_v0.9.5", input$downloadType)
    },
    content = function(file) {
      if(input$downloadType == ".csv"){
        write.csv(outputTable(), file, row.names = FALSE, quote=FALSE)
      }
      else if(input$downloadType == ".xlsx") {
        write_xlsx(outputTable(), file)
      }
      else if(input$downloadType == ".tsv") {
        write.table(outputTable(), file, quote = FALSE,
                    sep='\t', row.names = FALSE)
      }

    }


  )





  # server for marker search ----
  
  observeEvent(input$markerInfo,{
    showModal(modalDialog(
      title = "Inputs information",
      HTML("
<ul>
  <li><strong>Select species:</strong> currently Human and/or Mouse.</li>
  <li><strong>Insert marker genes:</strong> enter a list of marker genes, using <code>| , : ; ! ?</code> as delimiters to separate them.</li>
  <li><strong>Upload file with marker genes:</strong> provide a <code>.txt</code>, <code>.csv</code>, <code>.xlsx</code>, or <code>.tsv</code> file containing a list of marker genes. Use <code>| , : ; ! ?</code> as delimiters to separate them. Both the markers entered in the \"Insert marker genes\" box and those in the uploaded file will be considered.</li>
  <li><strong>Condition:</strong> healthy or multiple diseases.</li>
  <li><strong>Tissue:</strong> select one or multiple tissues from the list. When the <strong>'Tissue aware'</strong> button is enabled, tissue specificity is maintained (i.e., tissues remain separate). Otherwise, the selected tissues will be combined and analyzed together.</li>
</ul>
"

      )
    ))
  })
  
  observeEvent(input$helpM,{
    showModal(modalDialog(
      title = "Table Option Information",
      HTML("Marker genes can be ranked and selected by: <br>
             <ul><li> <strong> ECs</strong>: evidence consistency score, measuring the agreement of different annotation sources </li>
             <li> <strong> SPs</strong>: specificity score, indicating whether a gene is a marker for different cell types present in all the accordion database </li></ul>
             In addition: <br>
      <strong> Table type </strong> <br> <ul><li><strong>simple</strong>: retrieves a compact and easy table format</li><li> <strong>complete</strong>: additional information are added</li></ul>")
    ))
  })




  #table based on selected genes-marker
  tableMarkerInput<-reactive ({
    if(length(input$marker)!=0){
      marker_vec<-as.data.frame(unlist(strsplit(input$marker, "[\\|, +:;!?]+")))
      colnames(marker_vec)<-"marker_gene"
      marker_input<-tolower(marker_vec$marker_gene)
      marker_dt<-tolower(accordion_complete$marker)
      validate(need(marker_input %in% marker_dt, "Please insert valid marker genes!"))
      gene_marker<-accordion_complete[tolower(marker) %in% tolower(marker_vec$marker_gene)]
      #gene_marker<-gene_marker[species %in% input$speciesM]
      gene_marker
      # if(length(input$diseaseM>0)){
      #   gene_marker<- gene_marker[DO_diseasetype %in% input$diseaseM]
      #   gene_marker
      # }

    }
  })

  #table base on input file genes marker
  tableFileMarker<-reactive ({
    file_load<-input$markerfile
    fileName <- file_load$datapath
    mark<-readChar(fileName, file.info(fileName)$size)
    marker_vec<-as.data.frame(unlist(strsplit(mark, "[\\|, \t\n\r;:]+")))
    colnames(marker_vec)<-"marker_gene"
    table_marker_file<-accordion_complete[tolower(marker) %in% tolower(marker_vec$marker_gene)]
    #table_marker_file <- table_marker_file[species %in% input$speciesM]
    # if(length(input$diseaseM>0)){
    #   table_marker_file<- table_marker_file[DO_diseasetype %in% input$diseaseM]
    # 
    # }
    table_marker_file
  })

  # Final table output:
  #1. Marker insert by the user in the box
  #2. Marker add through a marker file
  #3. If both marker in the box and file marker are present add together

  tableInputComplete_nodis_notissue <- reactive({
    req(input$speciesM)  # Ensure input$speciesM exists
    
    if (nchar(input$marker) > 1 & !is.null(input$markerfile)) {
      table_tot <- rbind(tableMarkerInput(), tableFileMarker())
    } else if (nchar(input$marker) > 1 & is.null(input$markerfile)) {
      table_tot <- tableMarkerInput()
    } else if (nchar(input$marker) == 0 & !is.null(input$markerfile)) {
      table_tot <- tableFileMarker()
    } else {
      return(NULL)  # Return NULL if none of the conditions match
    }
    
    table_tot <- table_tot[species %in% input$speciesM]
    
    return(table_tot)
  })
  
  toListen_diseaseM <- reactive({
    list(input$speciesM, tableInputComplete_nodis_notissue())
  })
  
  observeEvent(toListen_diseaseM(), {
    if (is.null(tableInputComplete_nodis_notissue())) {
      # Set pickerInput to empty when table is NULL
      updatePickerInput(session, 'diseaseM',  selected = character(0), choices = character(0))
    } else {
      updatePickerInput(
        session,
        'diseaseM',
        selected = unique(tableInputComplete_nodis_notissue()$DO_diseasetype),
        choices = unique(tableInputComplete_nodis_notissue()$DO_diseasetype),
        options = list(`actions-box` = TRUE, style = "box-celltypes"),
        choicesOpt = list(style = rep("font-size: 18px; line-height: 1.6;", 
                                      uniqueN(tableInputComplete_nodis_notissue()$DO_diseasetype)))
      )
    }
  })
  
  tableInputComplete_notissue <- reactive({
    req(tableInputComplete_nodis_notissue(), input$diseaseM)
    table_filtered <- tableInputComplete_nodis_notissue()[DO_diseasetype %in% input$diseaseM]
    return(table_filtered)
  })
  
  toListen_tissueM <- reactive({
    list(input$speciesM, input$diseaseM, tableInputComplete_nodis_notissue())
  })
  
  observeEvent(toListen_tissueM(), {
    if (is.null(tableInputComplete_nodis_notissue())) {
      # Set pickerInput to empty when table is NULL
      updatePickerInput(session, 'tissueM', selected = character(0), choices = character(0))
    } else {
      updatePickerInput(
        session,
        'tissueM',
        selected = unique(tableInputComplete_notissue()$Uberon_tissue),
        choices = unique(tableInputComplete_notissue()$Uberon_tissue),
        options = list(`actions-box` = TRUE, style = "box-celltypes"),
        choicesOpt = list(style = rep("font-size: 18px; line-height: 1.6;", 
                                      uniqueN(tableInputComplete_notissue()$Uberon_tissue)))
      )
    }
  })
  tableInputComplete <- reactive({
    req(tableInputComplete_notissue(), input$tissueM)
    
    if (nrow(tableInputComplete_notissue()) > 0) {
      table_final <- tableInputComplete_notissue()[Uberon_tissue %in% input$tissueM]
      return(table_final)
    } else {
      return(NULL)
    }
  })

  PlotM <- reactive({
    healthy_ct<-tableInputComplete()[DO_diseasetype =="healthy"]$celltype_ID
    if (length(healthy_ct)>=1 & "healthy" %in% input$diseaseM) {
      cell_matching<-unique(tableInputComplete()$celltype_ID)
      if(length(cell_matching)>1){
        ontosubplot<-onto_plot2(cell_onto,cell_matching ,cex=0.8)
        nodes<-as.data.table(ontosubplot@nodes)
        nodes<-nodes[,V1:=tstrsplit(nodes$V1,"CL", fixed = TRUE, keep = 1)]
        nodes<-nodes[,V1:=tstrsplit(nodes$V1,"DOID", fixed = TRUE, keep = 1)]
        nodes<-nodes[, V1:=str_replace_all(V1,"\n"," ")]
        nodes<-nodes[, V1:=str_sub(V1,1,nchar(V1)-1)]
        ontosubplot@nodes<-nodes$V1
        ontosubplot
      } else{
        ontosubplot<-onto_plot(cell_onto,cell_matching ,cex=0.8)
        ontosubplot2<-make_graphNEL_from_ontology_plot(ontosubplot)
        dt_onto2<-as.data.table(ontosubplot2@nodes)
        label<-merge(dt_onto2,ontology_celltype,by.x="V1",by.y="celltype_ID")
        ontosubplot[["node_attributes"]][["label"]][[label$V1]]<-label$celltype
        ontosubplot2@nodes<-ontosubplot[["node_attributes"]][["label"]]
        ontosubplot2

      }

    }

  })


  
  # Reactive expression to store the Graphviz plot
  graphPlotM <- reactive({
    if(!is.null(tableInputComplete())){
      if(nrow(tableInputComplete())>0){
        healthy_ct<-tableInputComplete()[DO_diseasetype=="healthy"]$celltype_ID
        if(length(healthy_ct>=1)){
          return(hierac_plot2(accordion_complete, PlotM(), tableInputComplete(),input$cellidM, input$diseaseM))
        }
      }
    } else{
      return(NULL)
    }
  })
  
  output$plot1M <- renderGrViz({
    graphPlotM()  # Use the reactive graphPlot()
  })
  

  output$scalableplotM <- renderUI({
    tagList( tags$div(
      #div(grVizOutput('plot1',height = input$height, width = input$width)))
      grVizOutput('plot1M', height = "700px", width = "700px"),
      style = paste("transform: scale(", input$zoomM, "); transform-origin: center;")
    )
    )
  })
  
  output$save_plotM <- downloadHandler(
    filename = function() {
      if (input$file_formatM == "PNG") {
        "Ontology_plot.png"
      } else if(input$file_formatM == "PDF") {
        "Ontology_plot.pdf"
      }
    },
    content = function(file) {
      grViz_obj <- graphPlotM()  # Use the stored reactive Graphviz plot
      
      if (is.null(grViz_obj)) {
        stop("No plot available to save.")
      }
      
      svg_file <- tempfile(fileext = ".svg")
      
      # Convert Graphviz to SVG and save
      svg_content <- export_svg(grViz_obj)  # Export Graphviz as SVG
      writeLines(svg_content, svg_file)  # Write to temporary SVG file
      
      # Get actual width from UI
      default <-700
      
      # Scale dimensions based on zoom input
      scaled_zoom <- default * input$zoomM
      
      
      # Convert SVG to PNG with specified width & height
      if (input$file_format == "PNG") {
        suppressWarnings(
          rsvg_png(svg_file, file, width = scaled_zoom, height = scaled_zoom)
          
        )} else if(input$file_format == "PDF"){
          suppressWarnings(
            rsvg_pdf(svg_file, file, width = scaled_zoom / 100, height = scaled_zoom / 100) 
          )}
    }
  )
  
  
  

  click_plotM<-reactive ({
    click_node(accordion_complete, ontology_celltype, PlotM(), tableInputComplete(),input$cellidM, ontology_def, input$diseaseM)
  })


  txtM <- reactive({
    req(input$plot1M_click)
    nodeval <- input$plot1M_click$nodeValues[1:length(input$plot1M_click$nodeValues)]
    nodeval<-as.data.table(paste(nodeval,sep="", collapse=" "))
    return(  click_plotM()[label %in% nodeval$V1]$cell_definition)

  })
  output$celltype_defM <- renderText({
    req(txtM())
    txtM()
  })




  markerTableOutputM <- reactive ({
    if(!is.null(tableInputComplete())){
      if(nrow(tableInputComplete())){


      #keep tissue separated
      if(input$tissue_awareM == TRUE){
        #compute EC and specificity condition&tissue specific

        #EC score
        st_table_tissue_specific<-unique(tableInputComplete()[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID,celltype,celltype_ID,marker,marker_type),nrow)
        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(tableInputComplete(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]


        mark_spec<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID, marker,marker_type),nrow)

        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type", "specificity_score")
        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition", "original_tissue","Uberon_tissue","Uberon_ID","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      } else{ #not consider tissue
        #EC_score
        st_table_tissue_specific<-unique(tableInputComplete()[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]

        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,celltype,celltype_ID,marker,marker_type),nrow)

        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(tableInputComplete(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID, marker,marker_type),nrow)

        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","marker","marker_type", "specificity_score")

        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]
        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      }
      final_table<-unique(final_table)
      final_table
      }
    }
  })

  FinalTableM<-reactive({
    if(!is.null(markerTableOutputM())){
      if(nrow(markerTableOutputM())){
      if("celltype" %in% colnames(markerTableOutputM())){
        markerTableOutputM()[DO_diseasetype=="healthy",CL_celltype:=celltype][DO_diseasetype=="healthy",CL_ID:=celltype_ID][DO_diseasetype=="healthy",CL_cell_definition:=cell_definition]
        markerTableOutputM()[DO_diseasetype!="healthy",NCIT_celltype:=celltype][DO_diseasetype!="healthy",NCIT_ID:=celltype_ID][DO_diseasetype!="healthy",NCIT_cell_definition:=cell_definition]
        markerTableOutputM()[,c("celltype","celltype_ID","cell_definition"):=NULL]
      }
      if("EC_score" %in% colnames(markerTableOutputM())){
        setnames(markerTableOutputM(), "EC_score", "ECs")
        setnames(markerTableOutputM(), "specificity_score", "SPs")
      }

      if(nrow(markerTableOutputM())!=0 & input$tabletypeM=="Complete"){
        unique(markerTableOutputM()[ECs >= input$EC_scoreM & SPs >= input$specificityM])
      } else if(nrow(markerTableOutputM())!=0 & input$tabletypeM=="Simple" & input$tissue_awareM == TRUE){
        unique(markerTableOutputM()[ECs>= input$EC_scoreM & SPs >= input$specificityM
        ][,c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
      } else if(nrow(markerTableOutputM())!=0 & input$tabletypeM=="Simple" & input$tissue_awareM == FALSE){
        unique(markerTableOutputM()[ECs>= input$EC_scoreM & SPs >= input$specificityM
        ][,c("species","DO_diseasetype","DO_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])

      }


    }
    }
  })


  outputTable2M<-reactive({
    if(!is.null(FinalTableM())){
      if(nrow(FinalTableM()>0)){

      if(all(is.na(FinalTableM()$DO_ID))){
        suppressWarnings({
          FinalTableM()[,c("original_diseasetype","DO_diseasetype","DO_ID","DO_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition"):=NULL]
        })
      }
      if(all(is.na(FinalTableM()$CL_ID))){
        suppressWarnings({
          FinalTableM()[,c("CL_celltype","CL_ID","CL_cell_definition"):=NULL]
        })
      }

      FinalTableM()
      }
    } else{
      showNotification("No genes were found with the selected filters", type = "warning", duration = 5)

    }

  })

  output$table1M <- renderDataTable({
    if(!is.null(outputTable2M())){
      datatable(outputTable2M(),rownames = FALSE,options = list(
        pageLength=5, scrollX='400px',columnDefs = list(
          list(width = '130px', targets = "_all"))), filter = 'top')

    }
  })
  # Downloadable csv of selected dataset

  output$downloadDataM <- downloadHandler(
    filename = function() {
      paste0("TheCellMarkerAccordion_database_v0.9.5", input$downloadTypeM)
    },
    content = function(file) {
      if(input$downloadTypeM == ".csv"){
        write.csv(FinalTableM(), file, row.names = FALSE, quote=FALSE)
      }
      else if(input$downloadTypeM == ".xlsx") {
        write_xlsx(FinalTableM(), file)
      }
      else if(input$downloadTypeM == ".tsv") {
        write.table(FinalTableM(), file, quote = FALSE,
                    sep='\t', row.names = FALSE)
      }

    }

  )



  # server for custom markers integration with the accordion database -----

  observeEvent(input$intInfo,{
    showModal(modalDialog(
      title = "Inputs information",
      HTML("Upload a custom set of marker genes to be integrated with the Cell Marker Accordion database, either healthy or disease.
      Click the  <strong>'InputFile'</strong> button to access details on the input file format. <br>
      <strong> Select the Accordion database </strong> <br> <ul><li><strong>Healthy</strong>: integrate a custom set with the healthy Accordion database</li><li> <strong>Disease</strong>: integrate a custom set with the disease Accordion database</li></ul>
      Click the <strong>'Load demo example'</strong> button to load an example of custom marker gene set for integration. To download the example set, click the <strong> InputFile</strong> button and use the <strong> Download Example</strong> option.   <br>
      To start the integration press on the 'Start the integration!' button.
      Once the integration is performed you can explore and download the integrated table.
      <strong> Table type </strong> <br> <ul><li><strong>simple</strong>: retrieves a compact and easy table format</li><li> <strong>complete</strong>: additional information are added</li></ul>")
      #textInput('text2', 'You can also put UI elements here')
    ))
  })
  
  observeEvent(input$usermarkerinfo,{
    showModal(modalDialog(
      title = "Marker genes input file",
      HTML("Upload a custom set of marker genes to be integrated with the Cell Marker Accordion database, either healthy or disease. <br>
      The file must contain at least two columns:
      <ul>
        <li> <strong> cell_type </strong>:  specifies the cell type.
          <br>To ensure proper integration, cell types nomenclature should be standardized on the Cell Ontology or the NCI Thesaurus. <br>
          If non-standardized cell types are provided, they will be added as new cell types in the database.
        </li>
        <li> <strong> marker </strong>: lists the marker genes </li>
      </ul>
      Additional columns can also be included:
      <ul>
        <li> <strong> species </strong>: Specifies the species (default: Human).</li>
        <li> <strong> tissue </strong>: Specifies the related tissue. Standardization with Uberon Ontology is recommended for effective integration.<br>
        Non-standardized tissues will be added as new tissues. If omitted, integration will ignore tissue specificity.
        </li>
        <li> <strong> marker_type </strong>: Defines marker type (positive or negative; default: positive).</li>
        <li> <strong> resource </strong>: Indicates the data source. If omitted, markers are labeled as custom_set.</li>
        <li> <strong> disease </strong>: Required if the integration is performed with the disease database. Standardization with Disease Ontology is recommended. <br>
        Non-standardized diseases will be added as new diseases. If omitted, disease specificity is ignored.</li>
      </ul>
      <br><strong> Example </strong>"
      ),
      #HTML("<img src=input_file_example.jpg>"),
      downloadButton("download_ex_int", "Download Example"),
      HTML("<br>"),
      HTML('
        <table style="border-collapse: collapse; width: 100%; text-align: center; background-color: white;">
          <thead>
            <tr style="background-color: #f2f2f2; text-align: center;">
              <th style="border: 1px solid black; padding: 8px; text-align: center;">species</th>
              <th style="border: 1px solid black; padding: 8px; text-align: center;">tissue</th>
              <th style="border: 1px solid black; padding: 8px; text-align: center;">cell_type</th>
              <th style="border: 1px solid black; padding: 8px; text-align: center;">marker</th>
              <th style="border: 1px solid black; padding: 8px; text-align: center;">marker_type</th>
              <th style="border: 1px solid black; padding: 8px; text-align: center;">resource</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">glutamatergic neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Satb2</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">glutamatergic neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Satb2</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_2</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">glutamatergic neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Slc17a6</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">glutamatergic neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Slc17a7</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">glutamatergic neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Slc17a7</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_2</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">pyramidal neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Sv2b</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">pyramidal neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Calb1</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">pyramidal neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Pde1a</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_1</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">pyramidal neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Pde1a</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_2</td>
            </tr>
            <tr>
              <td style="border: 1px solid black; padding: 8px;">Mouse</td>
              <td style="border: 1px solid black; padding: 8px;">brain</td>
              <td style="border: 1px solid black; padding: 8px;">pyramidal neuron</td>
              <td style="border: 1px solid black; padding: 8px;">Pde1a</td>
              <td style="border: 1px solid black; padding: 8px;">positive</td>
              <td style="border: 1px solid black; padding: 8px;">custom_set_3</td>
            </tr>
          </tbody>
        </table>
      '),
      easyClose = TRUE,
      footer = NULL,
      size="l"

    ))
  })


  # Server logic for download
  output$download_ex_int <- downloadHandler(
    filename = function() {
      "Custom_set_integration_example.xlsx"
    },
    content = function(file) {
      table<-read_excel("data/Demo_example_integration.xlsx")
      write_xlsx(table, file)
    }
  )

  AccordionToIntegrate<-reactive({
    if("Healthy" %in% input$databaseInt){
      load("data/accordion_marker.rda")
      accordion_to_int<-accordion_marker
      accordion_to_int$SPs_tissue_specific<-NULL
      accordion_to_int$SPs_global<-NULL
      accordion_to_int$ECs_tissue_specific<-NULL
      accordion_to_int$ECs_global<-NULL

    } else{
      load("data/disease_accordion_marker.rda")
      accordion_to_int<-disease_accordion_marker
      accordion_to_int$SPs_disease_tissue_specific<-NULL
      accordion_to_int$SPs_disease_specific<-NULL
      accordion_to_int$SPs_global<-NULL
      accordion_to_int$ECs_disease_tissue_specific<-NULL
      accordion_to_int$ECs_disease_specific<-NULL
      accordion_to_int$ECs_global<-NULL
    }
    accordion_to_int
  })



  integratedDB<-reactive({
    if(!is.null(input$usermarker)){
      file_load_marker<-input$usermarker
      fileName_marker <- file_load_marker$datapath
      ext <-  tools::file_ext(fileName_marker)
      req(file_load_marker)
      validate(need(ext %in% c("csv", "xlsx","txt"), "Please upload a csv, txt or xlsx file"))
      if(ext=="xlsx"){
        user_inputfile<-read_excel(fileName_marker)
      }
      if(ext=="csv" | ext == "txt" | ext == "tsv"){
        user_inputfile<-fread(fileName_marker)

      }
    } else if(input$demo_ex >0 & is.null(input$usermarker)){
      user_inputfile<-read_excel("data/Demo_example_integration.xlsx")
    }
    if(!is.null(input$usermarker) | input$demo_ex >0){
      user_inputfile<-as.data.table(user_inputfile)
      user_inputfile<-unique(user_inputfile)
      validate(need(ncol(user_inputfile) >=2, "Insufficient number of columns! Need at least cell type and marker columns"))
      if("Healthy" %in% input$databaseInt){
        if(!"cell_type" %in% colnames(user_inputfile)){
          validate(need("cell_type" %in% colnames(user_inputfile),"\"cell_type\" column not found. Please specificy the \"cell_type\" columns contaning the list of cell types to integrate."))
        } else{
          setnames(user_inputfile, "cell_type", "CL_celltype")
        }
        if(!"marker" %in% colnames(user_inputfile)){
          validate(need("marker" %in% colnames(user_inputfile),"\"marker\" column not found.  Please specificy the \"marker\" column contaning the list of markers to integrate."))
        }
        if(!"species" %in% colnames(user_inputfile)){
          showNotification("\"species\" column not found. By default Human will be considered.", type = "warning", duration = 5)
          user_inputfile[,species:="Human"]
        }

        if(!"tissue" %in% colnames(user_inputfile)){
          showNotification("\"tissue\" column not found. The integration will not consider the tissues.", type = "warning", duration = 5)
          user_inputfile[,Uberon_tissue:=NA][,original_tissue:=NA]
          AccordionToIntegrate()[,Uberon_tissue:=NA][,Uberon_ID:=NA][,original_tissue:=NA]
        } else{
          setnames(user_inputfile, "tissue", "Uberon_tissue")
          user_inputfile[,original_tissue:=Uberon_tissue]

        }
        if(!"marker_type" %in% colnames(user_inputfile)){
          showNotification("\"marker_type\" column not found. By default genes will be considered as positive markers.", type = "warning", duration = 5)
          user_inputfile[,marker_type:="positive"]
        }
        if(!"resource" %in% colnames(user_inputfile)){
          showNotification("\"resource\" column not found. By default a unique resource will be associated with possible recurrence markers.", type = "warning", duration = 5)
          user_inputfile[,resource:="custom_set"]
        }
        user_inputfile[,original_celltype:=CL_celltype]
        user_inputfile[,log2FC:=NA]
        user_inputfile[,p.value:=NA]
        user_inputfile[,adjusted_p.value:=NA]
        user_inputfile[,pct1:=NA]
        #AccordionToIntegrate()[,c("SPs_tissue_specific","SPs_global","ECs_tissue_specific","ECs_global"):=NULL]
        user_inputfile<-unique(user_inputfile)

        #add cell_ID


        ontology_celltype<-as.data.frame(cell_onto[["name"]])
        colnames(ontology_celltype)<-"CL_celltype"
        ontology_celltype$CL_ID<-rownames(ontology_celltype)
        ontology_celltype$CL_cell_definition<-cell_onto[["def"]]
        ontology_celltype<-as.data.table(ontology_celltype)

        user_inputfile<-merge(user_inputfile, ontology_celltype, by="CL_celltype", all.x=TRUE)

        #add Uberon_ID
        load("data/uberon_onto.rda")
        ontology_tissue<-as.data.frame(uberon_onto[["name"]])
        colnames(ontology_tissue)<-"Uberon_tissue"
        ontology_tissue$Uberon_ID<-rownames(ontology_tissue)
        ontology_tissue$tissue_definition<-uberon_onto[["def"]]
        ontology_tissue<-as.data.table(ontology_tissue)
        user_inputfile<-merge(user_inputfile, ontology_tissue, by="Uberon_tissue", all.x=TRUE)

        #add gene description
        load("data/gene_description.rda")
        user_inputfile<-merge(user_inputfile, gene_description, by.x="marker",by.y="gene_symbol", all.x = TRUE)

        user_inputfile<-user_inputfile[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1")]
        user_inputfile<-unique(user_inputfile)
        user_inputfile[,CL_cell_definition:=str_remove_all(CL_cell_definition, '"')]
        user_inputfile[,CL_cell_definition:=str_remove_all(CL_cell_definition, "\\$")]
        user_inputfile[,CL_cell_definition:=str_remove_all(CL_cell_definition, r"(\\)")]
        user_inputfile[,CL_cell_definition:=tstrsplit(CL_cell_definition, "[", fixed = TRUE, keep = 1)]
        user_inputfile[,gene_description:=tstrsplit(gene_description,"[",fixed=TRUE,keep=1)]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, '"')]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, "\\$")]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, r"(\\)")]
        user_inputfile[,tissue_definition:=tstrsplit(tissue_definition, "[", fixed = TRUE, keep = 1)]
        marker_database<-rbind(AccordionToIntegrate(), user_inputfile)


      }
      if("Disease" %in% input$databaseInt){

        if(!"cell_type" %in% colnames(user_inputfile)){
          validate(need("cell_type" %in% colnames(user_inputfile),"\"cell_type\" column not found. Please specificy the \"cell_type\" columns contaning the list of cell types to integrate."))
        } else{
          setnames(user_inputfile, "cell_type", "NCIT_celltype")
        }
        if(!"marker" %in% colnames(user_inputfile)){
          validate(need("marker" %in% colnames(user_inputfile),"\"marker\" column not found. Please specificy the \"marker\" column contaning the list of markers to integrate."))
        }
        if(!"disease" %in% colnames(user_inputfile)){
          showNotification("\"disease\" column not found. The integration will not consider the disease.", type = "warning", duration = 5)
          user_inputfile[,DO_diseasetype:=NA][,DO_ID:=NA][,original_diseasetype:=NA][,DO_definition:=NA]
          AccordionToIntegrate()[,DO_diseasetype:=NA][,DO_ID:=NA][,original_diseasetype:=NA][,DO_definition:=NA]
        } else{
          setnames(user_inputfile, "disease", "DO_diseasetype")
          user_inputfile[,original_diseasetype:=DO_diseasetype][,DO_ID:=NA][,DO_definition:=NA]
        }
        if(!"species" %in% colnames(user_inputfile)){
          showNotification("\"species\" column not found. By default Human will be considered.", type = "warning", duration = 5)
          user_inputfile[,species:="Human"]
        }

        if(!"tissue" %in% colnames(user_inputfile)){
          showNotification("\"tissue\" column not found. The integration will not consider the tissues.", type = "warning", duration = 5)
          user_inputfile[,Uberon_tissue:=NA][,original_tissue:=NA]
          AccordionToIntegrate()[,Uberon_tissue:=NA][,Uberon_ID:=NA][,original_tissue:=NA]
        } else{
          setnames(user_inputfile, "tissue", "Uberon_tissue")
          user_inputfile[,original_tissue:=Uberon_tissue]

        }
        if(!"marker_type" %in% colnames(user_inputfile)){
          showNotification("\"marker_type\" column not found. By default genes will be considered as positive markers.", type = "warning", duration = 5)
          user_inputfile[,marker_type:="positive"]
        }
        if(!"resource" %in% colnames(user_inputfile)){
          showNotification("\"resource\" column not found. By default a unique resource will be associated with possible recurrence markers.", type = "warning", duration = 5)
          user_inputfile[,resource:="custom_set"]
        }

        user_inputfile[,NCIT_cell_definition:=NA]
        user_inputfile[,NCIT_ID:=NA]
        user_inputfile[,original_celltype:=NCIT_celltype]
        user_inputfile[,log2FC:=NA]
        user_inputfile[,p.value:=NA]
        user_inputfile[,adjusted_p.value:=NA]
        user_inputfile[,pct1:=NA]

        #AccordionToIntegrate()[,c("SPs_disease_tissue_specific","SPs_disease_specific","SPs_global","ECs_disease_tissue_specific","ECs_disease_specific","ECs_global"):=NULL]

        user_inputfile<-unique(user_inputfile)

        #Add Uberon_ID
        load("data/uberon_onto.rda")
        ontology_tissue<-as.data.frame(uberon_onto[["name"]])
        colnames(ontology_tissue)<-"Uberon_tissue"
        ontology_tissue$Uberon_ID<-rownames(ontology_tissue)
        ontology_tissue$tissue_definition<-uberon_onto[["def"]]
        ontology_tissue<-as.data.table(ontology_tissue)
        user_inputfile<-merge(user_inputfile, ontology_tissue, by="Uberon_tissue", all.x=TRUE)


        load("data/gene_description.rda")

        user_inputfile<-merge(user_inputfile, gene_description, by.x="marker",by.y="gene_symbol", all.x=TRUE)
        user_inputfile<-unique(user_inputfile)
        user_inputfile[,NCIT_cell_definition:=str_remove_all(NCIT_cell_definition, '"')]
        user_inputfile[,NCIT_cell_definition:=str_remove_all(NCIT_cell_definition, "\\$")]
        user_inputfile[,NCIT_cell_definition:=str_remove_all(NCIT_cell_definition, r"(\\)")]
        user_inputfile[,NCIT_cell_definition:=tstrsplit(NCIT_cell_definition, "[", fixed = TRUE, keep = 1)]
        user_inputfile[,gene_description:=tstrsplit(gene_description,"[",fixed=TRUE,keep=1)]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, '"')]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, "\\$")]
        user_inputfile[,tissue_definition:=str_remove_all(tissue_definition, r"(\\)")]
        user_inputfile[,tissue_definition:=tstrsplit(tissue_definition, "[", fixed = TRUE, keep = 1)]

        user_inputfile<-user_inputfile[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1")]
        user_inputfile<-unique(user_inputfile)

        marker_database<-rbind(AccordionToIntegrate(), user_inputfile)


      }
      marker_database<-unique(marker_database)
      marker_database

    }
  })

  demoLoaded_int <- reactiveVal(FALSE)  # Track if demo is loaded

  observeEvent(input$demo_ex, {
    showNotification("Demo example loaded!", type = "message", duration = 5)
    demoLoaded_int(TRUE)  # Mark demo as loaded
  })

  output$success_icon_int <- renderUI({
    if (demoLoaded_int()) {
      tags$span(
        icon("check-circle", class = "text-success", lib = "font-awesome"),
        style = "font-size: 20px; margin-left: 10px;"
      )
    }
  })
  #recalculate EC score and specificity-perform integration
  newTable<-reactive({
    if(!is.null(integratedDB())){
      if(input$start_int){
      tissue<-unique(integratedDB()$Uberon_tissue)
      withProgress(message ="Integrating custom markers with The Cell Marker Accordion database...", value = 0,{
      if(!all(is.na(tissue))){ #tissue aware
          if(input$databaseInt == "Healthy"){
            #compute EC and specificity condition&tissue specific
            #EC score
            st_table_tissue_specific<-unique(integratedDB()[,c("Uberon_tissue","Uberon_ID","species","CL_celltype","CL_ID","marker","marker_type","resource")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(CL_celltype)]
            incProgress(0.2, detail = "Calculating ECs...")
            ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,Uberon_tissue,Uberon_ID,CL_celltype,CL_ID,marker,marker_type),nrow)
            colnames(ec_score_tissue_specific)<-c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker", "marker_type","EC_score")
            accordion_ec_table<-merge(integratedDB(),ec_score_tissue_specific,by=c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker", "marker_type"), all.x = TRUE)

            #specificity
            incProgress(0.3, detail = "Calculating SPs...")
            st_table_tissue_specific<-unique(accordion_ec_table[,c("Uberon_tissue","Uberon_ID","species","CL_celltype","CL_ID","marker","marker_type","marker","marker_type")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(CL_celltype)]
            mark_spec<-ddply(st_table_tissue_specific,.(species,Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
            colnames(mark_spec)<-c("species","Uberon_tissue","Uberon_ID","marker","marker_type", "specificity_score")

            final_table<-merge(accordion_ec_table,mark_spec,by=c("species","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
            final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

            final_table<-final_table[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          } else{
            #compute EC and specificity condition&tissue specific

            #EC score
            st_table_tissue_specific<-unique(integratedDB()[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(NCIT_celltype)]
            incProgress(0.2, detail = "Calculating ECs...")

            ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
            colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","EC_score")
            accordion_ec_table<-merge(integratedDB(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker", "marker_type"), all.x = TRUE)

            #specificity
            incProgress(0.3, detail = "Calculating SPs...")

            st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","marker","marker_type")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(NCIT_celltype)]

            mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
            colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type", "specificity_score")

            final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
            final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

            final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition", "original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          }
        } else{
          if(input$databaseInt == "Healthy"){
            #compute EC and specificity condition&tissue specific

            #EC score
            st_table_tissue_specific<-unique(integratedDB()[,c("species","CL_celltype","CL_ID","marker","marker_type","resource")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(CL_celltype)]
            incProgress(0.2, detail = "Calculating ECs...")

            ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,CL_celltype,CL_ID,marker,marker_type),nrow)
            colnames(ec_score_tissue_specific)<-c("species","CL_celltype","CL_ID","marker", "marker_type","EC_score")
            accordion_ec_table<-merge(integratedDB(),ec_score_tissue_specific,by=c("species","CL_celltype","CL_ID","marker", "marker_type"), all.x = TRUE)

            #specificity
            incProgress(0.3, detail = "Calculating SPs...")

            st_table_tissue_specific<-unique(accordion_ec_table[,c("species","CL_celltype","CL_ID","marker","marker_type","marker","marker_type")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(CL_celltype)]

            mark_spec<-ddply(st_table_tissue_specific,.(species, marker,marker_type),nrow)
            colnames(mark_spec)<-c("species","marker","marker_type", "specificity_score")

            final_table<-merge(accordion_ec_table,mark_spec,by=c("species","marker","marker_type"),all.x = TRUE)
            final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

            final_table<-final_table[,c("species","original_celltype","CL_celltype","CL_ID","CL_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          } else{
            #compute EC and specificity condition&tissue specific

            #EC score
            st_table_tissue_specific<-unique(integratedDB()[,c("DO_diseasetype","DO_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","resource")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(NCIT_celltype)]
            incProgress(0.2, detail = "Calculating ECs...")

            ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,NCIT_celltype,NCIT_ID,marker,marker_type),nrow)
            colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","EC_score")
            accordion_ec_table<-merge(integratedDB(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","NCIT_celltype","NCIT_ID","marker", "marker_type"), all.x = TRUE)

            #specificity
            incProgress(0.3, detail = "Calculating SPs...")

            st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","species","NCIT_celltype","NCIT_ID","marker","marker_type","marker","marker_type")])
            st_table_tissue_specific<-st_table_tissue_specific[!is.na(NCIT_celltype)]

            mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID, marker,marker_type),nrow)
            colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","marker","marker_type", "specificity_score")

            final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","marker","marker_type"),all.x = TRUE)
            final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

            final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          }
        }
      incProgress(0.2, detail = "Integration completed!")
      final_table
      })
      }
    }
  })


  tableType<-reactive({
    if(!is.null(newTable())){
      tissue<-unique(integratedDB()$Uberon_tissue)
      if(input$databaseInt == "Healthy"){
          if("celltype" %in% colnames(newTable())){
          setnames(newTable(), "celltype","CL_celltype")
          setnames(newTable(), "celltype_ID","CL_ID")
          setnames(newTable(), "cell_definition","CL_cell_definition")
          }
        if(input$tabletypeInt0 == "Simple"){
          if(!all(is.na(tissue))){ #tissue aware
            table<-newTable()[,c("species","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","marker","marker_type","EC_score","specificity_score")]

          }else {
            table<-newTable()[,c("species","CL_celltype","CL_ID", "marker","marker_type","EC_score","specificity_score")]
          }
        } else {
          if(!all(is.na(tissue))){ #tissue aware
            table<-newTable()[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
            
          }else {
            table<-newTable()[,c("species","original_celltype","CL_celltype","CL_ID","CL_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          }
        }
      } else if (input$databaseInt == "Disease"){
        if("celltype" %in% colnames(newTable())){
          setnames(newTable(), "celltype","NCIT_celltype")
          setnames(newTable(), "celltype_ID","NCIT_ID")
          setnames(newTable(), "cell_definition","NCIT_cell_definition")
        }
        if(input$tabletypeInt0 == "Simple"){
          if(!all(is.na(tissue))){ #tissue aware
            #print(colnames(newTable()))
            table<-newTable()[,c("species","Uberon_tissue","Uberon_ID","NCIT_celltype","NCIT_ID","marker","marker_type","EC_score","specificity_score")]
            
          }else {
            table<-newTable()[,c("species","NCIT_celltype","NCIT_ID", "marker","marker_type","EC_score","specificity_score")]
          }
        } else {
          if(!all(is.na(tissue))){ #tissue aware
            #print(colnames(newTable()))
            table<-newTable()[,c("species","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
            
          }else {
            table<-newTable()[,c("species","original_celltype","NCIT_celltype","NCIT_ID","NCIT_cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
          }
        }
      }
          
      setnames(table, "EC_score","ECs")
      setnames(table, "specificity_score","SPs")
      table
      }
  })
  
  output$table1Int <- renderDataTable({
    if(!is.null(tableType())){
      datatable(tableType(),rownames = FALSE,options = list(
        pageLength=5, scrollX='400px',columnDefs = list(
          list(width = '130px', targets = "_all"))), filter = 'top')
    }
  })


  # Downloadable csv of selected dataset

  output$downloadDataInt <- downloadHandler(
    filename = function() {
      paste0("TheIntegratedCellMarkerAccordion_database_v0.9.5", input$downloadTypeInt)
    },
    content = function(file) {
      if(input$downloadTypeInt == ".csv"){
        write.csv(integratedDB(), file, row.names = FALSE, quote=FALSE)
      }
      else if(input$downloadTypeInt == ".xlsx") {
        write_xlsx(integratedDB(), file)
      }
      else if(input$downloadTypeInt == ".tsv") {
        write.table(integratedDB(), file, quote = FALSE,
                    sep='\t', row.names = FALSE)
      }

    })


  integratedDBmod <-reactive({
    if(!is.null(newTable())){
      if(input$databaseInt == "Healthy"){
        double_species<-unique(newTable()[,c("species","CL_celltype")])
        dup<-double_species[duplicated(CL_celltype)]
        table<-newTable()[, celltype_species:=ifelse(CL_celltype %in% dup$CL_celltype, paste0(CL_celltype, " (Hs, Mm)"),
                                                ifelse(species =="Human", paste0(CL_celltype, " (Hs)"),
                                                       paste0(CL_celltype, " (Mm)")))]



        setnames(table, "CL_celltype","celltype")
        setnames(table, "CL_ID","celltype_ID")
        setnames(table, "CL_cell_definition","cell_definition")

        table[,c("original_diseasetype","DO_ID","DO_definition"):=list(NA)]
        table[,DO_diseasetype:="healthy"]

      } else{

        double_species<-unique(newTable()[,c("species","NCIT_celltype")])
        dup<-double_species[duplicated(NCIT_celltype)]
        table<-newTable()[, celltype_species:=ifelse(NCIT_celltype %in% dup$NCIT_celltype, paste0(NCIT_celltype, " (Hs, Mm)"),
                                                  ifelse(species =="Human", paste0(NCIT_celltype, " (Hs)"),
                                                         paste0(NCIT_celltype, " (Mm)")))]

        setnames(table, "NCIT_celltype","celltype")
        setnames(table, "NCIT_ID","celltype_ID")
        setnames(table, "NCIT_cell_definition","cell_definition")

      }
   
      
      table<-table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition", "original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","celltype","celltype_species","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1")]
      table<-as.data.table(table)
      table
    }
  })



  observeEvent(input$helpInt,{
    showModal(modalDialog(
      title = "Table Option Information",
      HTML("Marker genes can be ranked and selected by: <br>
      <ul><li> <strong> ECs</strong>: evidence consistency score, measuring the agreement of different annotation sources </li>
       <li> <strong> SPs</strong>: specificity score, indicating whether a gene is a marker for different cell types present in all the accordion database </li>

      In addition:  <br>
      <strong> Merged subtypes </strong>:  if subtypes of at least one input cell type is displayed
      <br> <ul><li> yes: merge together the subtypes gene markers and assign them to selected cell type </li><li> no: merge of subtypes is not performed </li></ul>
             <strong> Table type </strong> <br> <ul><li>simple: provides a compact table with fewer columns for easier viewing </li><li> complete: displays the full database, including detailed mapping relationships for diseases, tissues, and cell types to the Disease Ontology, Uberon Ontology, Cell Ontology and NCIT.</li></ul>")
      #textInput('text2', 'You can also put UI elements here')
    ))
  })


  observeEvent(input$help_emptyInt,{
    showModal(modalDialog(
      title = "Table Option Information",
      HTML("Marker genes can be ranked and selected by: <br>
             <ul><li> <strong> ECs</strong>: evidence consistency score, measuring the agreement of different annotation sources </li>
             <li> <strong> SPs</strong>: specificity score, indicating whether a gene is a marker for different cell types present in all the accordion database </li>
             In addition: <br>
             <strong> Table type </strong> <br> <ul><li>simple: provides a compact table with fewer columns for easier viewing </li><li> complete: displays the full database, including detailed mapping relationships for diseases, tissues, and cell types to the Disease Ontology, Uberon Ontology, Cell Ontology and NCIT.</li></ul>")
      #textInput('text2', 'You can also put UI elements here')
    ))
  })

  observeEvent(integratedDBmod(),{
    updatePickerInput(session,'diseaseInt',
                      choices=unique(integratedDBmod()$DO_diseasetype),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(integratedDBmod()$DO_diseasetype))))
  })

  table1 <- reactive({
    if(!is.null(integratedDBmod())){
      validate(need(length(input$diseaseInt) >=1, "Please select at least one condition"))
      integratedDBmod()[DO_diseasetype %in% input$diseaseInt & species %in% input$speciesInt]
    }
  })


  observe({
    if(!is.null(integratedDBmod())){
    tissue_data <- integratedDBmod()[DO_diseasetype %in% input$diseaseInt & species %in% input$speciesInt]$Uberon_tissue
    unique_tissues <- unique(tissue_data)

    # Check if length > 1 and no NA values
    if (length(unique_tissues) == 1 && is.na(unique_tissues)) {
      output$showTissuePanel <- reactive(FALSE)
    } else {
      output$showTissuePanel <- reactive(TRUE)
    }
    # Make sure to invalidate the reactive when necessary
    outputOptions(output, 'showTissuePanel', suspendWhenHidden = FALSE)
    }
  })


  observeEvent(table1(),{
    updatePickerInput(session,'tissueInt',
                      choices=unique(table1()$Uberon_tissue),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(table1()$Uberon_tissue))))
  })


  table2<-reactive({
    if(!is.null(table1())){

      if(length(unique(table1()$Uberon_tissue))==1){
        if(is.na(unique(tabl1()$Uberon_tissue))) {
          filtered_data <- table1()
        }
      } else {
        # Use all cell types if no specific tissue is selected
        filtered_data <- table1()[Uberon_tissue %in% input$tissueInt]
      }
    filtered_data
    }
  })


  observeEvent(table2(), {
    # Get unique cell types
    unique_celltypes <- unique(table2()$celltype_species)

    # Update pickerInput for cell types
    updatePickerInput(
      session,
      'celltypeInt',
      choices = unique_celltypes,
      options = list(`actions-box` = TRUE, style = "box-celltypes"),
      choicesOpt = list(style = rep("font-size: 18px; line-height: 1.6;", length(unique_celltypes)))
    )
  })

  markerTableCompleteInt<-reactive({
    if(!is.null(table2())){

    table2()[celltype_species %in% input$celltypeInt]
    }
  })

  controlDescendantInt<-reactive({
    if(!is.null(markerTableCompleteInt())){
      if("healthy" %in% input$diseaseInt){
        onto_igraph<-graph_from_graphnel(onto_plot, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
        node<-as.data.table(unique(markerTableCompleteInt()[celltype_species %in% input$celltypeInt & Uberon_tissue %in% input$tissueInt & DO_diseasetype=="healthy"]$celltype))
        node_present<-V(onto_igraph)$name
        node_present <- gsub("CL:\\d+", "", gsub("\n", " ", node_present))
        # Trim extra spaces
        node_present <- trimws(node_present)      
        dt<- data.table(celltype=character(), distance=numeric())
        for (i in 1:nrow(node)){
          if (node[i]$V1 %in% node_present) {
            distan<-max(eccentricity(onto_igraph, vids = node[i]$V1, mode = c("out")))
            data<-as.data.table(t(c(node[i]$V1,distan)))
            colnames(data)<-c("celltype","distance")
            dt<-rbind(dt,data)
          }
        }  
        as.data.table(dt)
      } else{
        data.table(celltype="", distance =0)
      }
    }
    })


  observeEvent(input$celltypeInt,{
    updatePickerInput(session,'descendantsofInt',
                      choices=unique(markerTableCompleteInt()[DO_diseasetype %in% input$diseaseInt & species %in% input$speciesInt & Uberon_tissue %in% input$tissueInt & celltype %in% controlDescendantInt()[distance!=0]$celltype]$celltype_species),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), 105)))
  })

  descendantTableInt<-reactive ({
    if(length(input$descendantsofInt)>0){
      select_descendant(onto_plot, markerTableCompleteInt(), input$descendantsofInt)
    }

  })

  mergeDescendantTableInt<-reactive ({
    if(length(input$descendantsofInt)>0 & input$mergeDescendantInt=="Yes"){
      table_merge_descendant(onto_plot, cell_onto, markerTableCompleteInt(), input$descendantsofInt, input$speciesInt, input$tissueInt, colnames(markerTableCompleteInt()))
    }
  })

  markerTablePreOutputInt <- reactive ({
     if(!is.null(markerTableCompleteInt())){
    #table with descandants merged
    if(length(input$descendantsofInt)!=0 & input$mergeDescendantInt=="Yes"){
      table_merged_desc_other(markerTableCompleteInt(),input$celltypeInt,mergeDescendantTableInt(),input$speciesInt,input$tissueInt,colnames(markerTableCompleteInt()))
    }
    #table with descendants NOT merged
    else if (length(input$descendantsofInt)!=0 & input$mergeDescendantInt=="No"){
      table_desc_other(markerTableCompleteInt(),input$celltypeInt,descendantTableInt(),input$speciesInt,input$tissueInt,colnames(markerTableCompleteInt()))
    } #table without descendant
    else if(length(input$descendantsofInt)==0){
      table_input_celltypes(markerTableCompleteInt(),input$celltypeInt,input$speciesInt,input$tissueInt, input$diseaseInt)
    }
    }
  })

  markerTableOutputInt <- reactive ({
    if(!is.null(markerTableCompleteInt())){
    if(nrow(markerTablePreOutputInt() > 0)){
      #keep tissue separated
      if(input$tissue_awareInt == TRUE){
        #compute EC and specificity condition&tissue specific

        #EC score
        st_table_tissue_specific<-unique(markerTablePreOutputInt()[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID,celltype,celltype_ID,marker,marker_type),nrow)
        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(markerTablePreOutputInt(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","species","celltype","celltype_ID","marker","marker_type","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]

        mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID,Uberon_tissue,Uberon_ID, marker,marker_type),nrow)
        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type", "specificity_score")

        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]

        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition", "original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      } else{ #not consider tissue
        #EC_score
        st_table_tissue_specific<-unique(markerTablePreOutputInt()[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type","resource")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]

        ec_score_tissue_specific<-ddply(st_table_tissue_specific,.(species,DO_diseasetype,DO_ID,celltype,celltype_ID,marker,marker_type),nrow)

        colnames(ec_score_tissue_specific)<-c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type","EC_score")
        accordion_ec_table<-merge(markerTablePreOutputInt(),ec_score_tissue_specific,by=c("species","DO_diseasetype","DO_ID","celltype","celltype_ID","marker", "marker_type"), all.x = TRUE)

        #specificity
        st_table_tissue_specific<-unique(accordion_ec_table[,c("DO_diseasetype","DO_ID","species","celltype","celltype_ID","marker","marker_type","marker","marker_type")])
        st_table_tissue_specific<-st_table_tissue_specific[!is.na(celltype)]
        mark_spec<-ddply(st_table_tissue_specific,.(species, DO_diseasetype,DO_ID, marker,marker_type),nrow)

        colnames(mark_spec)<-c("species","DO_diseasetype","DO_ID","marker","marker_type", "specificity_score")

        final_table<-merge(accordion_ec_table,mark_spec,by=c("species","DO_diseasetype","DO_ID","marker","marker_type"),all.x = TRUE)
        final_table[,specificity_score:=format(round(1/specificity_score,2), nsmall=2)]
        final_table<-final_table[,c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","celltype","celltype_ID","cell_definition", "marker","gene_description","marker_type", "resource", "log2FC", "p.value", "adjusted_p.value","pct1", "EC_score","specificity_score")]
      }
    }
    }
  })

  markerTablePlotInt <- reactive ({
    if(!is.null(markerTableCompleteInt())){
      #plot with descendants
      if (length(input$descendantsofInt)!=0){
        table_desc_other(markerTableCompleteInt(),input$celltypeInt,descendantTableInt(),input$speciesInt,input$tissueInt, colnames(markerTableCompleteInt()))
      }
      else if(length(input$descendantsof)==0){
        table_input_celltypes(markerTableCompleteInt(),input$celltypeInt,input$speciesInt,input$tissueInt, input$diseaseInt)
      }
    }
  })

  PlotInt <- reactive({
    if(!is.null(markerTablePlotInt())){
      healthy_ct<-unique(markerTablePlotInt()[DO_diseasetype=="healthy"]$celltype_ID)
      if("healthy" %in% input$diseaseInt){
        if (length(healthy_ct) > 1 | length(input$descendantsofInt)>=1){
          ontosubplot<-onto_plot2(cell_onto,healthy_ct,cex=0.8)
          nodes<-as.data.table(ontosubplot@nodes)
          nodes<-nodes[,V1:=tstrsplit(nodes$V1,"CL", fixed = TRUE, keep = 1)]
          nodes<-nodes[,V1:=tstrsplit(nodes$V1,"DOID", fixed = TRUE, keep = 1)]

          nodes<-nodes[, V1:=str_replace_all(V1,"\n"," ")]
          nodes<-nodes[, V1:=str_sub(V1,1,nchar(V1)-1)]
          ontosubplot@nodes<-nodes$V1
          ontosubplot
        } else if(length(healthy_ct) == 1){
          ontosubplot<-onto_plot(cell_onto,healthy_ct,cex=0.8)
          ontosubplot2<-make_graphNEL_from_ontology_plot(ontosubplot)
          dt_onto2<-as.data.table(ontosubplot2@nodes)
          label<-merge(dt_onto2,ontology_celltype,by.x="V1",by.y="celltype_ID")
          ontosubplot[["node_attributes"]][["label"]][[label$V1]]<-label$celltype
          ontosubplot2@nodes<-ontosubplot[["node_attributes"]][["label"]]

          ontosubplot2

        }
      }
    }
  })






  

  # Reactive expression to store the Graphviz plot
  graphPlotInt <- reactive({
    if(!is.null(markerTablePlotInt())){
      healthy_ct<-markerTablePlotInt()[DO_diseasetype=="healthy"]$celltype_ID
      if(length(healthy_ct>=1)){
        if(length(input$descendantsofInt)>=1){
          return(hierac_plot1_desc(markerTableCompleteInt(),PlotInt(),healthy_ct,descendantTableInt(),input$cellidInt, input$diseaseInt))
        }else{
          return(hierac_plot1(markerTableCompleteInt(),PlotInt(),healthy_ct,input$cellidInt, input$diseaseInt))
        }
      }
    } else{
      return(NULL)  # Return NULL if no plot is available
      
    }
  })
  
  output$plot1Int <- renderGrViz({
    graphPlotInt()  # Use the reactive graphPlot()
  })
  
  
  
  output$scalableplotInt <- renderUI({
    tagList( tags$div(
      #div(grVizOutput('plot1',height = input$height, width = input$width)))
      grVizOutput('plot1Int', height = "700px", width = "700px"),
      style = paste("transform: scale(", input$zoomInt, "); transform-origin: center;")
    )
    )
  })
  
  output$save_plotInt <- downloadHandler(
    filename = function() {
      if (input$file_formatInt == "PNG") {
        "Ontology_plot.png"
      } else if(input$file_formatInt == "PDF") {
        "Ontology_plot.pdf"
      }
    },
    content = function(file) {
      grViz_obj <- graphPlotInt()  # Use the stored reactive Graphviz plot
      
      if (is.null(grViz_obj)) {
        stop("No plot available to save.")
      }
      
      svg_file <- tempfile(fileext = ".svg")
      
      # Convert Graphviz to SVG and save
      svg_content <- export_svg(grViz_obj)  # Export Graphviz as SVG
      writeLines(svg_content, svg_file)  # Write to temporary SVG file
      
      # Get actual width from UI
      default <-700
      
      # Scale dimensions based on zoom input
      scaled_zoom <- default * input$zoomInt
      
      
      # Convert SVG to PNG with specified width & height
      if (input$file_format == "PNG") {
        suppressWarnings(
          rsvg_png(svg_file, file, width = scaled_zoom, height = scaled_zoom)
          
        )} else if(input$file_format == "PDF"){
          suppressWarnings(
            rsvg_pdf(svg_file, file, width = scaled_zoom / 100, height = scaled_zoom / 100) 
          )}
    }
  )
  
  
  click_plotInt<-reactive ({
    click_node(integratedDBmod(),ontology_celltype,PlotInt(),markerTableCompleteInt(),input$cellidInt, ontology_def, input$diseaseInt)
  })

  txtInt <- reactive({
    req(input$plot1Int_click)
    nodeval <- input$plot1Int_click$nodeValues[1:length(input$plot1Int_click$nodeValues)]
    nodeval<-as.data.table(paste(nodeval,sep="", collapse=" "))
    return(click_plotInt()[label %in% nodeval$V1]$cell_definition)

  })

  output$celltype_defInt <- renderText({
    req(txtInt())
    txtInt()
  })


  outputTableInt<-reactive({
    if(!is.null(markerTableOutputInt())){
      if("celltype" %in% colnames(markerTableOutputInt())){
        markerTableOutputInt()[DO_diseasetype=="healthy",CL_celltype:=celltype][DO_diseasetype=="healthy",CL_ID:=celltype_ID][DO_diseasetype=="healthy",CL_cell_definition:=cell_definition]
        markerTableOutputInt()[DO_diseasetype!="healthy",NCIT_celltype:=celltype][DO_diseasetype!="healthy",NCIT_ID:=celltype_ID][DO_diseasetype!="healthy",NCIT_cell_definition:=cell_definition]
        markerTableOutputInt()[,c("celltype","celltype_ID","cell_definition"):=NULL]
      }
      if("EC_score" %in% colnames(markerTableOutputInt())){
        setnames(markerTableOutputInt(), "EC_score", "ECs")
        setnames(markerTableOutputInt(), "specificity_score", "SPs")
      }

      if(length(input$celltypeInt)!=0 & input$tabletypeInt=="Complete" & input$tissue_awareInt==TRUE){
        if(input$mergeDescendantInt=="Yes"){
          markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                              c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","celltype_ancestor","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }else{
          markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                              c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_tissue","Uberon_tissue","Uberon_ID","tissue_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }
      } else if(length(input$celltypeInt)!=0 & input$tabletypeInt=="Complete" & input$tissue_awareInt==FALSE){
        if(input$mergeDescendantInt=="Yes"){
          markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                              c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","celltype_ancestor","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }else{
          markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                              c("species","original_diseasetype","DO_diseasetype","DO_ID","DO_definition","original_celltype","CL_celltype","CL_ID","CL_cell_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition","marker", "marker_type","gene_description","ECs","SPs","resource","log2FC","p.value","adjusted_p.value","pct1")]
        }
      }  else if(length(input$celltypeInt)!=0 & input$tabletypeInt=="Simple" & input$tissue_awareInt==TRUE){
        if(input$mergeDescendantInt=="Yes"){
          unique(markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                                     c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","celltype_ancestor","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }else{
          unique(markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                                     c("species","DO_diseasetype","DO_ID","Uberon_tissue","Uberon_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }
      } else if(length(input$celltypeInt)!=0 & input$tabletypeInt=="Simple" & input$tissue_awareInt==FALSE){
        if(input$mergeDescendantInt=="Yes"){
          unique(markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                                     c("species","DO_diseasetype","DO_ID","celltype_ancestor","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }else{
          unique(markerTableOutputInt()[ECs>= input$EC_scoreInt & SPs >= input$specificityInt][,
                                                                                                                                                     c("species","DO_diseasetype","DO_ID","CL_celltype","CL_ID","NCIT_celltype","NCIT_ID","marker", "marker_type","ECs","SPs")])
        }
      }
    }

  })

  outputTableInt2<- reactive({
    if(!is.null(outputTableInt())){
      if(all(is.na(outputTableInt()$DO_ID))){
        suppressWarnings({
          table<-outputTableInt()[,c("original_diseasetype","DO_diseasetype","DO_ID","DO_definition","NCIT_celltype","NCIT_ID","NCIT_cell_definition"):=NULL]
        })
      }
      if(all(is.na(outputTableInt()$CL_ID))){
        suppressWarnings({
          table<-outputTableInt()[,c("CL_celltype","CL_ID","CL_cell_definition"):=NULL]
        })
      }
      return(table)
    }
  })

  output$filteredInt <- renderDataTable({
    if(!is.null(outputTableInt())){
      datatable(outputTableInt2(),rownames = FALSE)

    }
  })


  # Downloadable csv of selected dataset

  output$downloadDataInt2 <- downloadHandler(
    filename = function() {
      paste0("TheIntegratedCellMarkerAccordionFilt_database_v0.9.5", input$downloadTypeInt2)
    },
    content = function(file) {
      if(input$downloadTypeInt2 == ".csv"){
        write.csv(outputTableInt2(), file, row.names = FALSE, quote=FALSE)
      }
      else if(input$downloadTypeInt2 == ".xlsx") {
        write_xlsx(outputTableInt2(), file)
      }
      else if(input$downloadTypeInt2 == ".tsv") {
        write.table(outputTableInt2(), file, quote = FALSE,
                    sep='\t', row.names = FALSE)
      }

    }

  )



  # server for enrichment -----
  observeEvent(input$annoInfo,{
    showModal(modalDialog(
      title = "Cell type marker enrichment analysis",
      HTML("
  Users can upload a file containing markers for each cluster or related to a single entity, and the 
  <strong>Cell Marker Accordion</strong> will retrieve the respective cell type with the highest correlation.<br>
  <br>
<strong>IMPORTANT! Enrichment analysis is performed using <em>Fisher's exact test</em> to identify significant associations between the input gene list and cell type-specific markers in the <strong>Cell Marker Accordion</strong> database. For full access to the Cell Marker Accordion annotation algorithm, 
we recommend using our <a href=\"https://github.com/TebaldiLab/cellmarkeraccordion\" target=\"_blank\">R package</a></strong>.<br>
<p>
<br>
  Click the  <strong>'InputFile'</strong> button to access details on the input file format. <br>
  Click the <strong>'Load demo example'</strong> button to load an example set of <strong>FindAllMarkers</strong> output. To download the example set, click the <strong> InputFile</strong> button and use the <strong>Download Example</strong> option. <br>
  The number of positive and negative genes to retain for each cluster can be specified by entering the desired number in the box and clicking the <strong>'Add value'</strong> button.<br>
</p>

<p>
  In the <strong>FindAllMarkers</strong> output, genes are classified as positive if <code>avg_log2FC > 0</code> and negative if <code>avg_log2FC < 0</code>.  Genes will be ranked based on their <code>avg_log2FC < 0</code> values and then only the top N for each cluster will be used. <br>
  In a custom table, this classification can be defined in the <code>'gene_type'</code> column.
</p>


"

      )
    ))
  })
  observeEvent(input$userclusterfileinfo,{
    showModal(modalDialog(
      title = "Input file format",
      HTML("Upload marker genes file to perform cell type marker enrichment analysis. <br>
    The file can be one of the following types:
    <ul><li> <strong> The Output table generated by the <em>FindAllMarkers </em> function in the Seurat package,</strong> as txt, xlsx, csv or tsv file</li>
           <strong> Example </strong>"),
    downloadButton("download_ex_anno1", "Download Example"),

     # HTML("<img src=findallmarker_ex.png>"),
     HTML("<table border='1' style='width:100%; border-collapse: collapse; text-align: center;'>
          <thead>
            <tr>
                        <tr style='background-color: #f2f2f2;'>

                <th style='text-align: center;'>cluster</th>
                <th style='text-align: center;'>p_val</th>
                <th style='text-align: center;'>avg_log2FC</th>
                <th style='text-align: center;'>pct.1</th>
                <th style='text-align: center;'>pct.2</th>
                <th style='text-align: center;'>p_val_adj</th>
                <th style='text-align: center;'>gene</th>
            </tr>
          </thead>
          <tbody>
<tr><td>0</td><td>0</td><td>1.196</td><td>0.327</td><td>0.135</td><td>0</td><td>CCR7</td></tr>
    <tr><td>0</td><td>0</td><td>1.178</td><td>0.320</td><td>0.122</td><td>0</td><td>LEF1</td></tr>
    <tr><td>0</td><td>0</td><td>0.999</td><td>0.447</td><td>0.275</td><td>0</td><td>SELL</td></tr>
    <tr><td>0</td><td>0</td><td>0.990</td><td>0.789</td><td>0.612</td><td>0</td><td>TMEM66</td></tr>
    <tr><td>0</td><td>0</td><td>0.985</td><td>0.544</td><td>0.306</td><td>0</td><td>CD27</td></tr>
    <tr><td>0</td><td>0</td><td>0.972</td><td>0.765</td><td>0.466</td><td>0</td><td>CD3E</td></tr>
    <tr><td>0</td><td>0</td><td>0.957</td><td>0.359</td><td>0.200</td><td>0</td><td>PIK3IP1</td></tr>
    <tr><td>0</td><td>0</td><td>0.948</td><td>0.620</td><td>0.411</td><td>0</td><td>GIMAP7</td></tr>
    <tr><td>0</td><td>0</td><td>0.921</td><td>0.275</td><td>0.133</td><td>0</td><td>MAL</td></tr>
    <tr><td>0</td><td>0</td><td>0.894</td><td>0.641</td><td>0.447</td><td>0</td><td>NOSIP</td></tr>
    <tr><td>1</td><td>0</td><td>2.879</td><td>0.727</td><td>0.037</td><td>0</td><td>CD8B</td></tr>
    <tr><td>1</td><td>0</td><td>2.114</td><td>0.390</td><td>0.009</td><td>0</td><td>RP11-291B21.2</td></tr>
    <tr><td>1</td><td>0</td><td>1.997</td><td>0.316</td><td>0.020</td><td>0</td><td>S100B</td></tr>
    <tr><td>1</td><td>0</td><td>1.568</td><td>0.339</td><td>0.032</td><td>0</td><td>CD8A</td></tr>
    <tr><td>1</td><td>0</td><td>0.905</td><td>0.513</td><td>0.255</td><td>0</td><td>RGS10</td></tr>
    <tr><td>1</td><td>0</td><td>0.633</td><td>0.808</td><td>0.460</td><td>0</td><td>CD3D</td></tr>
    <tr><td>1</td><td>0</td><td>0.623</td><td>0.699</td><td>0.438</td><td>0</td><td>NOSIP</td></tr>
    <tr><td>1</td><td>0</td><td>0.593</td><td>1.000</td><td>0.989</td><td>0</td><td>RPL31</td></tr>
    <tr><td>1</td><td>0</td><td>0.591</td><td>0.999</td><td>0.973</td><td>0</td><td>RPL34</td></tr>
    <tr><td>1</td><td>0</td><td>0.583</td><td>0.999</td><td>0.981</td><td>0</td><td>RPS25</td></tr>
          </tbody>
        </table>"),


      HTML("<br>"),
      HTML("<li> Alternatively, you can provide a <strong>custom table</strong> with at least one column (<strong>gene</strong> column) containing the list of genes (one per row).
      By default, all genes are considered positive (high expression) and associated with a single identity class (one cluster only). </li>
      \nYou may also include additional columns:
             <ul><li> <strong> cluster </strong>: indicate the identity class of the markers
              <li> <strong>  gene_type </strong>: positive, whether the gene is positive (high expression) or negative, whether the gene is negative (low expression)</ul>"),
      HTML("<br>"),
      HTML("<strong> Example </strong>"),
      downloadButton("download_ex_anno2", "Download Example"),

    HTML("<table border='1' style='width:100%; border-collapse: collapse; text-align: center;'>
            <thead>
              <tr style='background-color: #f2f2f2;'>
                <th style='text-align: center;'>cluster</th>
                <th style='text-align: center;'>gene</th>
                <th style='text-align: center;'>gene_type</th>
              </tr>
            </thead>
            <tbody>
              <tr><td>0</td><td>CCR7</td><td>positive</td></tr>
              <tr><td>0</td><td>TMEM66</td><td>positive</td></tr>
              <tr><td>0</td><td>HLA-DRA</td><td>negative</td></tr>
              <tr><td>0</td><td>CD74</td><td>negative</td></tr>
              <tr><td>1</td><td>S100B</td><td>positive</td></tr>
              <tr><td>1</td><td>RGS10</td><td>positive</td></tr>
              <tr><td>1</td><td>KLF6</td><td>negative</td></tr>
              <tr><td>1</td><td>CLIC1</td><td>negative</td></tr>
              <tr><td>2</td><td>VIM</td><td>positive</td></tr>
              <tr><td>2</td><td>GSTK1</td><td>positive</td></tr>
              <tr><td>2</td><td>CTSW</td><td>negative</td></tr>
              <tr><td>2</td><td>ZFA51</td><td>negative</td></tr>
            </tbody>
          </table>
        "),

      easyClose = TRUE,
      footer = NULL,
      size="l"
    ))
  })

  # Server logic for download
  output$download_ex_anno1 <- downloadHandler(
    filename = function() {
      "Custom_set_annotation_example1.xlsx"
    },
    content = function(file) {
      table<-fread("data/FindAllMarkers_small_ex.tsv")
      write_xlsx(table, file)
    }
  )

  # Server logic for download
  output$download_ex_anno2 <- downloadHandler(
    filename = function() {
      "Custom_set_annotation_example2.xlsx"
    },
    content = function(file) {
      table<-read_xlsx("data/Demo_example_annotation2.xlsx")
      write_xlsx(table, file)
    }
  )


  observeEvent(input$filterhelpA,{
    showModal(modalDialog(
      title = "Filters Information",
      HTML("Once the custom file is loaded, you can specify filters for the <strong>Cell Marker Accordion</strong> database before running enrichment. In particular:
</p>

<ul>
  <li><strong>Select species:</strong> currently Human and/or Mouse.</li>
  <li><strong>Condition:</strong> healthy or multiple diseases.</li>
  <li>
    <strong>Tissue:</strong> select one or multiple tissues from the list. If no tissue is selected, the enrichment will be performed using all tissues combined. 
    If multiple tissues are selected, the enrichment will be performed by combining the selected tissues.
  </li>
  <li>
    <strong>Cell type:</strong> select one or multiple cell types from the list. If no cell type is selected, the enrichment will be performed using all cell types.
  </li>
  <li>
    <strong>ECs (evidence consistency score):</strong> measures the agreement of different enrichment sources. 
    Filter marker genes from the <strong>Cell Marker Accordion</strong> database with an ECs &ge; the selected value.
  </li>
  <li>
    <strong>SPs (specificity score):</strong> indicates whether a gene is a marker for different cell types present in the entire <strong>Accordion</strong> database. 
    Filter marker genes from the <strong>Accordion</strong> database with an SPs &ge; the selected value.
  </li>
  <li>
    <strong>Maximum number of markers to keep for each cell type:</strong> specify the top N marker genes to retain for each cell type during automatic enrichment 
    Markers are ranked based on their ECs and SPs. Default is <code>50</code>.
  </li>
</ul>")    ))
  })

  observeEvent(input$addpos, {
    req(input$nmarkerotherpos)
    otherVal <- "nmarkerotherpos"
    names(otherVal) <- input$nmarkerotherpos
    updatedValues <- c(nmarker_values, otherVal)
    updateRadioButtons(session, "nmarkerpos", choices = updatedValues)
  })

  observeEvent(input$addneg, {
    req(input$nmarkerotherneg)
    otherVal <- "nmarkerotherneg"
    names(otherVal) <- input$nmarkerotherneg
    updatedValues <- c(nmarker_values, otherVal)
    updateRadioButtons(session, "nmarkerneg", choices = updatedValues)
  })

  inputTable <- reactive({
    suppressWarnings({
    if(input$button){
    if (is.null(input$clusterfile) & input$demo_ex_anno > 0) {
      # Read the default example file
      table <- fread("data/FindAllMarkers_small_ex.tsv", verbose = FALSE)

      if (!is.na(as.numeric(input$nmarkerpos))) {
        positive_marker <- table %>%
          group_by(cluster) %>%
          slice_max(n = as.numeric(input$nmarkerpos), order_by = avg_log2FC)
        positive_marker <- as.data.table(positive_marker)
        positive_marker <- positive_marker[avg_log2FC > 0]
      } else { # ALL
        positive_marker <- as.data.table(table)
        positive_marker <- positive_marker[avg_log2FC > 0]
      }

      positive_marker<-positive_marker[gene %in% (accordion_complete[species %in% input$speciesA]$marker)]
      positive_marker[, pos_marker := paste(gene, collapse = " "), by = cluster]
      positive_marker <- unique(positive_marker[, c("cluster", "pos_marker")])

      if (!is.na(as.numeric(input$nmarkerneg))) { # Select number of markers
        negative_marker <- table %>%
          group_by(cluster) %>%
          slice_max(n = as.numeric(input$nmarkerneg), order_by = avg_log2FC)
        negative_marker <- as.data.table(negative_marker)
        negative_marker <- negative_marker[avg_log2FC < 0]
      } else { # ALL
        negative_marker <- as.data.table(table)
        negative_marker <- negative_marker[avg_log2FC < 0]
      }

      negative_marker<-negative_marker[gene %in% (accordion_complete[species %in% input$speciesA]$marker)]
      negative_marker[, neg_marker := paste(gene, collapse = " "), by = cluster]
      negative_marker <- unique(negative_marker[, c("cluster", "neg_marker")])

      input_marker <- merge(positive_marker, negative_marker, by = "cluster", all.x = TRUE, all.y = TRUE)
      return(input_marker)


    } else if (!is.null(input$clusterfile)) {
      file_load_cluster <- input$clusterfile
      fileName_cluster <- file_load_cluster$datapath
      ext <- tools::file_ext(fileName_cluster)

      req(file_load_cluster)
      validate(need(ext %in% c("csv", "xlsx", "txt", "tsv"), "Please upload a csv, txt, tsv or xlsx file"))

      if (ext == "xlsx") {
        user_inputfile <- read_excel(fileName_cluster)
      } else if (ext %in% c("csv", "txt", "tsv")) {
        user_inputfile <- fread(fileName_cluster)
      }

      user_inputfile <- as.data.table(user_inputfile)

      col_names_findmarkers <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene") # Expected column names
      if (all(col_names_findmarkers %in% colnames(user_inputfile))) {  # If input matches FindMarker output
        user_inputfile<-user_inputfile[,c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")]

        if (!is.na(as.numeric(input$nmarkerpos))) {
          positive_marker <- user_inputfile %>%
            group_by(cluster) %>%
            slice_max(n = as.numeric(input$nmarkerpos), order_by = avg_log2FC)
          positive_marker <- as.data.table(positive_marker)
          positive_marker <- positive_marker[avg_log2FC > 0]
        } else { # ALL
          positive_marker <- as.data.table(user_inputfile)
          positive_marker <- positive_marker[avg_log2FC > 0]
        }

        positive_marker<-positive_marker[gene %in% (accordion_complete[species %in% input$speciesA]$marker)]
        positive_marker[, pos_marker := paste(gene, collapse = " "), by = cluster]
        positive_marker <- unique(positive_marker[, c("cluster", "pos_marker")])

        if (!is.na(as.numeric(input$nmarkerneg))) { # Select number of markers
          negative_marker <- user_inputfile %>%
            group_by(cluster) %>%
            slice_max(n = as.numeric(input$nmarkerneg), order_by = -avg_log2FC)
          negative_marker <- as.data.table(negative_marker)
          negative_marker <- negative_marker[avg_log2FC < 0]
        } else { # ALL
          negative_marker <- as.data.table(user_inputfile)
          negative_marker <- negative_marker[avg_log2FC < 0]
        }

        negative_marker<-negative_marker[gene %in% (accordion_complete[species %in% input$speciesA]$marker)]
        negative_marker[, neg_marker := paste(gene, collapse = " "), by = cluster]
        negative_marker <- unique(negative_marker[, c("cluster", "neg_marker")])

        input_marker <- merge(positive_marker, negative_marker, by = "cluster", all.x = TRUE, all.y = TRUE)
        return(input_marker)

      } else { # Input different from FindMarkers output
        user_inputfile <- as.data.table(user_inputfile)

        if (!"gene" %in% colnames(user_inputfile)) {
          validate(need("gene" %in% colnames(user_inputfile), "\"gene\" column not found. Please specify the \"gene\" column containing the list of genes for enrichment analysis."))
        }

        if (!"gene_type" %in% colnames(user_inputfile)) {  # Corrected 'colnmames' typo
          user_inputfile[, gene_type := "positive"]
          showNotification("\"gene_type\" column not found. Genes will be considered as positive", type = "warning", duration = 5)
        }

        if (!"cluster" %in% colnames(user_inputfile)) {
          user_inputfile[, cluster := "custom_cluster"]
          showNotification("\"cluster\" column not found. Genes will be associated with a single identity class", type = "warning", duration = 5)
        }

        pos_marker <- user_inputfile[gene_type == "positive"][, pos_marker := paste(gene, collapse = " "), by = cluster]
        pos_marker <- unique(pos_marker[, c("cluster", "pos_marker")])

        neg_marker <- user_inputfile[gene_type == "negative"][, neg_marker := paste(gene, collapse = " "), by = cluster]
        neg_marker <- unique(neg_marker[, c("cluster", "neg_marker")])

        user_inputfile <- merge(pos_marker, neg_marker, by = "cluster", all.x = TRUE, all.y = TRUE)
        return(user_inputfile)
      }
    }
    }
      })
  })

  
  demoLoaded_anno <- reactiveVal(FALSE)  # Track if demo is loaded
  
  observeEvent(input$demo_ex_anno, {
    showNotification("Demo example loaded!", type = "message", duration = 5)
    demoLoaded_anno(TRUE)  # Mark demo as loaded
  })
  
  output$success_icon_anno <- renderUI({
    if (demoLoaded_anno()) {
      tags$span(
        icon("check-circle", class = "text-success", lib = "font-awesome"),
        style = "font-size: 20px; margin-left: 10px;"
      )
    }
  })
  

  output$exampletable <- renderDataTable({
    if(input$demo_ex_anno>0 & is.null(input$clusterfile)){
      table<-fread("data/FindAllMarkers_small_ex.tsv", verbose = F)
      table[,p_val:=format(p_val, digits=2)]
      table<-table[order(-avg_log2FC)]
      table[,avg_log2FC:=format(avg_log2FC, digits=2)]
      table<-datatable(table, options = list(scrollX='400px',
        dom = 't',      # Removes unnecessary UI elements
        pageLength = 5 # Controls number of rows displayed
      ),rownames = FALSE) %>%
        formatStyle(columns = names(table), fontSize = '12px')%>%
        formatStyle(columns = names(table), target="row", fontSize = '12px')
      return(table)
    }

  })


  markerTable<-reactive ({
    table<-accordion_complete[DO_diseasetype %in% input$diseaseA & species %in% input$speciesA]
    table
  })


  toListen_tissueA <- reactive({
    list(input$speciesA, input$diseaseA)
  })

  observeEvent(toListen_tissueA(),{
    updatePickerInput(session,'tissueA', selected = "blood",
                      choices=unique(markerTable()$Uberon_tissue),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(markerTable()$Uberon_tissue))))
  })

  toListen_celltypeA <- reactive({
    list(input$speciesA, input$diseaseA, input$tissueA)
  })
  
  observeEvent(toListen_celltypeA(),{
    updatePickerInput(session,'celltypeA', selected = unique(markerTable()$celltype),
                      choices=unique(markerTable()$celltype),
                      option=list(`actions-box` = TRUE,style="box-celltypes"),
                      choicesOpt = list(style = rep(("font-size: 18px; line-height: 1.6;"), uniqueN(markerTable()$celltype))))
  })

  markerTableTissue <- reactive ({
    if(is.null(input$tissueA)){ #if no tissues are selected perform annotation with all aggregated tissues
      markerTable()[, Uberon_ID:="ALL"]
      markerTable()[, Uberon_tissue:="ALL"]
      table<- markerTable()
    } else{
      table<-markerTable()[Uberon_tissue %in% input$tissueA]
    }
    tissue<-unique(table$Uberon_tissue) 
    if(length(tissue) >0){ #if more tissues are selected perform annotation with all aggregated tissues
      table[,Uberon_tissue:=paste(tissue,collapse=", ")]
    }
    if(!is.null(input$celltypeA)){ #if no cell types are selected consider all, otherwise consider only the selected cell type
      table<-table[celltype %in% input$celltypeA]
    }
    
      #calculate EC score based on filtered tissue
      accordion_marker<-unique(table[,c("species","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker","marker_type","resource")])
      ECs<-ddply(accordion_marker,.(species,Uberon_tissue,Uberon_ID,celltype,celltype_ID,marker,marker_type),nrow)
      colnames(ECs)<-c("species","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker", "marker_type","ECs")
      accordion_marker<-merge(accordion_marker,ECs,by=c("species","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker","marker_type"), all.x = TRUE)
      accordion_marker[,ECs_reg := log10(ECs)+1]

      #compute SPs for positive and negative markers
      mark_spec<-ddply(table,.(marker,marker_type),nrow)
      colnames(mark_spec)<-c("marker","marker_type","SPs")
      accordion_marker<-merge(accordion_marker,mark_spec,by=c("marker","marker_type"),all.x = TRUE)

      length_ct_pos<-uniqueN(accordion_marker[marker_type=="positive"]$celltype)
      length_ct_neg<-uniqueN(accordion_marker[marker_type=="negative"]$celltype)

      #scale and log transforme SPs
      accordion_marker<-accordion_marker[marker_type=="positive",SPs_reg := scales::rescale(as.numeric(SPs), to = c(1,length_ct_pos),from = c(length_ct_pos,1))
      ][marker_type=="negative",SPs_reg := scales::rescale(as.numeric(SPs), to = c(1,length_ct_neg),from = c(length_ct_neg,1))
      ][,c("species","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker","marker_type","SPs","SPs_reg","ECs_reg","ECs")]

      accordion_marker<-accordion_marker[marker_type=="positive",SPs_reg := scales::rescale(as.numeric(SPs_reg), to = c(min(accordion_marker[marker_type=="positive"]$ECs),max(accordion_marker[marker_type=="positive"]$ECs)))
      ][marker_type=="negative",SPs_reg := scales::rescale(as.numeric(SPs_reg), to = c(min(accordion_marker[marker_type=="negative"]$ECs),max(accordion_marker[marker_type=="negative"]$ECs)))
      ][,c("species","Uberon_tissue","Uberon_ID","celltype","celltype_ID","marker","marker_type","SPs","SPs_reg","ECs_reg","ECs")]

      accordion_marker[,SPs:=as.numeric(format(round(1/SPs,2), nsmall=2))]
      accordion_marker[,SPs_reg:=log10(SPs_reg)+1]

      accordion_marker[,combined_score := SPs_reg * ECs_reg]
      accordion_marker<-unique(accordion_marker)
      accordion_marker

  })


  accMarkerAnnoTable <- reactive ({
    accordion_complete_filt<-markerTableTissue()[ECs >= input$EC_scoreA & SPs >= input$specificityA]
    #req(accordion_complete_filt)
    validate(need(nrow(accordion_complete_filt)>0, "No marker genes found with the selected filters."))   
    if(nrow(accordion_complete_filt)>0){
      if(as.numeric(input$max_n_marker)){
        accordion_complete_filt<-accordion_complete_filt[order(-combined_score)]
        accordion_complete_filt<-accordion_complete_filt[,.SD[1:input$max_n_marker], by=c("celltype")]
      }
  
      positive_marker_acc<-accordion_complete_filt[marker_type=="positive" & species %in% input$speciesA][,pos_marker_acc:= paste(marker, collapse=" "), by= celltype]
      positive_marker_acc<-unique(positive_marker_acc[,c("celltype","pos_marker_acc")])
      negative_marker_acc<-accordion_complete_filt[marker_type=="negative" & species %in% input$speciesA][,neg_marker_acc:= paste(marker, collapse=" "), by= celltype]
      negative_marker_acc<-unique(negative_marker_acc[,c("celltype","neg_marker_acc")])
  
      table<-merge(positive_marker_acc,negative_marker_acc, by="celltype",all.x=TRUE, all.y=TRUE)
      as.data.table(table)
    }
  })


  inputTableLong <- reactive ({
    if(!is.null(accMarkerAnnoTable())){
    suppressWarnings({

    input_marker_long<-do.call("rbind", replicate(nrow(accMarkerAnnoTable()), inputTable(), simplify = FALSE))
    input_marker_long<-input_marker_long[order(cluster)]
    input_marker_long
    })
    }
  })
  accMarkerAnnoTableLong <- reactive ({
    if(!is.null(accMarkerAnnoTable())){
    suppressWarnings({
    do.call("rbind", replicate(nrow(inputTable()), accMarkerAnnoTable(), simplify = FALSE))
    })
    }
  })
  combineTable <- reactive({ #create table to perform Fisher Test
    if(!is.null(accMarkerAnnoTableLong())){
      
    dt_combine<-cbind(inputTableLong(),accMarkerAnnoTableLong())
    dt_combine<-as.data.table(dt_combine)
    #IY overlap between input and anno
    dt_combine[,positive_intersect:= mapply(function(x, y) paste0(str_sort(intersect(x, y)), collapse = " "), strsplit(dt_combine$pos_marker, ' '), strsplit(dt_combine$pos_marker_acc, ' '))]
    dt_combine[,negative_intersect:= mapply(function(x, y) paste0(str_sort(intersect(x, y)), collapse = " "), strsplit(dt_combine$neg_marker, ' '), strsplit(dt_combine$neg_marker_acc, ' '))]

    dt_combine[,tot_overlap:= paste0(positive_intersect, " ", negative_intersect)]
    dt_combine[,tot_overlap := str_trim(tot_overlap)]

    #IN  in input, not in anno
    dt_combine[,in_input_not_anno_pos := mapply(function(x, y) paste0(setdiff(x, y), collapse = " "), strsplit(dt_combine$pos_marker, ' '), strsplit(dt_combine$pos_marker_acc, ' '))]
    dt_combine[,in_input_not_anno_neg:= mapply(function(x, y) paste0(setdiff(x, y), collapse = " "), strsplit(dt_combine$neg_marker, ' '), strsplit(dt_combine$neg_marker_acc, ' '))]
    dt_combine[,in_input_not_anno_tot:= paste0(in_input_not_anno_pos, " ", in_input_not_anno_neg)]
    dt_combine[,in_input_not_anno_tot := str_trim(in_input_not_anno_tot)]

    # BY not in input, in anno
    dt_combine[,not_input_in_anno_pos := mapply(function(x, y) paste0(setdiff(x, y), collapse = " "), strsplit(dt_combine$pos_marker_acc, ' '), strsplit(dt_combine$pos_marker, ' '))]
    dt_combine[,not_input_in_anno_neg:= mapply(function(x, y) paste0(setdiff(x, y), collapse = " "), strsplit(dt_combine$neg_marker_acc, ' '), strsplit(dt_combine$neg_marker, ' '))]
    dt_combine[,not_input_in_anno_tot:= paste0(not_input_in_anno_pos, " ", not_input_in_anno_neg)]
    dt_combine[,not_input_in_anno_tot := str_trim(not_input_in_anno_tot)]
    dt_combine
    }
  })

  toListen_anno <- reactive({
    list(input$button, inputTable(),accMarkerAnnoTable())
  })

  observeEvent(toListen_anno, {
    annoResultsTable <- reactive({
      if(input$button >0){
        if(!is.null(inputTable())){
        if(nrow(inputTable())){

        out_df<-tibble()
        inputseppos<-separate_rows(as.data.frame(inputTable()), pos_marker, sep=" ")
        inputsepneg<-separate_rows(as.data.frame(inputTable()), neg_marker, sep=" ")
        univ<- uniqueN(unique(c(accordion_complete$marker, inputseppos$pos_marker, inputsepneg$neg_marker)))
        withProgress(message = 'Running enrichment analysis...', value = 0, {
          for (cl in unique(inputTable()$cluster)){
            sub_dt<-combineTable()[cluster==cl]
            # Increment the progress bar, and update the detail text.
            incProgress(1 / (nrow(inputTable())-1), detail = paste("Analyzing cluster", cl))
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            for (ct in (accMarkerAnnoTable()$celltype)){
              #print(ct)
              sub_ct_dt<-sub_dt[celltype==ct]
              IY<- count.fields(textConnection(sub_ct_dt$tot_overlap), sep = " ") # in input, in anno (overlap)
              IN<- count.fields(textConnection(sub_ct_dt$in_input_not_anno_tot), sep = " ") # in input, not in anno
              BY<- count.fields(textConnection(sub_ct_dt$not_input_in_anno_tot), sep = " ") # not in input, in anno
              BN<- univ - (IY+IN+BY) # not in input, not in anno
              fisher_data <- fisher.test(matrix(c(IY, BY, IN, BN), 2, 2), alternative="greater")
              # generation of the output (a tibble with 1 row)
              out_df <- rbind(out_df, tibble("cluster"=cl, # name of the input
                                             "celltype"=ct, # name of the annotation
                                             #"anno_class"=anno_class,
                                             "overlap_size"=IY, # size of the overlap
                                             "p_value"=fisher_data$p.value, # p-value of the enrichment (Fisher test)
                                             "odds_ratio"=fisher_data$estimate, # conditional MLE estimate from fisher.test function
                                             #"combined_score"=0, # default combined score (to match enrichr enrichment results)
                                             #"odds_ratio_sample"= (IY/IN)/(BY/BN), # sample odds ratio
                                             #"fold_enrichment"= (IY/(IY+IN))/((IY+BY)/(BY+IY+BY+BN)), # fold enrichment formula
                                             "input_size"=IY+IN, # size of processed input
                                             "anno_size"=IY+BY, # size of processed annotation
                                             "background_size"= IY+IN+BY+BN, # size of processed background
                                             "overlap_input_ratio"= IY/(IY+IN), # overlap ratio with respect to input size
                                             "overlap_anno_ratio"= IY/(IY+BY), # overlap ratio with respect to annotation size
                                             "positive_overlap_ids"= sub_ct_dt$positive_intersect, # string with positive overlap markers ids
                                             "negative_overlap_ids" = sub_ct_dt$negative_intersect, # string with negative overlap markers ids
                                             #"overlap_ids"=  sub_ct_dt$tot_overlap # string with overlap ids
              ))
            }
          }
          out_dt<-as.data.table(out_df)
          out_dt<-out_dt[order(cluster,p_value)]
          out_dt_max<-out_dt[!duplicated(out_dt$cluster)]
          out_dt_max[,odds_ratio:=format(odds_ratio, digits=2)]
          out_dt_max[,p_value:=format(p_value, digits=2)]
          out_dt_max[,overlap_input_ratio:=format(overlap_input_ratio, digits=2)]
          out_dt_max[,overlap_anno_ratio:=format(overlap_anno_ratio, digits=2)]

          out_dt_max

        })
      }
        }
      }
    })

    output$table1A <- renderDataTable({
      datatable(annoResultsTable(),rownames = FALSE,options = list(
        pageLength=5, scrollX='400px',columnDefs = list(
          list(width = '130px', targets = "_all"))), filter = 'top')
    })

    output$downloadDataA <- downloadHandler(
      filename = function() {
        paste0("TheCellMarkerAccordion_enrichment_results", input$downloadTypeA)
      },
      content = function(file) {
        if(input$downloadTypeA == ".csv"){
          write.csv(annoResultsTable(), file, row.names = FALSE, quote=FALSE)
        }
        else if(input$downloadTypeA == ".xlsx") {
          write_xlsx(annoResultsTable(), file)
        }
        else if(input$downloadTypeA == ".tsv") {
          write.table(annoResultsTable(), file, quote = FALSE,
                      sep='\t', row.names = FALSE)
        }

      }


    )


    tablePlotOntology <- reactive ({
      if(!is.null(annoResultsTable())){
        if(nrow(annoResultsTable())>0){
      dt_plot<-merge(annoResultsTable(),ontology_celltype, by="celltype")
      dt_plot<-dt_plot[,c("celltype","celltype_ID")]
      dt_plot
        }
      }
    })

    PlotA <- reactive({
      if (nrow(tableInputComplete())>=1) {
        cell_matching<-unique(tablePlotOntology()$celltype_ID)
        if(length(cell_matching)>1){
          ontosubplot<-onto_plot2(cell_onto,cell_matching ,cex=0.8)
          nodes<-as.data.table(ontosubplot@nodes)
          nodes<-nodes[,V1:=tstrsplit(nodes$V1,"CL", fixed = TRUE, keep = 1)]
          nodes<-nodes[, V1:=str_replace_all(V1,"\n"," ")]
          nodes<-nodes[, V1:=str_sub(V1,1,nchar(V1)-1)]
          ontosubplot@nodes<-nodes$V1
          ontosubplot
        }
        else{
          ontosubplot<-onto_plot(cell_onto,cell_matching ,cex=0.8)
          ontosubplot2<-make_graphNEL_from_ontology_plot(ontosubplot)
          dt_onto2<-as.data.table(ontosubplot2@nodes)
          label<-merge(dt_onto2,ontology_celltype,by.x="V1",by.y="celltype_ID")
          ontosubplot[["node_attributes"]][["label"]][[label$V1]]<-label$celltype
          ontosubplot2@nodes<-ontosubplot[["node_attributes"]][["label"]]
          ontosubplot2

        }

      }

    })


    
    # Reactive expression to store the Graphviz plot
    graphPlotA <- reactive({
      if(!is.null(tablePlotOntology())){
        if(nrow(tablePlotOntology())>0){
          return(hierac_plot2(accordion_complete, PlotA(), tablePlotOntology(),input$cellidA, input$diseaseA))
        }
      } else{
        return(NULL)  # Return NULL if no plot is available
        
      }
    })
    
    output$plot1A <- renderGrViz({
      graphPlotA()  # Use the reactive graphPlot()
    })
    
    
    output$scalableplotA <- renderUI({
      tagList( tags$div(
        #div(grVizOutput('plot1',height = input$height, width = input$width)))
        grVizOutput('plot1A', height = "700px", width = "700px"),
        style = paste("transform: scale(", input$zoomA, "); transform-origin: center;")
      )
      )
    })
    
    output$save_plotA <- downloadHandler(
      filename = function() {
        if (input$file_formatA == "PNG") {
          "Ontology_plot.png"
        } else if(input$file_formatA == "PDF") {
          "Ontology_plot.pdf"
        }
      },
      content = function(file) {
        grViz_obj <- graphPlotA()  # Use the stored reactive Graphviz plot
        
        if (is.null(grViz_obj)) {
          stop("No plot available to save.")
        }
        
        svg_file <- tempfile(fileext = ".svg")
        
        # Convert Graphviz to SVG and save
        svg_content <- export_svg(grViz_obj)  # Export Graphviz as SVG
        writeLines(svg_content, svg_file)  # Write to temporary SVG file
        
        # Get actual width from UI
        default <-700
        
        # Scale dimensions based on zoom input
        scaled_zoom <- default * input$zoomA
        
        
        # Convert SVG to PNG with specified width & height
        if (input$file_format == "PNG") {
          suppressWarnings(
            rsvg_png(svg_file, file, width = scaled_zoom, height = scaled_zoom)
            
          )} else if(input$file_format == "PDF"){
            suppressWarnings(
              rsvg_pdf(svg_file, file, width = scaled_zoom / 100, height = scaled_zoom / 100) 
            )}
      }
    )
    




    click_plotA<-reactive ({
      click_node(accordion_complete, ontology_celltype, PlotA(), tablePlotOntology(),input$cellidA, ontology_def, input$diseaseA)
    })


    txtA <- reactive({
      req(input$plot1A_click)
      nodeval <- input$plot1A_click$nodeValues[1:length(input$plot1A_click$nodeValues)]
      nodeval<-as.data.table(paste(nodeval,sep="", collapse=" "))
      return(  click_plotA()[label %in% nodeval$V1]$cell_definition)

    })
    output$celltype_defA <- renderText({
      req(txtA())
      txtA()
    })

  })



}



shinyApp(ui = ui, server = server, options=list(port=8200))

