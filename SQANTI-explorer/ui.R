#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(plotly)
library(igvShiny)
library(htmlwidgets)

# Define UI for application that draws a histogram
shinyUI(dashboardPage(skin="red",

    # Application title
    dashboardHeader(title="SQANTI Explorer"),

    # Sidebar with a slider input for number of bins
    dashboardSidebar(width=350,
        sidebarMenu(
            menuItem("Data", tabName="data", icon=icon("database")),
            menuItem("Plots", tabName="plots", icon=icon("chart-bar")),
            menuItem("Genome Browser", tabName = "browser", icon=icon("dna")),
            hr(),
            h4("Filter data: ", style = "margin-left: 1.3rem"),
            checkboxInput("polyexonic", "Polyexonic"),
            checkboxInput("monoexonic", "Monoexonic"),
            checkboxInput("noRTS", "No RTS"),
            checkboxInput("noIntraPriming", "No intra-priming"),
            checkboxInput("allCanonical", "All Canonical SJs"),
            checkboxInput("minCovGTZero", "min_cov > 0"),
            checkboxInput("onlyGenes", "Only Genes (does not accumulate values from transcripts)")
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName="data",
                fluidRow(
                    box(
                        title="Add Classification Files",
                        fileInput("classification_file", "Classification File: ", accept = c(".txt")),
                        fileInput("gtf_file", "GTF (GFF v2) File: ", accept = c(".gtf")),
                        textInput("name", "Name: "),
                        hr(),
                        h4("Genome (only required for genome browser)"),
                        selectInput("genome", label="Genome: ", choices=igvShiny::getSupportedGenomes()),
                        actionButton("addClassification", label = "Add Classification File")
                    ),
                    box(title = "Your Inputs", tableOutput('inputTable'), downloadButton("downloadData", "Download Data")),
                )
            ),
            tabItem(tabName="plots",
                fluidRow(
                    column(width=4,
                        box(width=NULL,
                            title="Group By",
                            selectInput("groupBy", "Group By:", list(
                                "SQANTI Filter" = "SQANTI_filter", "Novel Transcript" = "novel_transcript", "Novel Gene" = "novel_gene", "All (no grouping)" = "name")
                            )
                        ),
                        valueBoxOutput(width=NULL, "selected_transcript_count"),
                        box(width = NULL, title="Download",
                            downloadButton("downloadFilteredData", "Download Filtered Classification", style = "margin-bottom: 1rem"),
                            downloadButton("downloadSelectedData", "Download Selected Classification", style = "margin-bottom: 1rem"),
                            downloadButton("downloadFilteredGTF", "Download Filtered GTF", style = "margin-bottom: 1rem"),
                            downloadButton("downloadSelectedGTF", "Download Selected GTF", style = "margin-bottom: 1rem")
                        )
                    ),
                    column(width=8,
                        tabBox(width=NULL,title = "Plot",height = 500,
                               tabPanel("Summary", plotlyOutput("pie_chart")),
                               tabPanel("Basic Count", plotlyOutput("count_plot")),
                               tabPanel("Basic %", plotlyOutput("perc_plot")),
                               tabPanel("% Monoexonic", plotlyOutput("mono_plot")),
                               tabPanel("% Artifacts", plotlyOutput("arti_plot")),
                               tabPanel("% Novel Transcipts", plotlyOutput("novel_trans_plot")),
                               tabPanel("% Novel Genes", plotlyOutput("novel_genes_plot")),
                               tabPanel("Gene Expression", plotOutput("gene_expression"))
                        )
                    )
                )
            ),
            tabItem(tabName="browser",
                fluidRow(box(width=12,
                             uiOutput("selectGenomeData"),
                             actionButton("render_igv", "Render")
                )),
                fluidRow(box(width=12,title = "Genome Browser", igvShinyOutput("igv")))
            )
        ),
    )
))
