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

# Define UI for application that draws a histogram
shinyUI(dashboardPage(skin="red",

    # Application title
    dashboardHeader(title="SQANTI Explorer"),

    # Sidebar with a slider input for number of bins
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data", tabName="data", icon=icon("database")),
            menuItem("Plots", tabName="plots", icon=icon("chart-bar")),
            menuItem("Genome Browser", tabName = "browser", icon=icon("dna"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName="data",
                fluidRow(
                    box(
                        title="Add Classification Files",
                        fileInput("classification_file", "Classification File: ", accept = c(".txt")),
                        textInput("name", "Name: "),
                        actionButton("addClassification", label = "Add Classification File"),
                        hr(),
                        fileInput("BED_file", "BED File: ", accept = c(".bed"))
                    ),
                    box(title = "Your Inputs", tableOutput('inputTable'), downloadButton("downloadData", "Download Data")),
                )
            ),
            tabItem(tabName="plots",
                fluidRow(
                    column(width=4,
                        box(width=NULL,
                            title="Filters",
                            checkboxInput("polyexonic", "Polyexonic"),
                            checkboxInput("monoexonic", "Monoexonic"),
                            checkboxInput("noRTS", "No RTS"),
                            checkboxInput("noIntraPriming", "No intra-priming"),
                            checkboxInput("allCanonical", "All Canonical SJs"),
                            checkboxInput("minCovNotNa", "No NA in min_cov"),
                            checkboxInput("minCovGTZero", "min_cov > 0"),
                            checkboxInput("onlyGenes", "Only Genes (does not accumulate values from transcripts)"),
                        ),
                        box(width=NULL,
                            title="Group By",
                            selectInput("groupBy", "Group By:", list(
                                "SQANTI Filter" = "SQANTI_filter", "Novel Transcript" = "novel_transcript", "Novel Gene" = "novel_gene", "All (no grouping)" = "name")
                            ),
                        )
                    ),
                    column(width=8,
                        tabBox(width=NULL,title = "Plot", 
                               tabPanel("Basic Count", plotlyOutput("count_plot")),
                               tabPanel("Basic %", plotlyOutput("perc_plot")),
                               tabPanel("% Monoexonic", plotlyOutput("mono_plot")),
                               tabPanel("% Artifacts", plotlyOutput("arti_plot")),
                               tabPanel("% Novel Transcipts", plotlyOutput("novel_trans_plot")),
                               tabPanel("% Novel Genes", plotlyOutput("novel_genes_plot"))
                        ),
                        valueBoxOutput("selected_transcript_count"),
                        box(width = NULL, title="Selected Transcripts",
                            textOutput("selected_transcripts")    
                        ),
                        box(width = NULL, title="Download",
                            downloadButton("downloadFilteredData", "Download Filtered Data"),
                            downloadButton("downloadSelectedData", "Download Selected Data")
                        )
                    )
                )
            ),
            tabItem(tabName="browser",
                fluidRow(box(width=12,title = "Genome", uiOutput("epivizChart")))
            )
        ),
    )
))