#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

require(readr) # Read tables
require(purrr) # FP
require(tibble) # Better dfs
require(fs)
require(stringr)
library(magrittr)
require(tidyr)
require(ggplot2)
require(dplyr)
require(epivizrChart)
library(plotly)

options(shiny.maxRequestSize = 30*1024^2)

# Plot helpers
geom_perc_y <- scale_y_continuous(labels = scales::percent_format(), name="Percentage (%)")
geom_label_stacked = geom_text(size = 3, position = position_stack(vjust = 0.5))
geom_label_dodge = geom_text(size=3, position=position_dodge(width=0.9), vjust=1.5)

bar_plot = function(data, input) {
    ggplot(data, aes_string(x="name", y="n", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity")
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    classifications <- reactiveVal()
    makeReactiveBinding("classifications")
    
    observeEvent(input$addClassification, {
        req(input$classification_file, input$name)
        
        if(input$addClassification > 0) {
            isolate({
                datapath <- input$classification_file[["datapath"]]
                new_row <- tibble(
                    name = input$name,
                    file = input$classification_file[["name"]],
                    path = datapath,
                    classification = list(read_tsv(datapath))
                ) %>%
                    unnest(cols=c(classification)) %>%
                    mutate(polyexonic = if_else(exons > 1, "Polyexonic", "Monoexonic")) %>%
                    mutate(novel_transcript = if_else(associated_transcript == "novel", "Novel", "Annotated")) %>%
                    mutate(novel_gene = if_else(grepl("novelGene", associated_gene), "Novel", "Annotated")) %>%
                    mutate(log_gene_exp = log(gene_exp + 0.01))

                if (is.null(classifications())) {
                    classifications(new_row)
                }
                else {
                    classifications(classifications() %>% add_row(new_row))
                }
            })
        }
    })
    
    data_to_plot <- reactive({
        validate(
            need(classifications(), "Please add a classification file.")
        )
        val <- classifications() %>% group_by_at(c("name", input$groupBy))

        if (input$polyexonic) {
            val <- filter(val, exons > 1)
        }
        
        if (input$monoexonic) {
            val <- filter(val, exons == 1)
        }
        
        if(input$noRTS) {
            val <- filter(val, RTS_stage == FALSE)
        }
        
        if(input$noIntraPriming) {
            val <- filter(val, `intra-priming` == FALSE)
        }
        
        if(input$allCanonical) {
            val <- filter(val, all_canonical == "canonical")
        }
        
        if (input$minCovNotNa) {
            val <- filter(val, !is.na(min_cov))
        }
        
        if (input$minCovGTZero) {
            val <- filter(val, min_cov > 0)
        }
        
        if (input$onlyGenes) {
            val <- distinct(val, associated_gene, .keep_all = TRUE)
        }
        
        return(val)
    })
    
    selected_data <- reactive({
        d <- event_data(event = "plotly_click")
        if (is.null(d)) return(NULL)
        data_to_plot() %>%
            filter((!!as.symbol(input$groupBy)) == d$customdata) %>%
            filter(name %in% sort(unique(data_to_plot()$name))[d$pointNumber + 1])
    })

    output$inputTable <- renderTable({
        validate(
            need(classifications(), "Please add a classification file.")
        )
        classifications() %>% select(name, file) %>% distinct(name, .keep_all = TRUE)
    })
    
    output$count_plot <- renderPlotly({
        data_to_plot() %>% count() %>%
            ggplot(., aes_string(x="name", y="n", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            xlab("Name") + ylab("N") + ggtitle("Group Counts")
    })
    
    output$perc_plot <- renderPlotly({
        data_to_plot() %>% group_by(name) %>% count_(input$groupBy) %>% mutate(perc = n / sum(n)) %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="fill", stat="identity") +
            geom_perc_y + xlab("Name") + ggtitle("Group Percentages")
    })
    
    output$mono_plot <- renderPlotly({
        data_to_plot() %>% count(polyexonic) %>% mutate(perc = n/(sum(n))) %>%  filter(polyexonic == "Monoexonic") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y +xlab("Name") + ggtitle("Percentage Monoexonic")
    })
    
    output$arti_plot <- renderPlotly({
        data_to_plot() %>% count(SQANTI_filter) %>% mutate(perc = n / sum(n)) %>% filter(SQANTI_filter == "Artifact") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + ggtitle("Percentage Artifacts")
    })
    
    output$novel_trans_plot <- renderPlotly({
        data_to_plot() %>% count(novel_transcript) %>% mutate(perc = n / sum(n)) %>% filter(novel_transcript == "Novel") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + ggtitle("Percentage Novel Transcripts")
    })

    output$novel_genes_plot <- renderPlotly({
        data_to_plot() %>% count(novel_gene) %>% mutate(perc = n / sum(n)) %>% filter(novel_gene == "Novel") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + ggtitle("Percentage Novel Genes")
    })
    
    output$selected_transcript_count <- renderValueBox({
        text <- "0"
        if (!is.null(selected_data())) {
            text <- paste0(selected_data() %>% count() %>% first())
        }
        return(valueBox(
            text, "Selected Transcripts", icon=icon("abacus"), color="yellow", fill=TRUE
        ))
    })
    
    output$selected_transcripts <- renderText({
        if (is.null(selected_data())) return("Click a bar")
        selected_data() %>% ungroup() %>%
            select(isoform) %>%
            top_n(100) %>%
            paste0(",")
    })

    output$epivizChart <- renderUI({
        validate(
            need(input$BED_file, "Please upload a BED file")
        )
        file <- rtracklayer::BEDFile(input$BED_file[["datapath"]])

        igv <- epivizChart(
            file,
            datasource_name = "file1",
            chr="chr1", start=5870967, end=6357783)

        igv$render_component(shiny=TRUE)
    })
    
    # Download the filtered data
    output$downloadData <- downloadHandler(
        filename = "SQANTI_explorer.csv",
        content = function(file) {
            write.csv(
                classifications() %>% select(-file, -path),
                file, row.names = FALSE
            )
        }
    )
    
    output$downloadFilteredData <- downloadHandler(
        filename = "SQANTI_explorer_filtered.csv",
        content = function(file) {
            write.csv(
                data_to_plot() %>% select(-file, -path),
                file, row.names = FALSE
            )
        }
    )
    
    output$downloadSelectedData <- downloadHandler(
        filename = "SQANTI_explorer_selected.csv",
        content = function(file) {
            write.csv(
                selected_data() %>% select(-file, -path),
                file, row.names = FALSE
            )
        }
    )
})
