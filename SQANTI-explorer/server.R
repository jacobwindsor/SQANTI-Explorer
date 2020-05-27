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
library(rtracklayer)
library(ggsci)

options(shiny.maxRequestSize = 300*1024^2)

# Plot helpers
geom_perc_y <- scale_y_continuous(labels = scales::percent_format(), name="Percentage (%)")
geom_label_stacked = geom_text(size = 3, position = position_stack(vjust = 0.5))
geom_label_dodge = geom_text(size=3, position=position_dodge(width=0.9), vjust=1.5)

bar_plot = function(data, input) {
    ggplot(data, aes_string(x="name", y="n", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity")
}


# Download Transcripts
write_isoforms <- function(data, filename, type = "gtf", cut=FALSE) {
    gtf_file <- data %>% ungroup() %>% distinct(gtf_path)
    # Ensure data contains only one file
    if (gtf_file %>% count() %>% dplyr::first() > 1) {
        print("Has more than one GTF file in input")
        return()
    }
    orig_transcripts <- rtracklayer::import.gff2(gtf_file %>% dplyr::first())
    transcripts_in_data <- NULL
    if(cut) {
        transcripts_in_data <-  data %>% ungroup() %>% dplyr::slice(0:1000) %>% pull(isoform)
    }
    else {
        transcripts_in_data <- data %>% ungroup() %>% pull(isoform)
    }
    filtered_gtf <- subset(orig_transcripts, transcript_id %in% transcripts_in_data)
    if(type=="gtf"){
        return(rtracklayer::export.gff2(filtered_gtf, filename))   
    }
    if(type=="gff"){
        return(rtracklayer::export.gff3(filtered_gtf, filename))   
    }
    if(type=="bed"){
        return(rtracklayer::export.bed(filtered_gtf, filename))  
    }
}

# Write classification file
write_classification <- function(data, filename) {
    write.csv(
        data %>% select(-file, -path, -gtf_path, -genome),
        filename, row.names = FALSE
    )
}

write_igv_session <- function(data, filename) {
    # Ensure data only for one name
    if (data %>% ungroup() %>% distinct(name) %>% count() %>% dplyr::first() > 1) {
        print("Has more than one name in input")
        return()
    }
    ungrouped <- data %>% ungroup()
    name <- ungrouped %>% select(name) %>% dplyr::first()
    genome <- ungrouped %>% select(genome) %>% dyplr::first()
    
}

# Add the temp beds file to the resources so can be accessed by URL
addResourcePath('temp_beds', "temp_beds")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    classifications <- reactiveVal()
    makeReactiveBinding("classifications")
    
    observeEvent(input$addClassification, {
        req(input$classification_file, input$gtf_file, input$name)
        
        if(input$addClassification > 0) {
            isolate({
                datapath <- input$classification_file[["datapath"]]
                new_row <- tibble(
                    name = input$name,
                    file = input$classification_file[["name"]],
                    gtf_path = input$gtf_file[["datapath"]],
                    genome = input$genome,
                    path = datapath,
                    classification = list(read_tsv(datapath))) %>%
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
            xlab("Name") + ylab("N") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$perc_plot <- renderPlotly({
        data_to_plot() %>% group_by(name) %>% count_(input$groupBy) %>% mutate(perc = n / sum(n)) %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="fill", stat="identity") +
            geom_perc_y + xlab("Name") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$mono_plot <- renderPlotly({
        data_to_plot() %>% count(polyexonic) %>% mutate(perc = n/(sum(n))) %>%  filter(polyexonic == "Monoexonic") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y +xlab("Name") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$arti_plot <- renderPlotly({
        data_to_plot() %>% count(SQANTI_filter) %>% mutate(perc = n / sum(n)) %>% filter(SQANTI_filter == "Artifact") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$novel_trans_plot <- renderPlotly({
        data_to_plot() %>% count(novel_transcript) %>% mutate(perc = n / sum(n)) %>% filter(novel_transcript == "Novel") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + theme_bw() + scale_fill_locuszoom()
    })

    output$novel_genes_plot <- renderPlotly({
        data_to_plot() %>% count(novel_gene) %>% mutate(perc = n / sum(n)) %>% filter(novel_gene == "Novel") %>%
            ggplot(., aes_string(x="name", y="perc", fill=input$groupBy, customdata=input$groupBy)) + geom_bar(position="dodge", stat="identity") +
            geom_perc_y + xlab("Name") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$gene_expression <- renderPlot({
        data_to_plot() %>% distinct(associated_gene, .keep_all=TRUE) %>%
            ggplot(., aes_string(x="name", y="log_gene_exp", fill=input$groupBy, customdata=input$groupBy)) +
            geom_violin(position=position_dodge(0.9), stat="ydensity", trim = TRUE) +
            geom_boxplot(width = 0.15, position=position_dodge(0.9), outlier.shape=NA, ) +
            xlab("Species") + ylab("Log(TPM)") + theme_bw() + scale_fill_locuszoom()
    })
    
    output$selected_transcript_count <- renderValueBox({
        text <- "0"
        if (!is.null(selected_data())) {
            text <- paste(selected_data() %>% ungroup() %>% count() %>% dplyr::first())
        }
        valueBox(
            text, "Selected Transcripts", icon=icon("calculator"), color="yellow"
        )
    })

    output$selectGenomeData <- renderUI({
        validate(
            need(classifications(), "Please add at least one dataset.")
        )
        selectInput("dataset_choice", "Select Data To Visualize: ", choices = classifications() %>% distinct(name) %>% pull(name))
    })
    
    data_to_view <- reactive(data_to_plot() %>% ungroup() %>% filter(name == input$dataset_choice))
    genome_name <- reactive(data_to_view() %>% distinct(genome) %>% dplyr::first())
    
    output$igv <- renderIgvShiny({
        igvShiny(list(
            genomeName="hg19",
            initialLocus="chr1:7,063,368-14,852,449"
        ))
    })
    
    observeEvent(input$render_igv, {
        print("running this")
        if(input$render_igv > 0) {
            isolate({
                igvShiny::loadGenome(session, list(genomeName = genome_name()))
                bed_file_name <- "temp_beds/temp.gff"
                as_url <- bed_file_name
                print(as_url)
                write_isoforms(data_to_view(), bed_file_name, "gff", cut=TRUE)
                print("written")
                igvShiny::loadGffTrackUrl(session, trackName = input$dataset_choice, url = as_url, deleteTracksOfSameName=FALSE, color = "red")
            })   
        }
    })
    
    output$downloadData <- downloadHandler(
        filename = "SQANTI_explorer.csv",
        content = function(file) {
            write_classification(classifications(), file)
        }
    )
    
    output$downloadFilteredData <- downloadHandler(
        filename = "SQANTI_explorer_filtered.csv",
        content = function(file) {
           write_classification(data_to_plot(), file)
        }
    )
    
    output$downloadSelectedData <- downloadHandler(
        filename = "SQANTI_explorer_selected.csv",
        content = function(file) {
            write_classification(selected_data(), file)
        }
    )
    
    output$downloadSelectedGTF <- downloadHandler(
        filename = "SQANTI_explorer_selected.gtf",
        content = function(file) {
            write_isoforms(selected_data(), file)
        }
    )
    
    output$downloadFilteredGTF <- downloadHandler(
        filename = "SQANTI_explorer_filtered_gtf.zip",
        content = function(file){
            #go to a temp dir to avoid permission issues
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            files <- data_to_plot() %>% ungroup() %>% distinct(name) %>% pull(name) %>%
                walk(function(row) {
                    write_isoforms(data_to_plot() %>% filter(name == row), paste0(row, ".gtf") )
                }) %>%
                map_chr(function(row) {
                    paste0(row, ".gtf") 
                })

            #create the zip file
            zip(file, files)
        }
    )
    
    output$downloadSelectedBED <- downloadHandler(
        filename = "SQANTI_explorer_selected.bed",
        content = function(file) {
            write_isoforms(selected_data(), file)
        }
    )
})
