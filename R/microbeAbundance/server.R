#
# This is a Shiny app to display microbe abundances
#

library(shiny)
library(tidyverse)

microbiome_samples_genus_df <- # read.csv("microbiome_samples_genus.csv")
   read.csv(file.path(here::here(), "R","microbeAbundance","microbiome_samples_genus.csv"))


#people.healthy.gut.genus <- phyloseq::subset_samples(psmr::people.norm, Site == "gut" & Condition == "Healthy")

# Define server logic required to draw food and microbes
server <- function(input, output) {

  output$taxa_names <- renderUI(

    selectInput(inputId = "dataset",
                label = "Choose a Microbe:",
                choices = microbiome_samples_genus_df %>% dplyr::filter(abundance > 0 & condition == "Healthy") %>% pull(taxa)
                %>% unique() %>% sort() %>% as.vector()
    ))

  output$health_condition <- renderUI(
    selectInput(inputId = "health_status",
                label = "Health Status:",
                choices = levels(microbiome_samples_genus_df$condition))
  )



  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    taxa_name <- input$dataset
    condition <- input$health_status

    taxa_name <- ifelse(is.null(taxa_name),"Bifidobacterium",taxa_name)

    microbiome_samples_genus_df %>% dplyr::filter(taxa == taxa_name &
                                                    condition == ifelse(is.null(input$health_status),
                                                                        "Healthy",
                                                                        input$health_status) ) %>%
      ggplot(aes( x = abundance)) +
      geom_histogram(bins = input$bins) +
      labs(title = paste("Abundance of", taxa_name, "in", condition),
           x = "Abundance (%)",
           y = "Frequency") +
      theme(axis.text.y   = element_blank())


  })
}

