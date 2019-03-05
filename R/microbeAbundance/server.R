#
# This is a Shiny app to display microbe abundances
#

library(shiny)
library(tidyverse)

healthy_genus_df <- read.csv(file.path(here::here(), "R","microbeAbundance","healthy_genus.csv"))


#people.healthy.gut.genus <- phyloseq::subset_samples(psmr::people.norm, Site == "gut" & Condition == "Healthy")

# Define server logic required to draw food and microbes
server <- function(input, output) {

  output$taxa_names <- renderUI(

    selectInput(inputId = "dataset",
                label = "Choose a Microbe:",
                choices = healthy_genus_df %>% filter(abundance > 0) %>% pull(taxa) %>% unique() %>% as.vector()
    )

  )

  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    taxa_name <- input$dataset

  healthy_genus_df %>% filter(taxa == taxa_name) %>%
      ggplot(aes( x = abundance)) +
      geom_histogram(bins = input$bins) +
      labs(title = paste("Abundance of", taxa_name, "in (Healthy People)"), x = "Abundance (%)")


  })
}

