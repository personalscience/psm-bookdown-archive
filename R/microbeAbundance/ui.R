#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



library(shiny)
library(tidyverse)


# Define UI for application that draws the microbe abundances
shinyUI(fluidPage(

  # Application title
  titlePanel("How common is that microbe?"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      uiOutput("taxa_names"),
   #   uiOutput("health_condition"),
      sliderInput("bins", "Bins:", value = 5, min = 3, max = 30),
   h3("URL components"),
   verbatimTextOutput("urlText")
      )

  ,
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
