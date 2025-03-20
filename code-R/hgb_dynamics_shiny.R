# Install required packages if not already installed:
# install.packages(c("shiny", "deSolve", "ggplot2"))

library(shiny)
library(deSolve)
library(ggplot2)

# Define the ODE system for iron metabolism
ironODE <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Conversion: Effective hemoglobin concentration (g/dL)
    # Healthy state: y = 1.62 corresponds to Hb = 13 g/dL.
    k <- 13 / 1.62
    Hb <- y * k  # effective Hb concentration
    
    # Compute the logistic parameters for absorption and erythropoiesis
    ca <- y0 - Ba * log((amax - a0) / (a0 - amin))
    ch <- y0 - Bh * log((hmax - h0) / (h0 - hmin))
    
    # Logistic function for iron absorption fraction:
    absp <- amin + (amax - amin) / (1 + exp((Hb - ca) / Ba))
    
    # Logistic function for erythropoiesis rate:
    eryth <- hmin + (hmax - hmin) / (1 + exp((Hb - ch) / Bh))
    
    # Differential equations:
    # x: iron in "other body iron" (OBI)
    # y: iron bound in hemoglobin (Hb)
    dx <- I * absp + d * y - eryth * x - e2
    dy <- eryth * x - d * y - e1
    
    list(c(dx, dy), Hb = Hb, absp = absp, eryth = eryth)
  })
}

# Define the Shiny UI
ui <- fluidPage(
  titlePanel("Iron Metabolism Model Simulation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("I", "Daily Iron Intake (mg/day):", 
                  min = 10, max = 100, value = 55, step = 1),
      sliderInput("e1", "Menstrual Iron Loss (mg/day):", 
                  min = 0.1, max = 1.0, value = 0.60, step = 0.05),
      sliderInput("e2", "Non-menstrual Iron Loss (mg/day):", 
                  min = 0.5, max = 2.0, value = 1.06, step = 0.05),
      sliderInput("d", "RBC Death Rate (per day):", 
                  min = 0.001, max = 0.01, value = 0.0055, step = 0.0001),
      sliderInput("a0", "Normal Absorption Fraction (a0):", 
                  min = 0.01, max = 0.05, value = 0.03, step = 0.005),
      sliderInput("amin", "Minimum Absorption Fraction (amin):", 
                  min = 0.001, max = 0.02, value = 0.01, step = 0.001),
      sliderInput("amax", "Maximum Absorption Fraction (amax):", 
                  min = 0.1, max = 0.5, value = 0.2, step = 0.01),
      sliderInput("Ba", "Steepness of Absorption Curve (Ba):", 
                  min = 0.1, max = 2.0, value = 0.8, step = 0.1),
      sliderInput("h0", "Baseline Erythropoiesis Rate (h0, g/day):", 
                  min = 0.005, max = 0.05, value = 0.014, step = 0.001),
      sliderInput("hmin", "Minimum Erythropoiesis Rate (hmin):", 
                  min = 0.005, max = 0.02, value = 0.0098, step = 0.001),
      sliderInput("hmax", "Maximum Erythropoiesis Rate (hmax):", 
                  min = 0.02, max = 0.2, value = 0.07, step = 0.005),
      sliderInput("Bh", "Steepness of Erythropoiesis Curve (Bh):", 
                  min = 0.1, max = 2.0, value = 0.5, step = 0.1),
      sliderInput("time", "Simulation Time (days):", 
                  min = 30, max = 365, value = 120, step = 10)
    ),
    mainPanel(
      plotOutput("timePlot"),
      verbatimTextOutput("modelInfo")
    )
  )
)

# Define the Shiny server
server <- function(input, output) {
  parameters <- reactive({
    list(
      I = input$I,
      e1 = input$e1,
      e2 = input$e2,
      d  = input$d,
      a0 = input$a0,
      amin = input$amin,
      amax = input$amax,
      Ba = input$Ba,
      h0 = input$h0,
      hmin = input$hmin,
      hmax = input$hmax,
      Bh = input$Bh,
      y0 = 13  # set point for hemoglobin (g/dL)
    )
  })
  
  # Initial conditions:
  # For a healthy woman, 1.62 g Fe in hemoglobin and 0.7 g in OBI.
  state <- reactive({
    c(x = 0.7, y = 1.62)
  })
  
  times <- reactive({
    seq(0, input$time, by = 0.1)
  })
  
  sim <- reactive({
    ode(y = state(), times = times(), func = ironODE, parms = parameters())
  })
  
  output$timePlot <- renderPlot({
    sim_data <- as.data.frame(sim())
    # Plot effective hemoglobin (converted from y) and the two compartments
    p <- ggplot(sim_data, aes(x = time)) +
      geom_line(aes(y = Hb, color = "Effective Hb (g/dL)")) +
      geom_line(aes(y = x, color = "Other Body Iron (g)")) +
      geom_line(aes(y = y, color = "Hb Iron (g)")) +
      labs(x = "Time (days)", y = "Value", color = "Variable") +
      theme_minimal()
    print(p)
  })
  
  output$modelInfo <- renderPrint({
    tail(as.data.frame(sim()), 5)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)