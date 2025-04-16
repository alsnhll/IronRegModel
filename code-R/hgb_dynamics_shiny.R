# Install required packages if not already installed:
# install.packages(c("shiny", "deSolve", "ggplot2", "dplyr", "patchwork"))

{
  library(shiny)
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(patchwork)
}

###############################################################################
# 1) Fixed parameters from the supplement/paper
###############################################################################

{
  params <- list(
    d    = 0.0055,      
    a0   = 0.03,
    amin = 0.01,
    amax = 0.20,
    Ba   = 0.1,
    Bh   = 0.062, 
    weight = 55,               # [kg] female weight
    e1 = 0.00060,              # [g/day] baseline daily menstrual excretion, as long as body Fe > 0. (0.001)
    e2   = 0.00106,            # [g/day] baseline daily other excretion, as long as body Fe > 0. (0.001)
    y0 = 13                    # [g/dL] "Healthy" Hb levels
  )
  
  # Derived constants
  PV <- params$weight*0.2*0.2  # healthy plasma volume
  BV <- PV/(1-0.38)            # blood volume, 0.38 is healthy hematocrit
  conv <- 285/(10*BV)          # 285 is conversion of g Fe to g Hb, BV is blood volume
}

###############################################################################
# 2) Supporting functions
###############################################################################
{
  #Absorption Fraction
  absp <- function(hb) {
    ca <- params$y0/conv - params$Ba * 
      log((params$amax - params$a0) / (params$a0 - params$amin))
    a <- params$amin + (params$amax - params$amin) / (1 + exp((hb - ca) / base_params$Ba))
    return(a)
  }
  
  #Erythropoiesis Fraction
  eryth <- function(hb) {
    h0 <- (params$d * (params$y0 / conv) + params$e1) / 0.7
    hmin <- 0.7 * h0
    hmax <- 5 * h0
    ch <- params$y0 / conv - params$Bh * log((hmax - h0) / (h0 - hmin))
    e <- hmin + (hmax - hmin) / (1 + exp((hb - ch) / params$Bh))
    return(e)
  }
  
  #Initial Iron Stores in Other Body
  OBI0 <- function(hb_int, e1) {
    o <- ((hb_int / conv) * params$d + e1) / eryth(hb_int / conv)
    return(o)
  }
  
  #Steady-state intake to maintain anemia
  I <- function(dep_hb, e1) {
    i <- (e1 + params$e2) / absp(dep_hb / conv)
    return(i)
  }
}

###############################################################################
# 3) ODE system
###############################################################################

ironODE <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    x <- state["x"]  # Other Body Iron
    y <- state["y"]  # Iron in Hemoglobin
    
    dx <- I * absp(y) + d * y - eryth(y) * x - params$e2
    dy <- eryth(y) * x - d * y - e1
    
    return(list(c(dx, dy)))
  })
}

###############################################################################
# 4) Run the Simulation
###############################################################################

run_sim <- function(hb_int, Int, red_men, hev_per, Tend) {
  OBI <- OBI0(hb_int, params$e1 * hev_per)
  Int_ss <- I(hb_int, params$e1 * hev_per)
  times <- seq(0, Tend * 30, by = 1)
  parameters <- c(params, list(
    I = Int_ss + Int / 1000,
    e1 = params$e1 * (1 - red_men / 100) * hev_per
  ))
  state_init <- c(x = OBI, y = hb_int / conv)
  
  out <- ode(
    y = state_init,
    times = times,
    func = ironODE,
    parms = parameters
  )
  
  df <- as.data.frame(out)
  df$hb <- df$y * conv
  df$hb_int <- hb_int
  df$int <- Int
  df$red_men <- red_men
  df$hev_per <- hev_per
  return(df)
}

###############################################################################
# 5) Define the Shiny UI
###############################################################################

ui <- fluidPage(
  titlePanel("Iron Metabolism Simulation"),
  sidebarLayout(
    sidebarPanel(width = 3, 
       checkboxGroupInput("hb_int", "Initial Hemoglobin (g/dL):", 
                   choices = c(8, 9, 10, 11, 12), selected = 10),
       checkboxGroupInput("red_men", "Reduction of Menstruation (TXA Efficacy) (%):", 
                   choices = c(0, 25, 37, 50), selected = 25), 
       sliderInput("hev_per", "Period Heaviness (# times from normal):", 
                   min = 0.5, max = 5, value = 1, step = 0.5), 
       sliderInput("Int", "Daily Iron Supplement (mg/day):", 
                   min = 0, max = 60, value = 10, step = 1),
       sliderInput("time", "Time (Months):", 
                   min = 1, max = 30, value = 24, step = 1),
       actionButton("simulate", "Run Simulation")
    ),
    mainPanel(
      width = 9,
      plotOutput("simulationPlot", height = "600px")
    )
  ),
  
  fluidRow(
    column(12,
           DT::dataTableOutput("summary_table")
    )
  )
)

###############################################################################
# 6) Define the Shiny Server
###############################################################################

server <- function(input, output) {
  sim_data <- eventReactive(input$simulate, {
    req(length(input$hb_int), length(input$red_men), length(input$hev_per), length(input$Int))
    grid <- expand.grid(
      hb_int = as.numeric(input$hb_int),
      red_men = as.numeric(input$red_men),
      hev_per = as.numeric(input$hev_per),
      Int = as.numeric(input$Int),
      stringsAsFactors = FALSE
    )
    grid$row <- seq_len(nrow(grid))
    
    pmap_dfr(grid, function(hb_int, red_men, hev_per, Int, row) {
      df <- run_sim(hb_int, Int, red_men, hev_per, Tend = input$time)
      return(df)
    })
  })
 
  #Print a table which summarizes the scenarios and initial conditions 
  summary_table <- eventReactive(input$simulate, {
    req(length(input$hb_int), length(input$red_men), length(input$hev_per), length(input$Int))
    selected_inputs <- expand.grid(
      hb_int = as.numeric(input$hb_int),
      red_men = as.numeric(input$red_men),
      hev_per = as.numeric(input$hev_per),
      Int = as.numeric(input$Int),
      stringsAsFactors = FALSE
    )
    
    selected_inputs <- selected_inputs %>%
      mutate(
        treatment = case_when(
          Int == 0 & red_men == 0 ~ "No Treatment",
          Int == 0 & red_men > 0  ~ "TXA Only",
          Int > 0 & red_men == 0  ~ "Supplement Only",
          Int > 0 & red_men > 0   ~ "TXA + Supplement"
        ),
        treatment_label = case_when(
          treatment == "Supplement Only" ~ paste0("Supplement Only: ", Int, " mg"),
          treatment == "TXA + Supplement" ~ paste0("TXA + Supplement: ", Int, " mg"),
          TRUE ~ treatment
        )
      )
    
    selected_inputs %>%
      mutate(
        OBI = OBI0(hb_int, input_params$e1 * hev_per),
        Int_ss = I(hb_int, input_params$e1 * hev_per)
      ) %>%
      transmute(
        `Initial Hb (g/dL)` = hb_int,
        `Period Heaviness (x)` = hev_per,
        `TXA Efficacy (%)` = red_men,
        `Iron Supplement (mg/day)` = Int,
        `Treatment` = treatment_label,
        `Initial Body Iron (g)` = round(OBI, 3),
        `Steady-State Intake (mg/day)` = round(Int_ss*1000, 3)
      )
  })
  
  ###############################################################################
  # 7) Outputs
  ############################################################################### 
  
  output$simulationPlot <- renderPlot({
    req(sim_data())
    df <- sim_data()
    
    p1 <- ggplot(df, aes(x = time / 30, y = hb, color = factor(hb_int), linetype = factor(red_men))) +
      geom_line(linewidth = 0.5) +
      labs(
        x = "Time (months)", 
        y = "Hb (g/dL)", 
        color = "Initial Hb Level [g/dL]", 
        linetype = "% Reduction of Menstruation"
      ) +
      #scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
      scale_y_continuous(limits = c(8, 13), breaks = seq(8, 13, 0.5)) +
      theme_minimal(base_size = 16) 
    
    p2 <- ggplot(df, aes(x = time / 30, y = x, color = factor(hb_int), linetype = factor(red_men))) +
      geom_line(linewidth = 0.5) +
      labs(
        x = "Time (months)", 
        y = "Body Iron (g)",
        color = "Initial Hb Level [g/dL]", 
        linetype = "% Reduction of Menstruation"
      ) +
      #scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
      #scale_y_continuous(limits = c(0.05, 0.5), breaks = seq(0.05, 0.5, 0.05)) +
      theme_minimal(base_size = 16) 
    
    p1 / p2 + plot_layout(guides = "collect")  & 
      theme(legend.position = "bottom", legend.direction = "horizontal")
  })
  
  output$summary_table <- DT::renderDataTable({
    summary_table()
  }, options = list(
    pageLength = 10,
    scrollX = TRUE
  )
  )
}

###############################################################################
# 8) Run the Shiny App
############################################################################### 
shinyApp(ui = ui, server = server)