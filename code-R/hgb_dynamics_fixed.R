# If needed, install packages:
# install.packages(c("deSolve", "ggplot2", "dplyr"))

library(deSolve)
library(ggplot2)
library(dplyr)

###############################################################################
# 1) Define the ODE system
###############################################################################
ironODE <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Convert y (g Fe in Hb) to approximate Hb (g/dL).
    # 1.62 g Fe in Hb corresponds to 13 g/dL:
    weight <- 55 #Female weight in [kg]
    PV <- weight*0.2*0.2 # healthy plasma volume
    BV <- PV/(1-0.38) #blood volume, 0.38 is healthy hematocrit
    #####Adjusted to match conversion in matlab code
    #k <- 13 / 1.62 #previous expression for k
    k <- 285/(10*BV) #285 is conversion of g Fe to g Hb, BV is blood volume
    Hb <- y * k
    
    # Logistic parameters for absorption
    ca <- y0 - Ba * log((amax - a0) / (a0 - amin))
    
    # Logistic parameters for erythropoiesis
    ch <- y0 - Bh * log((hmax - h0) / (h0 - hmin))
    
    # Iron absorption fraction
    #absp <- amin + (amax - amin) / (1 + exp((Hb - ca) / Ba))
    absp <- amin + (amax - amin) / (1 + exp((Hb - ca) / Ba))
    
    # Erythropoiesis rate
    #eryth <- hmin + (hmax - hmin) / (1 + exp((Hb - ch) / Bh))
    eryth <- hmin + (hmax - hmin) / (1 + exp((Hb - ch) / Bh))
    
    # ODEs:
    dx <- I * absp + d * y - eryth * x - e2
    dy <- eryth * x - d * y - e1
    
    list(c(dx, dy),  Hb = Hb)  # Pass Hb along as output , Hb=Hb
  })
}

###############################################################################
# 2) Fixed parameters from the supplement/paper
###############################################################################
base_params <- list(
  d    = 0.0055,   # RBC death rate (half-life ~127 days)
  e2   = 0.00106,     # Non-menstrual iron loss (g/day)
  #e1   = 0.00060,     # Menstrual iron loss (g/day)
  a0   = 0.03,     # Absorption fraction at healthy Hb
  amin = 0.01,     # Minimum absorption fraction
  amax = 0.20,     # Maximum absorption fraction
  Ba   = 0.80,     # Steepness of absorption curve
  h0   = 0.014,    # Baseline erythropoiesis rate (g/day)
  hmin = 0.7*0.014,
  hmax = 5*0.014,
  Bh   = 0.50,     # Steepness of erythropoiesis curve
  y0   = 13        # Set-point for Hb (g/dL)
)

###############################################################################
# 3) Define two initial conditions: Healthy vs. Anemic
###############################################################################
# - Healthy: x=0.7 g OBI, y=1.62 g Fe in Hb (~13 g/dL)
# - Anemic:  x=0.2 g OBI, y=1.0  g Fe in Hb (~8 g/dL)
init_states <- data.frame(
  scenario = c("20% Depletion", "50 % Depletion"),
  x_init   = c(0.7/k, 0.7/k),
  y_init   = c(0.8*13/k, 0.5*13/k)
)

###############################################################################
# 4) Define parameter ranges for e1 (menstrual loss) and I (daily intake)
###############################################################################
# We'll do 5 e1 values and 5 I values:
e1_vals <- c(0.0003, 0.00048, 0.0006)   # g/day #0.05, 0.1, 0.2, 0.4, 
I_vals  <- c(0.01,0.020,0.050)  # g/day 5, 10, 15, 20, 25, 30

# We will create a 2 (scenarios) x 5 (e1) facet grid, 
# and within each panel we plot 5 different I values as lines.

###############################################################################
# 5) Simulation time: 90 days (about 3 cycles)
###############################################################################
times <- seq(0, 900, by = 1)

###############################################################################
# 6) Run simulations for each scenario/e1/I combination
###############################################################################
sim_data_list <- list()

for (scen in init_states$scenario) {
  # get the correct initial states for the scenario
  init_row <- init_states[init_states$scenario == scen, ]
  for (e1val in e1_vals) {
    for (Ival in I_vals) {
      # combine base parameters with the chosen e1 and I
      params_i <- c(base_params, list(e1 = e1val, I = Ival))
      
      # initial conditions
      state_i <- c(x = init_row$x_init, y = init_row$y_init)
      
      # solve ODE
      out <- ode(y = state_i, times = times, func = ironODE, parms = params_i)
      out_df <- as.data.frame(out)
      
      # annotate
      out_df$scenario <- scen
      out_df$e1       <- e1val*1000
      out_df$I        <- Ival*1000
      
      sim_data_list[[length(sim_data_list) + 1]] <- out_df
    }
  }
}

# Combine all runs into one data frame
sim_data <- do.call(rbind, sim_data_list)

###############################################################################
# 7) Plot: 2 (rows) x 5 (columns) facet grid: scenario vs e1
#    Within each panel, color by I
###############################################################################
p <- ggplot(sim_data, aes(x = time/30, y = Hb, color = factor(I))) +
  geom_line() +
  facet_grid(scenario ~ e1, labeller = label_both) +
  labs(
    title = "Hemoglobin Dynamics Over 90 Days (3 Menstrual Cycles)",
    subtitle = "Varying scenario (Healthy vs. Anemic), menstrual loss (e1), and daily intake (I)",
    x = "Time (months)",
    y = "Hemoglobin (g/dL)",
    color = "Daily Intake\n(mg/day)"
  ) +
  theme_minimal() + 
  lims(y = c(6, 14))

print(p)