# If needed, install packages:
# install.packages(c("deSolve", "ggplot2", "dplyr"))

library(deSolve)
library(ggplot2)
library(dplyr)

###############################################################################
# 0) Input Parameters 
###############################################################################
input_params <- list(
  weight = 55,                     # [kg] female weight
  e1   = 0.00060,                  # [g/day] baseline daily menstrual excretion, as long as body Fe > 0. (0.001)
  e2   = 0.00106,                  # [g/day] baseline daily other excretion, as long as body Fe > 0. (0.001)
  Tend   = 30,                     # [months] time to run simulation
  y0   = 13,                       # [g/dL] "Healthy" Hb levels
  dep_hb =  c(0, 20, 50),          # Percent depletion for iron in Hb
  red_men = c(0, 20, 50)           # % Percent reduction in Fe lost to menstruation
)

###############################################################################
# 1) Fixed parameters from the supplement/paper
###############################################################################
base_params <- list(
  d    = 0.0055,                     # RBC death rate (half-life ~127 days)
  a0   = 0.03,                       # Absorption fraction at healthy Hb
  amin = 0.01,                       # Minimum absorption fraction
  amax = 0.20,                       # Maximum absorption fraction
  Ba   = 0.1,                        # Steepness of absorption curve
  Bh   = 0.062                       # Steepness of erythropoiesis curve
)

###############################################################################
# 2) Define the ODE system
###############################################################################

PV <- input_params$weight*0.2*0.2            # healthy plasma volume
BV <- PV/(1-0.38)                            # blood volume, 0.38 is healthy hematocrit
k <- 285/(10*BV)                             # 285 is conversion of g Fe to g Hb, BV is blood volume

# Function to solve the system of ODEs
ironODE <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dx <- I * absp(y) + d * y - eryth(y) * x - input_params$e2
    dy <- eryth(y) * x - d * y - e1
    #Hb <- y*k
    
    list(c(dx, dy))#,  Hb = Hb)  
  })
}

# Function calculating the absorption rate
absp <- function(hb) {
  
  # Iron absorption fraction

  ca <- input_params$y0/k - base_params$Ba * log((base_params$amax - base_params$a0) / (base_params$a0 - base_params$amin))
  a <- base_params$amin + (base_params$amax - base_params$amin) / (1 + exp((hb - ca) / base_params$Ba))
  return(a)
}

# Function calculating the erythropoiesis rate
eryth <- function(hb) {
  
  # Iron absorption fraction
  
  h0 <- (base_params$d*(input_params$y0/k)+input_params$e1)/0.7 # 0.013 ## Fix this
  hmin <- 0.7*h0
  hmax <- 5*h0
  ch <- input_params$y0/k - base_params$Bh * log((hmax - h0) / (h0 - hmin))
  e <- hmin + (hmax - hmin) / (1 + exp((hb - ch) / base_params$Bh))
  return(e)
}

###############################################################################
# 3) Define initial conditions: Various levels of anemia
###############################################################################
# - Healthy: x=0.7 g OBI, y=1.62 g Fe in Hb (~13 g/dL)
# - Anemic:  x=0.2 g OBI, y=1.0  g Fe in Hb (~8 g/dL)

# Function calculating the initial OBI reserves
OBI0 <- function(dep_hb) {
  o <- (((1-(dep_hb/100))*input_params$y0/k)*base_params$d+input_params$e1)/eryth((1-(dep_hb/100))*input_params$y0/k)
  return(o)
}

# Function calculating the intake to maintain steady-state anemia conditions under "normal" menstruation
I <- function(dep_hb) {
  i <- (input_params$e1+input_params$e2)/absp((1-(dep_hb/100))*input_params$y0/k)

  return(i)
}

init_states <- data.frame(
  scenario = c("Healthy", "20% Depletion", "50 % Depletion"),
  x_init = OBI0(input_params$dep_hb),
  y_init = ((1-(input_params$dep_hb/100))*input_params$y0/k), 
  I_vals = I(input_params$dep_hb)
)

print(init_states)

###############################################################################
# 4) Run simulations for each scenario/e1
###############################################################################
e1_vals <- ((1-(input_params$red_men/100)))*input_params$e1
times <- seq(0, input_params$Tend*30, by = 1)

sim_data_list <- list()

for (i in 1:length(init_states$scenario)) {
  # get the correct initial states for the scenario
  scen <- init_states[i,1]
  intake <- init_states[i,4]
  for (e1val in e1_vals) {
     # combine base parameters with the chosen e1 and I
      params_i <- c(base_params, list(e1 = e1val, I = intake))
      
      # initial conditions
      state_i <- c(x = init_states[i,2], y = init_states[i,3])
      
      # solve ODE
      out <- ode(y = state_i, times = times, func = ironODE, parms = params_i)
      out_df <- as.data.frame(out)
      
      # annotate
      out_df$scenario <- scen
      out_df$e1       <- e1val*1000
      out_df$I        <- intake*1000
      
      sim_data_list[[length(sim_data_list) + 1]] <- out_df
  }
}

# Combine all runs into one data frame
sim_data <- do.call(rbind, sim_data_list)

###############################################################################
# 5) Plot: 2 (rows) x 5 (columns) facet grid: scenario vs e1
#    Within each panel, color by I
###############################################################################
p1 <- ggplot(sim_data, aes(x = time/30, y = y*k, color = factor(scenario), linetype = factor(e1))) +
  geom_line() +
  #facet_grid(scenario ~ e1, labeller = label_both) +
  labs(
    #title = "Hemoglobin Dynamics Over 3 months",
    subtitle = "Varying scenario (Healthy vs. Anemic) and menstrual loss (e1)",
    x = "Time (months)",
    y = "Hemoglobin (g/dL)",
    color = "% Anemia",
    linetype = 'Fe lost to menstruation \n(mg/day)'
  ) +
  theme_minimal() + 
  lims(y = c(6, 14))

print(p1)

p2 <- ggplot(sim_data, aes(x = time/30, y = x, color = factor(scenario),linetype = factor(e1))) +
  geom_line() +
  #facet_grid(scenario ~ e1, labeller = label_both) +
  labs(
    #title = "OBI Levels Over 30 months",
    subtitle = "Varying scenario (Healthy vs. Anemic) and menstrual loss (e1)",
    x = "Time (months)",
    y = "Other Body Fe (g)",
    color = "% Anemia",
    linetype = 'Fe lost to menstruation \n(mg/day)'
  ) +
  theme_minimal() + 
  lims(y = c(0, 1))

print(p2)
