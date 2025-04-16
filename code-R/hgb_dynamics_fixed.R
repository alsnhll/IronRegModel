# If needed, install packages:
# install.packages(c("deSolve", "ggplot2", "dplyr", "patchwork"))

{
library(deSolve)
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)
}

###############################################################################
# 0) Input Parameters 
###############################################################################

input_params <- list(
  weight = 55,                     # [kg] female weight
  e1   = 0.00060,                  # [g/day] baseline daily menstrual excretion, as long as body Fe > 0. (0.001)
  e2   = 0.00106,                  # [g/day] baseline daily other excretion, as long as body Fe > 0. (0.001)
  Tend   = 24,                     # [months] time to run simulation
  y0   = 13,                       # [g/dL] "Healthy" Hb levels
  hb_int =  c(8, 10, 12),          # [g/dL] Initial level of iron in Hb (anemia < 12)
  red_men = c(0, 25, 37, 50),      # % Percent reduction in Fe lost to menstruation
  hev_per = c(1, 2, 4),            # Number of times heavier menstruation
  Int = c(0, 5, 10, 20)            # Daily intake supplement [mg/day]
)

###############################################################################
# 1) Fixed parameters from the supplement/paper
###############################################################################

{
base_params <- list(
  d    = 0.0055,                   
  a0   = 0.03,                     
  amin = 0.01,                     
  amax = 0.20,                     
  Ba   = 0.1,                      
  Bh   = 0.062                     
)

# Derived constants
PV <- input_params$weight*0.2*0.2            # healthy plasma volume
BV <- PV/(1-0.38)                            # blood volume, 0.38 is healthy hematocrit
conv <- 285/(10*BV)                          # 285 is conversion of g Fe to g Hb, BV is blood volume
}

###############################################################################
# 2) Supporting functions
###############################################################################
{
#Absorption Fraction
absp <- function(hb) {
  ca <- input_params$y0/conv - base_params$Ba * 
    log((base_params$amax - base_params$a0) / (base_params$a0 - base_params$amin))
  a <- base_params$amin + (base_params$amax - base_params$amin) / 
    (1 + exp((hb - ca) / base_params$Ba))
  return(a)
}

#Erythropoiesis Fraction
eryth <- function(hb) {
  h0 <- (base_params$d * (input_params$y0 / conv) + input_params$e1) / 0.7
  hmin <- 0.7 * h0
  hmax <- 5 * h0
  ch <- input_params$y0 / conv - base_params$Bh * log((hmax - h0) / (h0 - hmin))
  e <- hmin + (hmax - hmin) / (1 + exp((hb - ch) / base_params$Bh))
  return(e)
}

#Initial Iron Stores in Other Body
OBI0 <- function(hb_int, e1) {
  o <- ((hb_int / conv) * base_params$d + e1) / eryth(hb_int / conv)
  return(o)
}

#Steady-state intake to maintain anemia
I <- function(dep_hb, e1) {
  i <- (e1 + input_params$e2) / absp(dep_hb / conv)
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
    
    dx <- I * absp(y) + d * y - eryth(y) * x - input_params$e2
    dy <- eryth(y) * x - d * y - e1
    
    return(list(c(dx, dy)))
  })
}

###############################################################################
# 4) Run Simulation
###############################################################################

run_sim <- function(hb_int, Int, red_men, hev_per) {
  
  #Simulation parameters
  OBI <- OBI0(hb_int, input_params$e1 * hev_per)
  Int_ss <- I(hb_int, input_params$e1 * hev_per)
  times <- seq(0, input_params$Tend * 30, by = 1)
  parameters <- c(base_params, list(
    I = Int_ss + Int / 1000,
    e1 = input_params$e1 * (1 - red_men / 100) * hev_per
  ))
  
  #Initial State of ODE
  state_init <- c(x = OBI, y = hb_int / conv)
  
  #Output of ODE
  out <- ode(
    y = state_init,
    times = times,
    func = ironODE,
    parms = parameters
  )
  
  out_df <- as.data.frame(out)
  colnames(out_df) <- c("time", "x", "y")
  
  out_df <- out_df %>%
    mutate(
      hb = y * conv,
      hb_int = hb_int,
      intake = Int,
      red_men = red_men,
      hev_per = hev_per
    )
  
  return(out_df)
}

###############################################################################
# 4) Display Results
###############################################################################

{
#Print Initial States
init_states <- data.frame(
  scenario = paste0(input_params$hb_int, " g/dL"),
  x_init = OBI0(input_params$hb_int, input_params$e1),
  y_init = (input_params$hb_int)/conv, 
  I_SS = I(input_params$hb_int, input_params$e1)
)

print(init_states)
}

# All scenarios tested
{
scenarios <- expand.grid(
  hb_int = input_params$hb_int,
  hev_per = input_params$hev_per,
  red_men = input_params$red_men,
  int = input_params$Int,
  treatment = c("No Treatment", "TXA Only", "Supplement Only", "TXA + Supplement"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    int = case_when(
      treatment %in% c("No Treatment", "TXA Only") ~ 0,
      TRUE ~ int
    ),
    red_men = case_when(
      treatment %in% c("No Treatment", "Supplement Only") ~ 0,
      TRUE ~ red_men
    )
  ) %>%
  filter(
    !(treatment == "TXA Only" & red_men == 0),
    !(treatment == "Supplement Only" & int == 0), 
    !(treatment == "TXA + Supplement" & int == 0),
    !(treatment == "TXA + Supplement" & red_men == 0),
  ) %>%
  select(hb_int, hev_per, int, red_men, treatment) %>%
  distinct()

scenarios <- scenarios %>%
  mutate(
    treatment = case_when(
      treatment == "Supplement Only" ~ paste0("Supplement Only: ", int, " mg"),
      treatment == "TXA + Supplement" ~ paste0("TXA + Supplement: ", int, " mg"),
      TRUE ~ treatment
    )
  )

print(scenarios)
}

###############################################################################
# 5) Plotting Results
###############################################################################

# Run all simulations
{scenarios$row <- seq_len(nrow(scenarios))

all_results <- pmap_dfr(
  list(
    hb_int = scenarios$hb_int,
    Int = scenarios$int,
    red_men = scenarios$red_men,
    hev_per = scenarios$hev_per,
    row = scenarios$row
  ),
  function(hb_int, Int, red_men, hev_per, row) {
    sim <- run_sim(hb_int, Int, red_men, hev_per)
    sim$treatment <- scenarios$treatment[[row]]
    sim$hb_int <- hb_int
    sim$hev_per <- hev_per
    return(sim)
  }
)
}

#Filter to normal vs. heavy periods
unique_hev <- unique(all_results$hev_per)

#Conditional plotting for normal periods vs. heavy periods
if (length(unique_hev) == 1 && unique_hev == 1) {

normal_periods_data <- all_results %>%
  filter(hev_per == 1) %>%  # or hev_per != 1 for heavy periods
  mutate(facet_label = ifelse(red_men == 0, NA, paste0("TXA efficacy: ", red_men, "%")))

txa_levels <- sort(unique(all_results$red_men))
txa_levels <- txa_levels[txa_levels != 0]

baseline_rows <- normal_periods_data %>% filter(is.na(facet_label))
txa_rows <- normal_periods_data %>% filter(!is.na(facet_label))

# Duplicate baseline into each facet
baseline_duplicated <- map_dfr(txa_levels, function(eff) {
  baseline_rows %>%
    mutate(facet_label = paste0("TXA efficacy: ", eff, "%"))
})

normal_periods <- bind_rows(baseline_duplicated, txa_rows)
print(normal_periods)

# ------------------------------------------------------------------------------
# Plot: Hb over time
# ------------------------------------------------------------------------------

p_hb <- ggplot(normal_periods, aes(x = time / 30, y = hb, color = factor(hb_int), linetype = treatment)) +
  geom_line(linewidth = 0.5) +
  facet_grid(~ facet_label, scales = "free_y") +
  labs(
    x = "Time (months)",
    y = "Hemoglobin (g/dL)",
    color = "Initial Hb Level",
    linetype = "Treatment"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  theme(strip.text = element_text(size = 10, face = "bold"))

# ------------------------------------------------------------------------------
# Plot: Body iron over time
# ------------------------------------------------------------------------------

p_iron <- ggplot(normal_periods, aes(x = time / 30, y = x, color = factor(hb_int), linetype = treatment)) +
  geom_line(linewidth = 0.5) +
  facet_grid(~ facet_label, scales = "free_y") +
  labs(
    x = "Time (months)",
    y = "Body Iron (g)",
    color = "Initial Hb Level",
    linetype = "Treatment"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  theme(strip.text = element_text(size = 10, face = "bold"))

# ------------------------------------------------------------------------------
# Combine plots vertically
# ------------------------------------------------------------------------------

p_hb / p_iron + plot_layout(guides = "collect") #& theme(legend.position = "bottom")

} else {
  
  heavy_periods_data <- all_results %>%
    filter(hb_int == input_params$hb_int[1]) %>%  # or hev_per != 1 for heavy periods
    mutate(facet_label = ifelse(red_men == 0, NA, paste0("TXA efficacy: ", red_men, "%")))
  
  txa_levels <- sort(unique(all_results$red_men))
  txa_levels <- txa_levels[txa_levels != 0]
  
  baseline_rows <- heavy_periods_data %>% filter(is.na(facet_label))
  txa_rows <- heavy_periods_data %>% filter(!is.na(facet_label))
  
  # Duplicate baseline into each facet
  baseline_duplicated <- map_dfr(txa_levels, function(eff) {
    baseline_rows %>%
      mutate(facet_label = paste0("TXA efficacy: ", eff, "%"))
  })
  
  heavy_periods <- bind_rows(baseline_duplicated, txa_rows)
  print(heavy_periods)
  
  # ------------------------------------------------------------------------------
  # Plot: Hb over time
  # ------------------------------------------------------------------------------
  
  p_hb <- ggplot(heavy_periods, aes(x = time / 30, y = hb, color = factor(hev_per), linetype = treatment)) +
    geom_line(linewidth = 0.5) +
    facet_grid(~ facet_label, scales = "free_y") +
    labs(
      x = "Time (months)",
      y = "Hemoglobin (g/dL)",
      color = "#x Heavier Periods",
      linetype = "Treatment"
    ) +
    theme_minimal() +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
    theme(strip.text = element_text(size = 10, face = "bold"))
  
  # ------------------------------------------------------------------------------
  # Plot: Body iron over time
  # ------------------------------------------------------------------------------
  
  p_iron <- ggplot(heavy_periods, aes(x = time / 30, y = x, color = factor(hev_per), linetype = treatment)) +
    geom_line(linewidth = 0.5) +
    facet_grid(~ facet_label, scales = "free_y") +
    labs(
      x = "Time (months)",
      y = "Body Iron (g)",
      color = "#x Heavier Periods",
      linetype = "Treatment"
    ) +
    theme_minimal() +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
    theme(strip.text = element_text(size = 10, face = "bold"))
  
  # ------------------------------------------------------------------------------
  # Combine plots vertically
  # ------------------------------------------------------------------------------
  
  p_hb / p_iron + plot_layout(guides = "collect") #& theme(legend.position = "bottom")
}
