# 2 Feb 2024
# KACollier
# Create stacked barplot of our AllInds dataset for K2-K9

# Load libraries
library(tidyverse)
library(patchwork)
library(rlist)

# Helper functions
parse_q_and_fam <- function(q_loc, fam_loc) {
  ### READ IN Q FILE ###
  message("Reading in Q-file:")
  message("")
  
  q_data <- readr::read_table(file = q_loc, col_names = FALSE)
  num_K <- ncol(q_data)
  
  ### READ IN FAM FILE ###
  message("Reading in .fam file:")
  message("")
  
  fam_cols <- c("population","ind")
  fam_data <- readr::read_table(file = fam_loc, 
                                col_names = fam_cols)
  
  ### BIND Q AND FAM TIBBLES BY COLUMN ###
  Kpop <- bind_cols(fam_data, q_data) %>%
    arrange(population)
  return(Kpop)
}
parse_pop_order <- function(order_loc) {
  message("Reading in population order:")
  message("")
  pop_levels_data <- readr::read_table(order_loc, col_names = "population")
  pop_levels <- pop_levels_data$population
  return(pop_levels)
}
label_translate <- function(level_string) {
  # uses stringr to do basic transformations on population levels to make them nice labels
  label_string <- str_replace(level_string, "N_", "North ") %>%
    str_replace("S_", "South ") %>%
    str_replace("E_", "East ") %>%
    str_replace("W_", "West ") %>%
    str_replace("C_", "Central ")
  
  return(label_string)
}

get_long_K_tibble <- function(tibble = Kpop, pop_order, pop_labels) {
  ### CREATE K_LIST - ONE ELEMENT FOR EVERY VALUE OF K ###
  k_val <- (ncol(tibble) - 2)
  new_list <- vector(mode = 'list', length = k_val)
  
  # loop to fill all slots in 'new_list', using get_K_tibble() helper function
  for (i in 1:k_val) {
    colname <- str_c("X", i)
    new_list[[i]] <- get_K_tibble(tibble, eval(colname), i)
  }
  
  # Reshape 'new_list' into long data frame
  long_K_tibble <- bind_rows(new_list)
  
  # Mutate/reformat columns as factors
  long_K_tibble <- long_K_tibble %>%
    mutate(K = factor(as_factor(long_K_tibble$K), levels = seq(from = k_val, to = 1, by = -1)), # recode and relevel K to act like a factor
           population = factor(as_factor(long_K_tibble$population), # Recode <chr> population as factor
                               levels = pop_order,
                               labels = pop_labels)) # Relevel <fct> population to correctly order populations in plot
  return(long_K_tibble)
}
get_K_tibble <- function(tibble = Kpop, colname, k_val) {
  # This function pulls out the values for a single K, along with population/individual columns
  # It then writes them to a new tibble and returns that tibble
  single_K <- tibble %>%
    select(population, ind, colname) %>%
    rename(percent = colname) %>%
    mutate(K = k_val) %>%
    select(population, ind, K, percent)
  return(single_K)
}

get_ind_order <- function(tibble = Kpop) {
  #Helper function for 'mass_pop_plot_tile'
  # This function takes your Kpop tibble and returns a list of individuals sorted by population
  # Output serves as the ind_levels argument in that function
  ind_order <- tibble %>%
    filter(K == 1) %>% # removes all higher Ks, guaranteeing a single value for each ind
    arrange(population) %>% # sorts inds by pop. Could probably sort further, if desired
    pull(ind) %>% # pulls 'ind' column from tibble and converts to a list
    as_factor() # converts list to factor
  
  return(ind_order) # new factor object
}

# plotting functions
get_top_plot <- function(tibble) {
  
  # Set theme
  top_theme <- theme(axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.ticks.length = unit(0, "pt"), #length of tick marks
                     axis.title = element_blank(),
                     axis.line.y = element_blank(),
                     axis.line.x.bottom = element_blank(),
                     plot.margin = margin(r = 30, b = 0, l = 30, ),
                     legend.position = "none")
  
  # Get order of individuals
  ind_order <- get_ind_order(tibble)
  
  # Pipe Kpop to ggplot
  top_plot <- tibble %>%
    mutate(ind = factor(ind, levels = ind_order)) %>% #format data frame to order individuals
    # Begin plot command
    ggplot(aes(x = reorder(ind, population), y = percent, fill = K)) +
    geom_col(position = "fill") +
    theme_classic() +
    facet_grid(~population, space = "free", scales = "free_x") +
    top_theme
  
  # return output plot
  return(top_plot)
}
get_mid_plot <- function(tibble) {
  
  # Set theme
  mid_theme <- theme(axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.ticks.length = unit(0, "pt"), #length of tick marks
                     axis.title = element_blank(),
                     axis.line.y = element_blank(),
                     axis.line.x.bottom = element_blank(),
                     strip.text.x.top = element_blank(),
                     strip.background.x = element_blank(),
                     plot.margin = margin(r = 30, b = 0, l = 30, ),
                     legend.position = "none")
  
  # Get order of individuals
  ind_order <- get_ind_order(tibble)
  
  # Pipe Kpop to ggplot
  mid_plot <- tibble %>%
    mutate(ind = factor(ind, levels = ind_order)) %>% #format data frame to order individuals
    # Begin plot command
    ggplot(aes(x = reorder(ind, population), y = percent, fill = K)) +
    geom_col(position = "fill") +
    theme_classic() +
    facet_grid(~population, space = "free", scales = "free_x") +
    mid_theme
  
  # return output plot
  return(mid_plot)
}
get_bottom_plot <- function(tibble) {
  
  # Set theme
  bottom_theme <- theme(axis.text.y = element_blank(),
                        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank(),
                        axis.line.y = element_blank(),
                        strip.text.x.top = element_blank(),
                        strip.background.x = element_blank(),
                        plot.margin = margin(r = 30, b = 5, l = 30, ),
                        legend.position = "none")
  
  # Get order of individuals
  ind_order <- get_ind_order(tibble)
  
  # Pipe Kpop to ggplot
  bottom_plot <- tibble %>%
    mutate(ind = factor(ind, levels = ind_order)) %>% #format data frame to order individuals
    # Begin plot command
    ggplot(aes(x = reorder(ind, population), y = percent, fill = K)) +
    geom_col(position = "fill") +
    theme_classic() +
    facet_grid(~population, space = "free", scales = "free_x") +
    bottom_theme
  
  # return output plot
  return(bottom_plot)
}

#FIXME Read in variables:
# it could take in a fofn list, a fam file, and a pop_order file
# read these in as a list?
fam_loc <- "HoubaraFeb23_noUndulata_hiqual_admixture.fam"
order_loc <- "pop_order_allInds.txt"
q_loc_K2 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K2.rep1.Q"
q_loc_K3 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K3.rep1.Q"
q_loc_K4 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K4.rep1.Q"
q_loc_K5 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K5.rep1.Q"
q_loc_K6 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K6.rep1.Q"
q_loc_K7 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K7.rep1.Q"
q_loc_K8 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K8.rep1.Q"
q_loc_K9 <- "HoubaraFeb23_noUndulata_hiqual_admixture.K9.rep1.Q"

Kpop_list <- list(q_loc_K2, q_loc_K3, q_loc_K4, q_loc_K5, q_loc_K6, q_loc_K7, q_loc_K8, q_loc_K9)

# Parse variables to Kpop objects
output_prefix <- "Asian Houbara AllInds"
pop_levels <- parse_pop_order(order_loc)
pop_labels <- label_translate(pop_levels)



Kpop_list <- map(.x = Kpop_list, .f = parse_q_and_fam, fam_loc = fam_loc)
Kpop_list <- map(.x = Kpop_list, .f = get_long_K_tibble, pop_order = pop_levels, pop_labels = pop_labels)
  # save all Kpop_objects as a .csv:
rlist::list.save(Kpop_list, "Kpop_list_AllInds.rds")
rlist::list.save(Kpop_list, "Kpop_list_AllInds.yaml") # this is an awful yaml file???

# Make mass plots
  # get ind_order for all plots. This will be the same for each value of K
ind_order <- get_ind_order(Kpop_list[[1]]) # doesn't really matter what entry we run it on.
num_plots <- length(Kpop_list)
plot_list <- vector("list", length = num_plots)

plot_list <- map(.x = Kpop_list, .f = get_mass_pop_plot_tilepatch)
rlist::list.save(plot_list, "plot_list_AllInds.rds")

#### CREATE STACKED PLOT WITH PATCHWORK AND TILING ####
#wrap_plots((top / mid) / bottom, ncol = 1) # text case for how patchwork works

K2 <- get_top_plot(Kpop_list[[1]])
K3 <- get_mid_plot(Kpop_list[[2]])
K4 <- get_mid_plot(Kpop_list[[3]])
K5 <- get_mid_plot(Kpop_list[[4]])
K6 <- get_mid_plot(Kpop_list[[5]])
K7 <- get_mid_plot(Kpop_list[[6]])
K8  <- get_mid_plot(Kpop_list[[7]])
K9 <- get_bottom_plot(Kpop_list[[length(Kpop_list)]])
output_prefix = "Asian Houbara (all populations)"


stack_plot_title <- str_c("Admixture: ", output_prefix)
stack_plot_subtitle <- str_c("K2 - K9 ", "")

stacked_plot <- wrap_plots(((((((K2 / K3) / K4) / K5) / K6) / K7) / K8) / K9)
stacked_plot +
  plot_annotation(title = stack_plot_title,
                subtitle = stack_plot_subtitle) & 
  theme(title = element_text(face = "bold", size = 14))


#  
filename <- "AllInds_K2-9_stacked_plot.pdf"
ggsave(
  filename = filename,
  plot = stacked_plot,
  device = "pdf",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300
)
