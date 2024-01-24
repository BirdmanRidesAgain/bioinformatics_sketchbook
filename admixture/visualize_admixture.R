#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# visualize_admixture.R
# KCollier 19 Jan 2024

# This is a command-line-usable script that takes Q scores from Admixture and uses them to plot a stacked bargraph.
# tidyverse is a dependency (specifically, ggplot2 and stringr)

###################################################
### PARSE ARGUMENTS ###
###################################################
# Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  message("No .bam file was provided.")
  stop("Usage: ./visualize_contig_lengths.R <input.Q> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  message("No .fam file was provided.")
  stop("Usage: ./visualize_contig_lengths.R <input.Q> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==2)
{
  message("No pop_order.txt was provided.")
  stop("Usage: ./visualize_contig_lengths.R <input.Q> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==3)
{
  ARGS[4] = "output"
}


###################################################
### CHECKS FOR ALL NECESSARY PACKAGES ###
###################################################
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)

message("Parameters interpreted as:")
q_loc <- ARGS[1]
message(paste(" Input Q-file:", q_loc))
fam_loc <- ARGS[2]
message(paste(" Input .fam file:", fam_loc))
order_loc <- ARGS[3]
message(paste(" Input pop_order.txt:", order_loc))
output_prefix <- ARGS[4]
message(paste(" Output prefix:", output_prefix))
message("")

###################################################
### HELPER FUNCTIONS ###
###################################################
# K tibbles
### WRITE DATAFRAME FOR EACH VALUE OF K TO K_LIST ###
parse_pop_order <- function(order_loc) {
  message("Reading in population order:")
  message("")
  pop_levels_data <- readr::read_table(order_loc, col_names = "population")
  pop_levels <- pop_levels_data$population
  return(pop_levels)
}
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
  
get_long_K_tibble <- function(tibble, k_val, pop_order) {
  ### CREATE K_LIST - ONE ELEMENT FOR EVERY VALUE OF K ###
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
    mutate(K = factor(as.factor(long_K_tibble$K), levels = seq(from = num_K, to = 1, by = -1)), # recode and relevel K to act like a factor
      population = factor(as.factor(long_K_tibble$population), # Recode <chr> population as factor
                               levels = pop_order)) # Relevel <fct> population to correctly order populations in plot
  return(long_K_tibble)
}
get_K_tibble <- function(tibble, colname, k_val) {
# This function pulls out the values for a single K, along with population/individual columns
# It then writes them to a new tibble and returns that tibble
  single_K <- tibble %>%
    select(population, ind, colname) %>%
    rename(percent = colname) %>%
    mutate(K = k_val) %>%
    select(population, ind, K, percent)
  return(single_K)
}

# population tibbles
get_pop_tibble <- function(tibble, pop_index) {
# This function creates a new population-specific tibble for the given data frame and index.
  # It is sensitive to the sorting of your files, so it's not recommended to run outside of iterating through all pops.
  single_pop <- tibble %>%
    filter(population == eval(levels(tibble$population)[pop_index]))
  return(single_pop)
}

# plotting functions
get_pop_plots <- function(tibble, type_of_plot) {
  num_pops <- length(levels(tibble$population))
  pop_list <- vector(mode = 'list', length = num_pops)
  
  # write population-specific tibbles to 'pop_list'
  for (i in 1:num_pops) {
    pop_list[[i]] <- get_pop_tibble(tibble, i)
  }
  
  # iterate over list to plot all the values
  pop_plot_list <- vector(mode = 'list', length = num_pops)
  
  for (i in 1:num_pops) {
    if (type_of_plot == "ind") {
      pop_plot_list[[i]] <- get_ind_pop_plot(pop_list[[i]], eval(levels(tibble$population)[i]))
    } else if (type_of_plot == "mass") {
      pop_plot_list[[i]] <- get_mass_pop_plot(pop_list[[i]], eval(levels(tibble$population)[i]))
    } else {
      message("Cannot determine if the plot type is \"ind\" or \"mass\". Reverting to \"ind\".")
      pop_plot_list[[i]] <- get_ind_pop_plot(pop_list[[i]], eval(levels(tibble$population)[i]))
    }
  }

  # return the list of individual plots
  return(pop_plot_list)
}
get_ind_pop_plot <- function(tibble, pop_string) {
  
  #Define theme
  ind_pop_theme <- theme(
    title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))
  
  #Get title/subtitle strings
  plot_title <- str_replace(string = pop_string, pattern = "_", replacement = ". ")
  plot_subtitle <- str_c("K = ", num_K)
  
  #Plot:
  ind_pop_plot <- ggplot(tibble, aes(x = ind, y = percent, fill = K)) +
    labs(title = plot_title,
         subtitle = plot_subtitle) +
    xlab(label = "Individuals") +
    ylab(label = "Percent identity by K") +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) + 
    theme_classic() +
    ind_pop_theme
  
  #Return output plot
  return(ind_pop_plot)
}
get_mass_pop_plot <- function(tibble, pop_string) {
  # Identical to 'get_ind_pop_plot', except 'mass_pop_theme' is used and no subtitle is included
  # Define theme:
  mass_pop_theme <- theme(
    title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
    legend.position = "none")
  
  # Get title string
  #plot_title <- str_replace(string = pop_string, pattern = "_", replacement = ". ")
  # Plot:
  mass_pop_plot <- ggplot(tibble, aes(x = ind, y = percent, fill = K)) +
    #labs(title = plot_title) +
    xlab(label = "") +
    ylab(label = "") +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) + 
    theme_classic() + 
    facet_grid(~population) +
    mass_pop_theme
  
  #Return output plot
  return(mass_pop_plot)
}

#####################################################################
### PARSE INPUTS AND BIND Q+FAM FILES TO CREATE INITIAL DATAFRAME ###
#####################################################################
pop_levels <- parse_pop_order(order_loc)
num_pops <- length(pop_levels)

Kpop1 <- parse_q_and_fam(q_loc, fam_loc)
num_K <- ncol(Kpop1) -2

#########################################################
### REFORMAT KPOP TO HAVE ONE OBSERVATION PER K VALUE ###
#########################################################
Kpop <- get_long_K_tibble(Kpop1, num_K, pop_levels)






###################################################
### PLOT STACKED BAR OUTPUT ###
###################################################
### Get list of all individual plots
ind_pop_list <- vector(mode = 'list', length = num_pops)
ind_pop_list <- get_pop_plots(Kpop, "ind")

###### SECTION BELOW IS HARDCODED: NEEDS TO BE FIXED 
### Patchwork mass plot:

### Get list of all mass plots
mass_pop_list <- vector(mode = 'list', length = num_pops)
mass_pop_list <- get_pop_plots(Kpop, "mass")
mass_pop_list[[2]]


### OUTPUT A PLOT FOR ALL POPULATIONS ###
# 1. PATCHWORK APPROACH:
  # HARDCODED - FIXME
patch_plot_title <- str_replace(string = output_prefix, pattern = "_", replacement = " ")
patch_plot_subtitle <- str_c("Admixture: K = ", num_K)

# Generate patchwork plot
mass_patchwork_plot <-   
  mass_pop_list[[1]] + 
  mass_pop_list[[2]] + 
  mass_pop_list[[3]] + 
  mass_pop_list[[4]] +
  mass_pop_list[[5]] +
  mass_pop_list[[6]] +
  mass_pop_list[[7]] +
  mass_pop_list[[8]] +
  mass_pop_list[[9]] +
  mass_pop_list[[10]] +
  plot_annotation(title = patch_plot_title,
                  subtitle = patch_plot_subtitle,
                  tag_levels = 'A') & theme(text = element_text(face = "bold", size = 14))


###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
### Individual plots:
for (i in 1:num_pops) { 
  ind_plotname = str_c(pop_levels[i],"_admixture_K",num_K,".pdf")
  ggsave(
    filename = ind_plotname,
    plot = ind_pop_list[[i]],
    device = "pdf",
    width = 500,
    height = 250,
    units = "mm",
    dpi = 300
  )
}
### Patchwork mass plot:
patch_plotname = str_c(output_prefix,"_allpops_admixture_K",num_K,".pdf")
ggsave(
  filename = patch_plotname,
  plot = mass_patchwork_plot,
  device = "pdf",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300
)

message("Program end.")
