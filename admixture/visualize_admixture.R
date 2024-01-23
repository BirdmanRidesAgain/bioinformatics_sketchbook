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
  stop("Usage: ./visualize_contig_lengths.R <input.Q> <input.fam> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  message("No .fam file was provided. All individuals will be treated as one population.")
} else if(length(ARGS)==2)
{
  ARGS[3] = "output"
}

### WE NEED TO ADD THIS AS AN ARGUMENT SOMEHOW? ###
pop_levels <- c("Yemen",
               "Israel",
               "S_Iran",
               "N_Iran",
               "W_Uzbekistan",
               "W_Kazakhstan",
               "C_Uzbekistan",
               "C_Kazakhstan",
               "E_Kazakhstan",
               "Mongolia")

###################################################
### CHECKS FOR ALL NECESSARY PACKAGES ###
###################################################
# Check if tidyverse is installed. Install if it is not present.
#if (!requireNamespace("tidyverse", quietly = FALSE))
#{
#  install.packages("tidyverse")
#}
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(mapmixture)

message("All necessary packages loaded.")
message("")

message("Parameters interpreted as:")
message(paste(" Input Q-file:", ARGS[1]))
message(paste(" Input .fam file:", ARGS[2]))
message(paste(" Output prefix:", ARGS[3]))
message("")

###################################################
### HELPER FUNCTIONS ###
###################################################
# K tibbles
### WRITE DATAFRAME FOR EACH VALUE OF K TO K_LIST ###
get_long_K_tibble <- function(tibble, k_val, pop_order) {
  new_list <- vector(mode = 'list', length = k_val)
  
  # loop to fill all slots in 'new_list', using get_K_tibble() helper function
  for (i in 1:k_val) {
    colname <- str_c("X", i)
    new_list[[i]] <- get_K_tibble(tibble, eval(colname), i)
  }
  
  # Reshape 'new_list' into long data frame
  long_K_tibble <- bind_rows(new_list) #%>%
  long_K_tibble <- long_K_tibble %>%
    mutate(K = as.factor(K), # Recode <int> K to act like factor to ease plotting later
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
get_ind_pop_plots <- function(tibble) {
  num_pops <- length(levels(tibble$population))
  pop_list <- vector(mode = 'list', length = num_pops)
  
  # write population-specific tibbles to 'pop_list'
  for (i in 1:num_pops) {
    pop_list[[i]] <- get_pop_tibble(tibble, i)
  }
  
  # iterate over list to plot all the values
  pop_plot_list <- vector(mode = 'list', length = num_pops)
  for (i in 1:num_pops) {
    pop_plot_list[[i]] <- get_ind_pop_plot(pop_list[[i]], eval(levels(tibble$population)[i]))
  }
  
  # return the list of individual plots
  return(pop_plot_list)
}
get_ind_pop_plot <- function(tibble, pop_string) {
  #Get title/subtitle strings
  plot_title <- str_replace(string = pop_string , pattern = "_", replacement = ". ")
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

###################################################
### CODE BELOW THIS POINT IS WIP/TEST-ONLY ###
###################################################
library(tidyverse)

################################################################
### READ IN AND BIND Q+FAM FILES TO CREATE INITIAL DATAFRAME ###
################################################################
### READ IN Q FILE ###
message("Reading in Q-file:")
message("")

table_loc <- ARGS[1]
q_data <- readr::read_table(file = table_loc, col_names = FALSE)
num_K <- ncol(q_data)

### READ IN FAM FILE ###
message("Reading in .fam file:")
message("")

fam_loc <- "HoubaraFeb23_noUndulata_hiqual_admixture.fam"
fam_cols <- c("population","ind")
fam_data <- readr::read_table(file = fam_loc, 
                              col_names = fam_cols)

### BIND Q AND FAM TIBBLES BY COLUMN ###
Kpop <- bind_cols(fam_data,q_data) %>%
  arrange(population)

#########################################################
### REFORMAT KPOP TO HAVE ONE OBSERVATION PER K VALUE ###
#########################################################
### CREATE K_LIST - ONE ELEMENT FOR EVERY VALUE OF K ###
Kpop <- get_long_K_tibble(Kpop, num_K, pop_levels) # test of new function


# Prior code we knew that worked - remove once the script has successfully run
#for (i in 1:num_K) {
#  colname <- str_c("X",i)
#  K_list[[i]] <- get_K_tibble(Kpop, eval(colname), i)
#}
#
#Kpop <- bind_rows(K_list) %>%
#  mutate(K = as.factor(K), # Recode <int> K to act like factor
#         population = factor(as.factor(Kpop$population), # Recode and relevel <chr> population as factor
#                             levels = pop_levels))

###################################################
### PLOT STACKED BAR OUTPUT ###
###################################################
# Set themes
ind_pop_theme <- theme(
  title = element_text(size = 18, face = "bold"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16))

mass_pop_theme <- theme(
  title = element_text(size = 18, face = "bold"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  legend.position = "none")

### Get list of all 
test_list <- get_pop_tibbles(Kpop)
test_list
test_list[[1]]

# Prior code we know worked - remove once the script has run successfully
#num_pops <- length(levels(Kpop$population))
#pop_list <- vector(mode = 'list', length = num_pops)
#
#for (i in 1:num_pops) {
#  pop_list[[i]] <- get_pop_tibble(Kpop, i)
#}
#
## iterate over list to plot all the values
#pop_plot_list <- vector(mode = 'list', length = num_pops)
#for (i in 1:num_pops) {
#  pop_plot_list[[i]] <- get_ind_pop_plot(pop_list[[i]], eval(levels(Kpop$population)[i]))
#}

### OUTPUT A PLOT FOR ALL POPULATIONS ###
ggplot(Kpop, aes(x = ind, y = percent, fill = K)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  mass_pop_theme




###################################################
### SAVE PLOT TO LOCAL FILESYSTEM ###
###################################################
ggsave(
  str_c(prefix_filename,"_descending_contigs.tiff"),
  plot = plot1,
  device = "tiff",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300,
)
message("Program end.")


