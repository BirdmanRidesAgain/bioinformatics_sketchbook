#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# visualize_stacked_admixture.R
# KCollier 2 Feb 2024

# This is a command-line-usable script that takes a list of Q scores from Admixture,
# and uses them to plot all values of K for a single population of animals.
# Tidyverse and patchwork are dependencies.

###################################################
### PARSE ARGUMENTS ###
###################################################
# Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  message("No list of Q-files was provided.")
  stop("Usage: ./visualize_stacked_admixture.R <input_Qs.txt> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  message("No .fam file was provided.")
  stop("Usage: ./visualize_stacked_admixture.R <input_Qs.txt> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==2)
{
  message("No pop_order.txt was provided.")
  stop("Usage: ./visualize_stacked_admixture.R <input_Qs.txt> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==3)
{
  ARGS[4] = "output"
}

###################################################
### LOAD PACKAGES AND DEFINE VARIABLES ###
###################################################
# Dependencies
library(tidyverse)
library(patchwork)
library(rlist)

# Parse arguments
message("Parameters interpreted as:")
qs_loc <- ARGS[1]
message(paste(" Input list of Q-files:", qs_loc))
fam_loc <- ARGS[2]
message(paste(" Input .fam file:", fam_loc))
order_loc <- ARGS[3]
message(paste(" Input pop_order.txt:", order_loc))
output_prefix <- ARGS[4]
message(paste(" Output prefix:", output_prefix))
message("")

###################################################
### HELPER FUNCTIONS - NON-PLOTTING ###
###################################################
# Parsing functions:
get_list_from_1_col_txt <- function(list_loc) {
  message(stringr::str_c("Reading in: ", list_loc))
  message("")
  list_data <- readr::read_table(list_loc, col_names = "X1")
  list_data <- as.list(list_data$X1)
  
  return(list_data)
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
parse_q_and_fam <- function(q_loc, fam_loc) {
  ### READ IN Q FILE ###
  message("Reading in Q-file:")
  message("")
  
  q_data <- readr::read_table(file = q_loc, col_names = FALSE)
  num_K <- ncol(q_data)
  
  ### READ IN FAM FILE ###
  suppressWarnings({
    message("Reading in .fam file:")
    message("")
    
    fam_cols <- c("population","ind")
    fam_data <- readr::read_table(file = fam_loc, 
                                  col_names = fam_cols)
  })
  ### BIND Q AND FAM TIBBLES BY COLUMN ###
  Kpop <- bind_cols(fam_data, q_data) %>%
    arrange(population)
  return(Kpop)
}

# Data wrangling functions
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

###################################################
### MAIN ###
###################################################
### PARSE AND WRANGLE DATA ###
# get pop_level/label objects. Necessary for formatting Kpop objects.
pop_levels <- get_list_from_1_col_txt(order_loc)
pop_labels <- label_translate(pop_levels)


# Parse input Q files to Kpop objects 
Kpop_list <- get_list_from_1_col_txt(qs_loc) # get list of all Kpops
Kpop_list <- map(.x = Kpop_list, .f = parse_q_and_fam, fam_loc = fam_loc)
Kpop_list <- map(.x = Kpop_list, .f = get_long_K_tibble, pop_order = pop_levels, pop_labels = pop_labels)

# Get min and max K values. This assumes you aren't plotting K=1, because you aren't stupid.
min_k <- 2 # min_K will always be 2
max_k <- length(Kpop_list) + 1 # Because K=1 is not plotted, max_k is always 1 higher than the list size.

# save all Kpop_objects as a .csv:
message("Outputting formatted Kpop data to disc: ")
Kpop_list_filename <- str_c(output_prefix, "_K", min_k,"-", max_k, "_Kpop_data")
rlist::list.save(Kpop_list, str_c(Kpop_list_filename,".rds"))
rlist::list.save(Kpop_list, str_c(Kpop_list_filename,".yaml")) # this is an awful yaml file???
message(str_c("  ",Kpop_list_filename,".rds"," written to ", getwd()))
message(str_c("  ",Kpop_list_filename,".yaml"," written to ", getwd()))



###################################################
### HELPER FUNCTIONS - PLOTTING ###
###################################################
# Plotting functions
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
get_positional_plot <- function(tibble, pos_theme) {
  # Get order of individuals
  ind_order <- get_ind_order(tibble)
  
  # Pipe Kpop to ggplot
  pos_plot <- tibble %>%
    mutate(ind = factor(ind, levels = ind_order)) %>% #format data frame to order individuals
    # Begin plot command
    ggplot(aes(x = reorder(ind, population), y = percent, fill = K)) +
    geom_col(position = "fill") +
    theme_classic() +
    facet_grid(~population, space = "free", scales = "free_x") +
    pos_theme +
    # Shared theme components:
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(r = 30, b = 0, l = 30))
  
  # return output plot
  return(pos_plot)
}

K2 <- get_positional_plot(Kpop_list[[1]], top_theme)
K2
###################################################
### PLOT KPOP OBJECTS ###
###################################################
message("Now plotting Kpop objects: ")

### Define positional themes ###
# FIXME : these themes can be compressed to remove redundant lines

top_theme <- theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.ticks.length = unit(0, "pt"), #length of tick marks
                   axis.line.x.bottom = element_blank())

mid_theme <- theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.ticks.length = unit(0, "pt"), #length of tick marks
                   axis.line.x.bottom = element_blank(),
                   strip.text.x.top = element_blank(),
                   strip.background.x = element_blank())

bottom_theme <- theme(axis.text.y = element_blank(),
                      axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
                      axis.ticks.y = element_blank(),
                      strip.text.x.top = element_blank(),
                      strip.background.x = element_blank())

### Plot the first and last population with 'top' and 'bottom' themes
K2 <- get_positional_plot(Kpop_list[[1]], top_theme)
KN <- get_positional_plot(Kpop_list[[length(Kpop_list)]], bottom_theme)

### Plot all K values with the 'middle' theme
Kplot_middle_list <- map(.x = Kpop_list, .f = get_positional_plot, pos_theme = mid_theme)

  # Get the start and end indices of the desired 'middle' plots
start_middle_index <- 2
end_middle_index <- length(Kpop_list) - 1

# Use patchwork and the subset of the Kplot_middle_list var to stack the plots
stacked_plot <- wrap_plots(K2 / (Kplot_middle_list)[start_middle_index:end_middle_index] / KN)

  # Set titles:
stack_plot_title <- str_c(output_prefix)
stack_plot_subtitle <- str_c("K",min_k, " - K",max_k)

stacked_plot_title <- stacked_plot +
  plot_annotation(title = stack_plot_title,
                  subtitle = stack_plot_subtitle) & 
  theme(title = element_text(face = "bold", size = 14))

###################################################
### SAVE PLOTS TO FILESYSTEM ###
###################################################
filename <- str_replace(str_c(output_prefix, "_K",min_k,"-",max_k, "_stacked_plot.pdf"), " ", "_")
ggsave(
  filename = filename,
  plot = stacked_plot,
  device = "pdf",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300
)

# Plot version without a title
filename <- str_replace(str_c(output_prefix, "_K",min_k,"-",max_k, "_stacked_plot_title.pdf"), " ", "_")
ggsave(
  filename = filename,
  plot = stacked_plot_title,
  device = "pdf",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300
)

message(str_c("Plots written to", getwd()))
message("Program finished.")
