#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# plot_pairwise_fst.R
# KCollier 7 Feb 2024
# This is a command-line-usable script designed to create pairwise FST plots from a .vcf file and a .fam file
# Some code is not my own; most of the new additions are just to force it to work from command line.

# CID is 12738464

###################################################
### PARSE ARGUMENTS ###
###################################################
# Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  stop("Usage: ./plot_pairwise_fst.R <input.vcf> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==3)
{
  message("No prefix was provided. 'output' used for all outfiles.")
  ARGS[4] <- "output"
}
# Load libraries
library(tidyverse)
library(adegenet)
library(hierfstat)
library(vcfR)

message("All necessary packages loaded.")
message("")

message("Parameters interpreted as:")
file_loc <- ARGS[1]
message(paste(" Input .vcf: ", file_loc))

fam_loc <- ARGS[2]
message(paste(" Input .fam: ", fam_loc))

order_loc <- ARGS[3]
message(paste(" Input pop_order.txt: ", order_loc))

output_prefix <- ARGS[4]
message(paste(" Output prefix: ", output_prefix))
message("")

###################################################
### HELPER FUNCTIONS ###
###################################################
# Parsing functions
get_genind_from_vcf <- function(vcf) {
  message("Reading in .vcf:")
  data1 <- vcfR::read.vcfR(vcf)
  message("Converting vcf to genind object:")
  data <- vcfR::vcfR2genind(data1)
  return(data)
}
parse_fam_file <- function(fam_file) {
  ### READ IN FAM FILE ###
  message("Reading in .fam file:")
  
  fam_cols <- c("population","ind")
  fam_data <- readr::read_table(file = fam_file, 
                                col_names = fam_cols) %>%
    mutate(population = as.factor(population))
  return(fam_data)
}
match_ind_and_pop <- function(vcf, parse_fam) {
  # This function takes 2 args: 
  # 1. the individuals in a genind object/parsed vcf, and
  # 2. A data frame consisting of "ind, population". This is assumed to come from a parsed .fam file, but you can assemble one yourself.
  # It left-joins them to attempt to find matching populations for all individuals in the genind object. Then, it outputs the relevant data frame. Warnings are supplied if there are any NAs in the data frame.
  
  data_inds <- as_tibble(adegenet::indNames(vcf)) %>% # Get list of individuals from the vcf
    rename(ind = value)
  ind_pop_matches <- left_join(data_inds, parse_fam, by = "ind") # left join the two dataframes
  return(ind_pop_matches)
}
create_poplist <- function(matched_ind_pop_df) {
  # This function takes the matched ind/pop pairs and formats them to a list.
  # The list can be used as input to add populations to the genind object.
  # Any missing populations are given as "NA" 
  poplist <- matched_ind_pop_df %>%
    select(population) %>%
    pull(var = population)
  return(poplist)
}
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

# Plotting functions
get_fst_tibble <- function(fst_object, pop_order) {
  # generate two matrices according to population order
  fst_mat <- as.matrix(fst_object)
  fst_mat_order <- fst_mat[unlist(pop_order), unlist(pop_order)]
  
  # convert matrices to tibble
  ind <- which(upper.tri(fst_mat_order), arr.ind = TRUE)
  fst_tibble <- tibble(Pop1 = dimnames(fst_mat_order)[[2]][ind[, 2]],
                       Pop2 = dimnames(fst_mat_order)[[1]][ind[, 1]],
                       Fst = fst_mat_order[ind])
  
  # Get new population levels
  Pop1_levels <- pop_order[-c(1)]
  Pop2_levels <- pop_order[-c(length(pop_order))]
  
  # Get new population labels
  Pop1_labels <- label_translate(Pop1_levels)
  Pop2_labels <- label_translate(Pop2_levels)
  
  # Re-level populations and convert Fst scores
  fst_tibble <- fst_tibble %>%
    mutate(Pop1 = factor(Pop1,
                         levels = Pop1_levels,
                         labels = Pop1_labels),
           Pop2 = factor(Pop2, 
                         levels = pop_levels,
                         labels = pop_labels))
  fst_tibble$Fst[fst_tibble$Fst < 0] = 0 # set all scores below zero to zero
  
  return(fst_tibble)
}
plot_fst_tibble <- function(fst_tibble) {
  # display structure of tibble
  fst_tibble %>% str
  
  # Get requires plotting arguments
  fst_label <- expression(italic("F")[ST])
  max_fst <- max(fst_tibble$Fst)
  mid <- max_fst / 2
  
  fst_theme <- theme(axis.text = element_text(color = "black",
                                              size = 15,
                                              face = "bold"),
                     axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0),
                     axis.title = element_blank(),
                     panel.grid = element_blank(),
                     panel.background = element_blank(),
                     legend.position = "right",
                     legend.title = element_text(size = 14, face = "bold"),
                     legend.text = element_text(size = 10))
  
  # Begin plotting command
  fst_plot <- ggplot(fst_tibble, aes(x = Pop1, y = Pop2, fill = Fst)) +
    geom_tile(color = "black") +
    geom_text(aes(label = Fst), color = "black", size = 9) +
    scale_fill_gradient2(low = "blue", mid = "pink", high = "red",
                         midpoint = mid,
                         name = fst_label,
                         limits = c(0, max_fst),
                         breaks = c(0, 0.05, 0.10, 0.15)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0), position = "right") +
    theme_classic() +
    fst_theme
  
  return(fst_plot)
}
save_plot <- function(plot = "plot", out_dir = "./", output_prefix) {
  # Get filename from output prefix
  output_prefix_clean <- str_replace_all(output_prefix, pattern = " ", replacement = "_")
  filename <- str_c(output_prefix_clean,".pdf")
  filepath <- str_c(out_dir,filename)
  
  # check for existence of output directory
  if(dir.exists(out_dir) == FALSE) {
    message(str_c("Output directory '", out_dir, "' created"))
    dir.create(out_dir)
  } else if (out_dir == "./") {
    message(str_c("No output directory given. Working directory assumed."))
  }
  
  # save plot to output directory
  ggsave(
    filename = filepath,
    plot = plot,
    device = "pdf",
    width = 500,
    height = 250,
    units = "mm",
    dpi = 300
  )
  if (out_dir == "./") {
    message(str_c("Output plot '", filename, "' written to '", getwd(),"/",out_dir,"'"))
  } else {
    message(str_c("Output plot '", filename, "' written to '", getwd(),"'"))
  }
}

###################################################
### MAIN ###
###################################################
### PARSE AND WRANGLE DATA ###
# get pop_level/label objects.
pop_levels <- get_list_from_1_col_txt(order_loc)
pop_labels <- label_translate(pop_levels)

# get genind and population info
data <- get_genind_from_vcf(file_loc)
pop_info <- parse_fam_file(fam_loc)
# add population info to genind object
matched_pop_inds <- match_ind_and_pop(data, pop_info)
ordered_poplist <- create_poplist(matched_pop_inds)
pop(data) <- ordered_poplist # assign individuals to populations in the genind object



### CALCULATE PAIRWISE FST SCORES ###
message("Beginning to calculate pairwise FST. This will take a long time.")
fst <- genet.dist(data, method = "WC84") %>%
  round(digits = 3)
saveRDS(fst)
message("FST table saved to working directory.")

# now we generate and plot the fst matrix
fst_tibble <- get_fst_tibble(fst, pop_levels)
newplot <- plot_fst_tibble(fst_tibble)
save_plot(newplot, getwd(), output_prefix)

message("Program terminating.")
