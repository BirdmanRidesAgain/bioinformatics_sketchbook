#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# DAPC.R
# KCollier 25 Jan 2024
# This is a command-line-usable script designed to create DAPC and PCA figures from a vcf and a .fam file.
# Some code is not my own; most of the new additions are just to force it to work from command line.

###################################################
### PARSE ARGUMENTS ###
###################################################
# Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  stop("Usage: ./DAPC.R <input.vcf> <input.fam> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==2)
{
  message("No prefix was provided. 'output' used for all outfiles.")
  ARGS[2] <- "output"
}
# Load libraries
library(tidyverse)
library(adegenet)
library(scales)
library(vcfR)

message("All necessary packages loaded.")
message("")

message("Parameters interpreted as:")
file_loc <- ARGS[1]
message(paste(" Input .vcf: ", file_loc))

fam_loc <- ARGS[2]
message(paste(" Input .fam: ", fam_loc))

output_prefix <- ARGS[3]
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
remove_NA_from_genind <- function(matched_pop_inds, data) {
  # Produce list of individuals with no population assigned
  NA_list <- matched_pop_inds %>%
    filter(is.na(population)) %>%
    pull(var = ind)
  
  # filter those individuals from the genind object
  for (i in 1:length(NA_list)) {
    data <- data[indNames(data) != NA_list[[i]]]
    message(str_c(NA_list[[i]]), " removed")
  }
  return(data)
}

# Clustering functions
run_dapc <- function(data, num_PCs = 10, num_DAs = 3) {
  # Runs a DAPC, with data as a mandatory input. num_PCs and num_DAs have sensible default values.
  
  # Running the DAPC commands. Your number of PCs is at most argument K-1 (Joshua Thia 2022)
  num_PCs <- length(levels(pop_info$population)) - 1
  num_DAs <- 3 # FIXME : chosen after viewing the F-statistics for successive values of eigenvectors; hardcoded. Modifying this will require modifying the 'colnames' assignment statement in plot_dapc()
  dapc_data <- dapc(data, data$pop, n.pca = num_PCs, n.da = num_DAs)
  return(dapc_data)
}
plot_cluster_analysis <- function(ind_coords, colors = cols, title_string) {
  # Function used to plot clustering analyses
  # Define theme and title variables:
  label_theme <- theme(axis.title.x=element_text(size=14),
                       axis.title.y=element_text(size=14), 
                       plot.title=element_text(hjust=0.5, 
                                               size=16, 
                                               face="bold"))
  
  # Generate the plot:
  plot <- ggplot(ind_coords, aes(x = Axis1, y = Axis2)) +
    # Set guidelines marking the X and Y intercepts
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    
    # Add points
    geom_point(aes(fill = Pop), shape = 21, size = 3) +
    
    # Add coloring
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    # Add custom labels:
    label_theme +
    xlab("Axis 1") +
    ylab("Axis 2") +
    ggtitle(title_string)
  
  return(plot) # returns clustered plot
}


###################################################
### READ IN VCF AND FAM ###
###################################################
library(tidyverse)
library(adegenet)
library(scales)
library(vcfR)

file_loc <- "Houbara_NoYemen.recode.vcf"
data <- get_genind_from_vcf(file_loc)
fam_loc <- "Houbara_NoYemen.fam"
pop_info <- parse_fam_file(fam_loc)
output_prefix <- "No Yemen"
#######################################################
### ADD POPULATION INFORMATION TO VCF/GENIND OBJECT ###
#######################################################
matched_pop_inds <- match_ind_and_pop(data, pop_info)
ordered_poplist <- create_poplist(matched_pop_inds)
pop(data) <- ordered_poplist # assign individuals to populations in the genind object



#######################################################
### RUN PCA / DAPC ###
#######################################################
dapc_data <- run_dapc(data) # default "num_DAs" hardcoded to 3.


#######################################################
### PLOT PCA / DAPC ###
#######################################################
# FIXME hardcoded pop. order
pop_order <- c("Mongolia", "E_Kazakhstan", "C_Kazakhstan", "C_Uzbekistan", "W_Kazakhstan", "W_Uzbekistan", "N_Iran", "S_Iran", "Israel")
label_order <- c('Israel','South Iran','North Iran','West Uzbekistan','Central Uzbekistan','West Kazakhstan','Central Kazakhstan','East Kazakhstan','Mongolia')
# FIXME hardcoded colors. # The order is : "Yemen", "Israel", "S Iran", "N Iran", "W Uzbekistan", "W Kazakhstan", "C Uzbekistan", "C Kazakhstan", "E Kazakhstan", "Mongolia"
cols <- c("#B15928","#CAB2D6","#6A3D9A","#FB9A99","#E31A1C","#33A02C","#B2DF8A","#A6CEE3","#1F78B4")

# Prep data frames for plotting
ind_coords_dapc <- prep_dapc_ind_coords(dapc_data, data, pop_order, label_order) # DAPC
plot_title_dapc <- str_c(output_prefix, " DAPC")
output_dapc <- plot_cluster_analysis(ind_coords_dapc, cols, plot_title_dapc)


###################################################
### SAVE INDIVIDUAL COORDINATES AS CSV ###
###################################################
filename_dapc_coords <- stringr::str_c(output_prefix,"_dapc_ind_coords.csv")
readr::write_csv(ind_coords_dapc, file = filename_dapc_coords)


###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
# Printing functions
define_output_dir <- function(output_prefix, project_type) {
  # Function to check for the existence of an output directory and create one in WD if it does not exist.
  out_dir <- str_c("./", output_prefix, "_", project_type, "_figures/")
  
  # check for the existence of the directory
  if (dir.exists(out_dir) == FALSE) {
    message(str_c("Output directory '",out_dir,"' created"))
    dir.create(out_dir)
  }
  return(out_dir)
}
save_plot <- function(plot, filename) {
  filetype = "pdf"
  ggsave(
    filename = str_c(out_dir,filename),
    plot = plot,
    device = filetype,
    width = 500,
    height = 250,
    units = "mm",
    dpi = 300
  )
}
# Create an outfile for all graphic output
project_type <- "DAPC_PCA"
out_dir <- define_output_dir(output_prefix, project_type) # you need to run this, even if the file is already created. We use 'out_dir' as a variable later.

# Save plots in new out directory
filename_pca_plot <- str_c(output_prefix,"_PCA",".pdf")
save_plot(output_pca, filename_pca_plot)
filename_dapc_plot <- str_c(output_prefix,"_DAPC",".pdf")
save_plot(output_dapc, filename_dapc_plot)


message(str_c("Output files written to '", out_dir, "'"))
message("Program end.")
