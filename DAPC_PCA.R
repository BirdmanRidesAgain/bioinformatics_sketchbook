#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# DAPC_PCA.R
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
  message("No .vcf file was provided.")
  stop("Usage: ./DAPC_PCA.R <input.vcf> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  message("No .fam file was provided.")
  stop("Usage: ./DAPC_PCA.R <input.vcf> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==2)
{
  message("No pop_order.txt was provided.")
  stop("Usage: ./DAPC_PCA.R <input.vcf> <input.fam> <pop_order.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==3)
{
  message("No prefix was provided. 'output' used for all outfiles.")
  ARGS[4] = "output"
}


###################################################
### LOAD LIBRARIES AND PRINT PARAMS ###
###################################################
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
message(paste(" Input .fam file:", fam_loc))
order_loc <- ARGS[3]
message(paste(" Input pop_order.txt:", order_loc))
output_prefix <- ARGS[4]
message(paste(" Output prefix:", output_prefix))
message("")

###################################################
### HELPER FUNCTIONS ###
###################################################
### Parsing functions
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
} # this function is designed to exclude individuals without pop. data - superfluous otherwise
create_poplist <- function(matched_ind_pop_df) {
  # This function takes the matched ind/pop pairs and formats them to a list.
  # The list can be used as input to add populations to the genind object.
  # Any missing populations are given as "NA" 
  poplist <- matched_ind_pop_df %>%
    select(population) %>%
    pull(var = population)
  return(poplist)
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

### Clustering functions
# DAPCs
run_dapc <- function(data, num_PCs = 10, num_DAs = 3) {
  # Runs a DAPC, with data as a mandatory input. num_PCs and num_DAs have sensible default values.
  
  # Running the DAPC commands. Your number of PCs is at most argument K-1 (Joshua Thia 2022)
  num_PCs <- length(levels(pop_info$population)) - 1
  num_DAs <- 3 # FIXME : chosen after viewing the F-statistics for successive values of eigenvectors; hardcoded. Modifying this will require modifying the 'colnames' assignment statement in plot_dapc()
  dapc_data <- dapc(data, data$pop, n.pca = num_PCs, n.da = num_DAs)
  return(dapc_data)
}
prep_dapc_ind_coords <- function(dapc_data, data, pop_order, label_order) {
  # This command formats the input PCA data for the plotting command. It is equivalent, albeit simpler than the corresponding DAPC function.
  ind_coords <- as_tibble(dapc_data$ind.coord)
  colnames(ind_coords) <- c("Axis1","Axis2","Axis3")
  ind_coords$Ind = indNames(data)
  ind_coords$Pop <- factor(data$pop, 
                           levels=pop_order, 
                           labels=label_order)
  
  # Sorts the data frame by the population order given above
  #ind_coords <- left_join(data.frame(Pop = label_order), ind_coords, by = "Pop", multiple = "all") # This magically sorts my data frame my the vector given above
  
  return(ind_coords)
}

# PCAs
run_pca <- function(data) {
  # This function masks the PCA running from the user. We're implementing it mostly to avoid titanic data matrices.
  data_fillmissing <- tab(data, NA.method = "mean") # not sure what this does. I believe it averages out missing values.
  
  # now we run the actual pca
  pca_object <- dudi.pca(data_fillmissing, scannf = FALSE, scale = FALSE, nf = 3)
  return(pca_object)
}
prep_pca_ind_coords <- function(pca_data, data, pop_order, label_order) {
  # This command formats the input PCA data for the plotting command. It is equivalent, albeit simpler than the corresponding DAPC function.
  ind_coords <- as_tibble(pca_data$li)
  colnames(ind_coords) <- c("Axis1","Axis2","Axis3")
  ind_coords$Ind = indNames(data)
  ind_coords$Pop <- factor(data$pop, 
                           levels=pop_order, 
                           labels=label_order)
  
  # Sorts the data frame by the population order given above
  #ind_coords <- left_join(data.frame(Pop = pop_order), ind_coords, by = "Pop", multiple = "all") # This magically sorts my data frame my the vector given above
  
  return(ind_coords)
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

### Printing functions
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


###################################################
### READ IN VCF AND FAM ###
###################################################
data <- get_genind_from_vcf(file_loc)
pop_info <- parse_fam_file(fam_loc)

#######################################################
### ADD POPULATION INFORMATION TO VCF/GENIND OBJECT ###
#######################################################
matched_pop_inds <- match_ind_and_pop(data, pop_info)
ordered_poplist <- create_poplist(matched_pop_inds) # this is included to check for NAs
pop(data) <- ordered_poplist # assign individuals to populations in the genind object


#######################################################
### RUN PCA / DAPC ###
#######################################################
dapc_data <- run_dapc(data) # default "num_DAs" hardcoded to 3.
pca_data <- run_pca(data)

#######################################################
### PLOT PCA / DAPC ###
#######################################################
### Create levels from the given populations
pop_order <- parse_pop_order(order_loc)
label_order <- label_translate(pop_order)

# FIXME hardcoded colors. # The order is : "Mongolia", "E_Kazakhstan", "C_Kazakhstan", "C_Uzbekistan", "W_Kazakhstan", "W_Uzbekistan", "N_Iran", "S_Iran", "Israel", "Yemen"
cols <- c("#1F78B4","#A6CEE3","#B2DF8A","#33A02C","#E31A1C","#FB9A99","#6A3D9A","#CAB2D6","#B15928","#FFFF99")

# Prep data frames for plotting
ind_coords_dapc <- prep_dapc_ind_coords(dapc_data, data, pop_order, label_order) # DAPC
ind_coords_pca <- prep_pca_ind_coords(pca_data, data, pop_order, label_order) # PCA

plot_title_dapc <- str_c(output_prefix, " DAPC")
plot_title_pca <- str_c(output_prefix, " PCA")

output_dapc <- plot_cluster_analysis(ind_coords_dapc, cols, plot_title_dapc)
output_pca <- plot_cluster_analysis(ind_coords_pca, cols, plot_title_pca)


###################################################
### SAVE INDIVIDUAL COORDINATES AS CSV ###
###################################################
filename_dapc_coords <- stringr::str_c(output_prefix,"_dapc_ind_coords.csv")
filename_pca_coords <- stringr::str_c(output_prefix,"_pca_ind_coords.csv")

readr::write_csv(ind_coords_dapc, file = filename_dapc_coords)
readr::write_csv(ind_coords_pca, file = filename_pca_coords)


###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
# Create an outfile for all graphic output
project_type <- "DAPC_PCA"
out_dir <- define_output_dir(output_prefix, project_type) # you need to run this, even if the file is already created. We use 'out_dir' as a variable later.

# Save plots in new out directory
filename_dapc_plot <- str_c(output_prefix,"_DAPC",".pdf")
save_plot(output_dapc, filename_dapc_plot)
filename_pca_plot <- str_c(output_prefix,"_PCA",".pdf")
save_plot(output_pca, filename_pca_plot)


message(str_c("Output files written to '", out_dir, "'"))
message("Program end.")
