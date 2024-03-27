#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# DAPC_PCA.R
# KCollier 25 Jan 2024
# This is a command-line-usable script designed to create DAPC and PCA figures from a vcf and a .fam file.
# Some code is not my own; most of the new additions are just to force it to work from command line.
# This is adapted from an earlier script meant to run a single houbara DAPC/PCA. 
# The latest update (26 Mar 2024) edits this behavior to run round of both clustering algs. per 'unknown' individual.

###################################################
### PARSE ARGUMENTS ###
###################################################
default_args <- list(help = FALSE,
                     vcf = "NA",
                     pop = "NA", 
                     classify = FALSE,
                     classified_pop = "unknown",
                     prefix = "output", 
                     pop_order = "NA", 
                     cols = "NA")
args <- R.utils::commandArgs(trailingOnly = TRUE, 
                             asValues = TRUE,
                             defaults = default_args)
print_help <- function() {
  message("\nDAPC_PCA.R, Keiler Collier, 25 Jan 2024, v.0.1.2
  Suggested usage (standard): ./DAPC_PCA.R --args --vcf <input.vcf> --pop <input.txt> --prefix <output_prefix>
  Suggested usage (classify): ./DAPC_PCA.R --args --vcf <input.vcf> --pop <input.txt> --classify --classified_pop <string> --prefix <output_prefix>

  
  STANDARD USE:
  This program runs two clustering alogrithms (PCA and DAPC) and plots them with GGplot2.
    Mandatory Inputs:
      1. A filtered .vcf.
      2. A .tsv with the population of each individual in column 1, and the individual's ID in column 2. These must match the vcf.
        Example formatting is as follows; population and individual ordering is irrelevant:
          pop_1 ind_1
          pop_1 ind_2
          pop_2 ind_3
          pop_2 ind_4
          
      3. A prefix for all outputs. If not defined, this is set to \"output\" by default.
      
    Default Outputs:
      1. High-quality .pngs of the PCA/DAPC plots. These are placed into a subdirectory (<prefix>_plots) titled with your prefix.
      2. One .csv apiece of the PCA and DAPC objects with associated population information.
      
    Optional Arguments:
      --pop_order <pop_order.txt>
        Controls the ordering of your population in the plots; orders populations alphabetically otherwise.
          This argument takes a text file of population names, each on a new line. 
          If the desired order were \"pop_a\", \"pop_b\", \"pop_c\", the text file would be formatted:
            pop_a
            pop_b
            pop_c
          
      --cols <cols.txt>
        Controls the colors used for each population.
          This argument takes a text file of population names, each on a new line. All colors included in the RColorBrewer library are valid inputs.
          Colors will be applied in the order that the populations are ordered - changing --pop_order will result in shuffling colors.
          For the above population order, pop_a would be colored 'red', pop_b 'steelblue', and pop_c 'salmon'
            red
            steelblue
            salmon
      
  CLASSIFY MODE:
  Suggested usage: ./DAPC_PCA.R --args --vcf <input.vcf> --pop <input.txt> --classify --classified_pop <string> --prefix <output_prefix>
    If we have multiple unknown individuals, we may want to produce individual PCA/DAPCs for each one. This is best done when you have a well-characterized reference pool.
    'Classify mode' is an experimental function that takes all individuals in a particular population (\"unknown\",by default), and produces PCAs with all reference inds and one unknown at a time.
      
    Mandatory Inputs:
      1. A filtered .vcf, as above.
      2. A .tsv with the population of each individual, as above. However, some subset of individuals needs to match the string specified by the --classified_pop arugment.
        For the following example, ind_5 and ind_6 are specified:
          pop_1 ind_1
          pop_1 ind_2
          pop_2 ind_3
          pop_2 ind_4
          unknown ind_5
          unknown ind_6
          
      3. A prefix for all outputs. If not defined, this is set to \"output\" by default.
      
    Default Outputs:
      1. n high-quality .pngs of the PCA/DAPC plots in a separate directory (<prefix>_classified_plots), where n is the number of inds in the specified population.
      2. 2*n .csvs of the PCA and DAPC objects with associated population information.
          \n")
  quit()
}
if ( args$help == TRUE ) {
  print_help()
}

# Sanity check for required inputs:
if ( args$vcf == default_args$vcf || args$pop == default_args$pop ) {
  message("VCF or population file missing. See --help for more detailed information.")
  stop("./DAPC_PCA.R --args --vcf <input.vcf> --pop <pop_order.txt> <output_prefix>", call. = FALSE)
} 
# If you didn't add in a col or pop_order file, set those values to NULL internally
if (args$pop_order == "NA") {
  args$pop_order <- NULL
} 
if (args$cols == "NA") {
  args$cols <- NULL
}

###################################################
### LOAD LIBRARIES AND PRINT PARAMS ###
###################################################
message("Loading libraries: ")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(scales))

message("All necessary packages loaded.")
message("")

suppressPackageStartupMessages(library(tidyverse))

### PRINT CRITICAL PARAMS FOR USER:
message(str_c("Input vcf is: ", toString(args$vcf)))
message(str_c("Input population file is: ", toString(args$pop)))
message(str_c("Output prefix is: ", toString(args$prefix)))


###################################################
### HELPER FUNCTIONS ###
###################################################
### Parsing functions
get_genind_from_vcf <- function(vcf) {
  message("Reading in .vcf:")
  data1 <- vcfR::read.vcfR(vcf)
  message("Converting vcf to genind object:")
  data <- vcfR::vcfR2genind(data1)
  message("Genind successfully generated.\n")
  return(data)
}
parse_pop_file <- function(pop_file) {
  ### READ IN FAM FILE ###
  message("\nReading in pop file:")
  pop_cols <- c("population","ind")
  pop_data <- suppressMessages(read_table(file = pop_file, 
                                col_names = pop_cols) %>%
    mutate(population = as.factor(population)))
  message("Population data parsed.\n")
  return(pop_data)
  }
match_ind_and_pop <- function(vcf, parse_fam) {
  # This function takes 2 args:
  # 1. the individuals in a genind object/parsed vcf, and
  # 2. A data frame consisting of "population, ind". This is assumed to come from a parsed .fam file, but you can assemble one yourself.
  # It left-joins them to attempt to find matching populations for all individuals in the genind object. Then, it outputs the relevant data frame.

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
run_DAPC <- function(data, num_PCs = 10, num_DAs = 3) {
  # Running the DAPC commands. Your number of PCs is at most argument K-1 (Joshua Thia 2022)
  num_PCs <- length(levels(pop_info$population)) - 1
  num_DAs <- 3 # FIXME : chosen after viewing the F-statistics for successive values of eigenvectors; hardcoded. Modifying this will require modifying the 'colnames' assignment statement in plot_dapc()
  DAPC_data <- dapc(data, data$pop, n.pca = num_PCs, n.da = num_DAs)
  return(DAPC_data)
}
prep_DAPC_ind_coords <- function(DAPC_data, data, pop_order, label_order) {
  # This command formats the input PCA data for the plotting command. It is equivalent, albeit simpler than the corresponding DAPC function.
  ind_coords <- as_tibble(DAPC_data$ind.coord)
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
run_PCA <- function(data) {
  # This function masks the PCA running from the user. We're implementing it mostly to avoid titanic data matrices.
  data_fillmissing <- tab(data, NA.method = "mean") # not sure what this does. I believe it averages out missing values.

  # now we run the actual pca
  PCA_object <- dudi.pca(data_fillmissing, scannf = FALSE, scale = FALSE, nf = 3)
  return(PCA_object)
}
prep_PCA_ind_coords <- function(PCA_data, data, pop_order, label_order) {
  # This command formats the input PCA data for the plotting command. It is equivalent, albeit simpler than the corresponding DAPC function.
  ind_coords <- as_tibble(PCA_data$li)
  colnames(ind_coords) <- c("Axis1","Axis2","Axis3")
  ind_coords$Ind = indNames(data)
  ind_coords$Pop <- factor(data$pop,
                           levels=pop_order,
                           labels=label_order)

  # Sorts the data frame by the population order given above
  #ind_coords <- left_join(data.frame(Pop = pop_order), ind_coords, by = "Pop", multiple = "all") # This magically sorts my data frame my the vector given above

  return(ind_coords)
}

plot_cluster_analysis <- function(ind_coords, title_string, colors = cols) {
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
define_output_dir <- function(prefix, project_type) {
  # Function to check for the existence of an output directory and create one in WD if it does not exist.
  out_dir <- str_c("./", prefix, "_", project_type, "_figures/")

  # check for the existence of the directory
  if (dir.exists(out_dir) == FALSE) {
    message(str_c("Output directory '",out_dir,"' created"))
    dir.create(out_dir)
  }
  return(out_dir)
}
save_plot <- function(plot, filename, out_dir) {
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
data <- get_genind_from_vcf(args$vcf)
pop_info <- parse_pop_file(args$pop)

#######################################################
### ADD POPULATION INFORMATION TO VCF/GENIND OBJECT ###
#######################################################
matched_pop_inds <- match_ind_and_pop(data, pop_info)
ordered_poplist <- create_poplist(matched_pop_inds) # this is included to check for NAs
pop(data) <- ordered_poplist # assign individuals to populations in the genind object


#######################################################
### RUN SINGULAR PCA / DAPC ###
#######################################################
DAPC_data <- run_DAPC(data) # default "num_DAs" hardcoded to 3.
PCA_data <- run_PCA(data)

#####################
### CLASSIFY MODE ###
#####################
#This entire thing is going to take place in a giant fucking if statement
if (args$classify) { 
  message("Running classify mode:")
  message("Designated population is: ", toString(args$classified_pop))
  
  ######################################
  ### CLASSIFY MODE HELPER FUNCTIONS ###
  ######################################
  subsample_genind <- function(list, data) {
    new_list <- data[-list]
    return(new_list)
  }
  
  ######################################
  ### CLASSIFY MODE MAIN FUNCTION ###
  ######################################
  ### Get named list of indices of 'unknown' inds
  classify_idx <- which(pop(data) == args$classified_pop)

  ### Use for loop to create a list of arrays of all n-1 combinations of classify_idx
    # FIXME - translate to .map function later
  classify_idx_list <- list()
  for (i in 1:length(classify_idx)) {
    permutation <- list(classify_idx[-i])
    classify_idx_list <- append(classify_idx_list, permutation)
  }
  ### Use a second loop to subset data object and assign each subset to a new genind
  data_subsets <- map(.x = classify_idx_list, .f = subsample_genind, data)
  
  ### Run DAPCs/PCAs of all data subsets
  classify_DAPC <- map(.x = data_subsets, .f = run_DAPC)
  classify_PCA <- map(.x = data_subsets, .f = run_PCA)
}

#######################################################
### PLOT PCA / DAPC ###
#######################################################
### FIXME - Everything below here is fucking spaghetti code. Please reformat into a data frame or something useful
  # you shouldn't have ten different vectors sitting around; this is DUMB

### If pop_order is defined as a param, re-order the factor according to that file. Otherwise, leave as-is.
  # FIXME - parse_pop_order and your color read_tsv() call can be one helper function.

if (!is.null(args$pop_order)) {
  pop_order <- parse_pop_order(args$pop_order)
} else { pop_order = levels(pop(data)) }
label_order <- label_translate(pop_order)

### If colors are provided as an input param, read them in here. Otherwise, pick colors with rainbow().
if (!is.null(args$cols)) {
  cols <- read_csv(args$cols, col_names = "colors") 
} else {
  cols <- rainbow(length(pop_order))
}

#Prep data frames for plotting
ind_coords_DAPC <- prep_DAPC_ind_coords(DAPC_data, data, pop_order, label_order) # DAPC
ind_coords_PCA <- prep_PCA_ind_coords(PCA_data, data, pop_order, label_order) # PCA

plot_title_DAPC <- str_c(args$prefix, " DAPC")
plot_title_PCA <- str_c(args$prefix, " PCA")

output_DAPC <- plot_cluster_analysis(ind_coords_DAPC, plot_title_DAPC, cols)
output_PCA <- plot_cluster_analysis(ind_coords_PCA, plot_title_PCA, cols)


# IF CLASSIFY MODE HAS BEEN ENABLED:
if (args$classify) {
  classify_DAPC_coords <- map2(.x = classify_DAPC, .y = data_subsets, .f = prep_DAPC_ind_coords, pop_order, label_order)
  classify_PCA_coords <- map2(.x = classify_PCA, .y = data_subsets, .f = prep_PCA_ind_coords, pop_order, label_order)
  
  # Get the names of all the individuals to generate plot titles
  indnames_unknown <- rownames(data[classify_idx]$tab)
  classify_plot_title_DAPC <- map(.x = indnames_unknown, .f = str_c, "_", args$prefix, "_DAPC")
  classify_plot_title_PCA <- map(.x = indnames_unknown, .f = str_c, "_", args$prefix, "_PCA")
  classify_output_DAPC <- map2(.x = classify_DAPC_coords, .y = classify_plot_title_DAPC, .f = plot_cluster_analysis, cols)
  classify_output_PCA <- map2(.x = classify_PCA_coords, .y = classify_plot_title_PCA, .f = plot_cluster_analysis, cols)
}


###################################################
### SAVE INDIVIDUAL COORDINATES AS CSV ###
###################################################
filename_DAPC_coords <- stringr::str_c(args$prefix,"_DAPC_ind_coords.csv")
filename_PCA_coords <- stringr::str_c(args$prefix,"_PCA_ind_coords.csv")

readr::write_csv(ind_coords_DAPC, file = filename_DAPC_coords)
readr::write_csv(ind_coords_PCA, file = filename_PCA_coords)

### IF CLASSIFY MODE ENABLED
# helper function to save classify_csvs
save_classify_csvs <- function(coord, indname, prefix, plot_type) {
  filename_coord <- str_c(prefix, "_", indname, "_", plot_type, ".csv")
  write_csv(coord, file = filename_coord)
}
if (args$classify) {
  map2(.x = classify_DAPC_coords, .y = indnames_unknown, .f = save_classify_csvs, prefix = args$prefix, plot_type = "DAPC")
  map2(.x = classify_PCA_coords, .y = indnames_unknown, .f = save_classify_csvs, prefix = args$prefix, plot_type = "PCA")
}


###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
# Create an outfile for all graphic output
project_type <- "DAPC_PCA"
out_dir <- define_output_dir(args$prefix, project_type) # you need to run this, even if the file is already created. We use 'out_dir' as a variable later.

# Save plots in new out directory
filename_DAPC_plot <- str_c(args$prefix,"_DAPC",".pdf")
save_plot(output_DAPC, filename_DAPC_plot, out_dir)
filename_PCA_plot <- str_c(args$prefix,"_PCA",".pdf")
save_plot(output_PCA, filename_PCA_plot, out_dir)

### IF CLASSIFY MODE ENABLED
if (args$classify) {
  save_classify_plot <- function(plot, indname, prefix, out_dir, plot_type) {
    filestring <- str_c(prefix, "_", indname, "_", plot_type, "_classify.pdf")
    save_plot(plot, filestring, out_dir)
  }
  map2(.x = classify_output_DAPC, .y = classify_plot_title_DAPC, .f = save_classify_plot, args$prefix, out_dir, "DAPC")
  map2(.x = classify_output_PCA, .y = classify_plot_title_PCA, .f = save_classify_plot, args$prefix, out_dir, "PCA")
}

message(str_c("Output files written to '", out_dir, "'"))
message("Program end.")
