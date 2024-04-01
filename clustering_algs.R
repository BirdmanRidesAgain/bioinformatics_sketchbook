#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# clustering_algs.R
# KCollier 25 Jan 2024
# This is a command-line-usable script designed to create DAPC and PCA figures from a vcf and a .fam file.
# This is adapted from an earlier script meant to run a single houbara DAPC/PCA. 
# The latest update (26 Mar 2024) edits this behavior to run round of both clustering algs. per 'unknown' individual.

###################################################
### PARSE ARGUMENTS AND SANITY CHECKS ###
###################################################
default_args <- list(help = FALSE,
                     vcf = NULL,
                     pop = NULL, 
                     classify = FALSE,
                     classified_pop = "unknown",
                     pca = FALSE,
                     mclust_num = 0,
                     prefix = "output", 
                     pop_order = NULL, 
                     cols = NULL,
                     filetype = "pdf")
args <- R.utils::commandArgs(
  asValues = TRUE,
  defaults = default_args)

# Sanity check for --help
# FIXME - please update help function
print_help <- function() {
  sink(stdout(), type = "message")
  message("\nDAPC_PCA.R, Keiler Collier, 25 Jan 2024, v.0.1.2
  Suggested usage (standard): ./clustering_algs.R --vcf <input.vcf> --pop <input.txt> --prefix <output_prefix>
  Suggested usage (classify): ./clustering_algs.R --vcf <input.vcf> --pop <input.txt> --classify --classified_pop <string> --prefix <output_prefix>
  
  ALL OPTIONS:
                     --vcf : <input vcf; mandatory>
                     --pop : <input tsv file linking pops with inds; mandatory> 
                     --classify : <runs classify mode when set>
                     --classified_pop : <population in --pop to run individual DAPCs from; \"unknown\" by default>
                     --pca : <runs PCAs in addition to DAPCs when set>
                     --mclust_num : <runs mclust a second time, using the value provided as an a priori number of populations (K)>
                     --prefix : <prefix used for all output materials. \"output\" by default> 
                     --pop_order : <text file supplying the order of populations for GGplot. Optional.>
                     --cols : <text file supplying the colors used for each population in GGplot. Optional.>
                     --filetype : <output format of plots. Options of 'jpeg', 'tiff', 'pdf'; default 'pdf'>
  
  QUICK START:
  STANDARD (SINGULAR) MODE:
  This program runs one discriminant analysis of principal components (DAPC) analysis on all individuals, and another on all reference individuals - ie, all individuals not included in the 'classified_pop' option. It then plots the coordinates with GGplot2, estimates group identity with Mclust, and saves the coordinates and plots to your filesyste. You can also optionally run a vanilla principal components analysis (PCA) by adding the --pca flag to the command-line invocation.
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
  sink(NULL, type = "message")
  quit()
}
if (args$help) {
  print_help()
}

# Sanity check for required inputs:
if ( is.null(args$vcf) || is.null(args$pop) ) {
  message("VCF or population file missing. See --help for more detailed information.")
  stop("./clustering_algs.R --vcf <input.vcf> --pop <pop_order.txt> --prefix <output_prefix>", call. = FALSE)
}

###################################################
### LOAD LIBRARIES AND PRINT PARAMS ###
###################################################
### CHECK AND LOAD LIBS
check_installed_libs <- function(libs) {
  message("Checking for uninstalled libraries...")
  new_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
  if(length(new_libs)) {
    message("Missing packages found. Installing.")
    install.packages(new_libs)
  }
  message("")
}
lib_load <- function(libs) {
  # Function to smooth out library loading process
  message("Loading libraries:")
  for (i in 1:length(libs)) { 
    suppressPackageStartupMessages(library(libs[i], character.only = TRUE))
    message(paste("\t",toString(libs[i]),"successfully loaded"))
  }
  message("All packages loaded\n")
}

libs <- c("tidyverse", "adegenet", "vcfR", "scales", "mclust") # list of libs
# Load libraries present in system
check_installed_libs(libs)
lib_load(libs)

### PRINT CRITICAL PARAMS FOR USER:
  # FIXME - should only write the end-user useful params.
print_params <- function(args) {
  message("User-defined parameters:")
  
  for (i in 1:length(args)) {
    name <- names(args[i])
    if (!name == "") {
      message(str_c("\t",name, ": ", toString(args[i])))
      }
    }
  message("")
}
print_params(args)

###################################################
### PREP GENIND OBJECT FOR ANALYSES ###
###################################################
message("BEGIN ANALYSIS:")
### HELPER FUNCTIONS
get_genind_from_vcf <- function(vcf) {
  message("READING IN VCF:")
  data <- read.vcfR(vcf) %>%
    vcfR2genind()
  message("Genind successfully generated.\n")
  return(data)
}
parse_pop_file <- function(pop_file) {
  ### READ IN FAM FILE ###
  message("\nParsing pop file:")
  pop_cols <- c("population","ind")
  pop_data <- suppressMessages(read_table(file = pop_file, 
                                          col_names = pop_cols) %>%
                                 mutate(population = as.factor(population)))
  message("Population data parsed.\n")
  return(pop_data)
}
get_ordered_pop_vector <- function(genind, pop_info) {
  # Gets df of individuals
  inds <- as_tibble(indNames(genind)) %>%
    rename(ind = value)
  
  # Associates population data with inds, and gets ordered vect. of populations
  poplist <- left_join(inds, pop_info, by = "ind") %>%
    select(population) %>%
    pull(var = population)
  return(poplist)
}
parse_pop_order <- function(order_loc) {
  message("Reading in population order:")
  message("")
  pop_levels_data <- read_table(order_loc, col_names = "population")
  pop_levels <- pop_levels_data$population
  return(pop_levels)
}

### PARSE INPUT FILES
genind <- get_genind_from_vcf(args$vcf)
pop_info <- parse_pop_file(args$pop)
### ADD POP INFO TO GENIND
ordered_popvect <- get_ordered_pop_vector(genind, pop_info)
pop(genind) <- ordered_popvect # assign individuals to populations in the genind object

  #get pop_order
if (!is.null(args$pop_order)) {
  pop_order <- parse_pop_order(args$pop_order)
} else { pop_order = levels(pop(genind)) }

#get colors
if (!is.null(args$cols)) {
  cols <- as.vector(read.csv(args$cols, header = FALSE))
  cols <- unlist(unname(cols))
} else { 
  cols <- rainbow(length(pop_order))
}


#################################
### ANALYSIS HELPER FUNCTIONS ###
#################################
### CLUSTERING ALGORITHMS
prep_ind_coords <- function(clust_data, genind = genind, pop_order = pop_order, clust_type = "DAPC") {
  # This command formats the input clustering alg. data into a usable data frame
  if (clust_type == "DAPC") {
    coord <- as_tibble(clust_data$ind.coord)
  }
  else { coord <- as_tibble(clust_data$li) }

  
  colnames(coord) <- c("Axis1","Axis2","Axis3")
  coord$Ind = indNames(genind)
  coord$Pop <- factor(genind$pop,
                           levels=pop_order)
  
  # Sorts the data frame by the population order given above
  #ind_coords <- left_join(data.frame(Pop = label_order), ind_coords, by = "Pop", multiple = "all") # This magically sorts my data frame my the vector given above
  
  return(coord)
}
run_DAPC <- function(genind, num_PCs = 10, num_DAs = 3, pop_order = pop_order) {
  # Running the DAPC commands. Your number of PCs is at most argument K-1 (Joshua Thia 2022)
  num_PCs <- length(levels(pop_info$population)) - 1
  num_DAs <- 3 # FIXME : chosen after viewing the F-statistics for successive values of eigenvectors; hardcoded. Modifying this will require modifying the 'colnames' assignment statement in plot_dapc()
  DAPC_data <- dapc(genind, genind$pop, n.pca = num_PCs, n.da = num_DAs)
  
  ind_coords <- prep_ind_coords(clust_data = DAPC_data, genind = genind, pop_order = pop_order, clust_type = "DAPC")
  return(ind_coords)
}
run_PCA <- function(genind, pop_order = pop_order) {
  # This function masks the PCA running from the user. We're implementing it mostly to avoid titanic genind matrices.
  genind_fillmissing <- tab(genind, NA.method = "mean") # not sure what this does. I believe it averages out missing values.
  
  # now we run the actual pca
  PCA_data <- dudi.pca(genind_fillmissing, scannf = FALSE, scale = FALSE, nf = 3)
  ind_coords <- prep_ind_coords(clust_data = PCA_data, genind = genind, pop_order = pop_order, clust_type = "PCA")
  return(ind_coords)
}

### SUBSAMPLING
subsample_genind <- function(list, genind) {
  new_genind <- genind[-list]
  return(new_genind)
}

### MCLUST
run_mclust <- function(coord, num_clust = 0) {
  if (num_clust == 0) {
    clust <- Mclust(coord)
  }
  else {
    clust <- Mclust(coord, num_clust)
  }
  # We want $classification only
  filt_clust <- clust %>%
    keep_at(at = c("classification")) %>%
  return(filt_clust)
}
mclust_wrapper <- function(coord, num_clust = 0) {
  # Get input tibble to work on
  tibble <- coord %>% 
    select(starts_with("Axis"))
  
  out_list <- list() # initialize list
  out_list$optim_clust <- run_mclust(tibble) #Run mclust to get the optimal number of clusters
  coord <- cbind(coord,out_list$optim_clust) %>% #bind and rename mclust results to coord
    rename(optim_classification = classification)
  
  if (num_clust > 0) { #If a-priori value was supplied, run with that
    prefix <- str_c("a_priori_K",num_clust,"_")
    rename_to_var <- function(colname) { return(colname) }
    
    out_list$a_priori_clust <- run_mclust(tibble, num_clust)
    coord <- cbind(coord, out_list$a_priori_clust) %>%
      rename_with(.fn = ~paste0(prefix,.), starts_with("classification"))
  }
  return(coord)
}

### PLOTTING
get_plot_title <- function(prefix, plot_type = "DAPC") {
  plot_title <- str_c(prefix, " ", plot_type)
  return(plot_title)
}
plot_cluster_analysis <- function(coord, plot_title = plot_title, cols = cols) {
  # Function used to plot clustering analyses
  # Define theme and title variables:
  label_theme <- theme(axis.title.x=element_text(size=14),
                       axis.title.y=element_text(size=14),
                       plot.title=element_text(hjust=0.5,
                                               size=16,
                                               face="bold"))
  # Generate the plot:
  plot <- ggplot(coord, aes(x = Axis1, y = Axis2)) +
    # Set guidelines marking the X and Y intercepts
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    # Add points
    geom_point(aes(fill = Pop), shape = 21, size = 3) +
    # Add coloring
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    # Add custom labels:
    label_theme +
    xlab("Axis 1") +
    ylab("Axis 2") +
    ggtitle(plot_title)
  
  return(plot) # returns clustered plot
}


##############################
### SINGULAR ANALYSIS MODE ###
##############################
if (singular<-TRUE) {
  message("Starting singular mode:")
  
  # FIXME - this is a function to run csv and plot simultaneously. It incorrectly names your plots, though - fix it.
  get_csv_plot_pair <- function(genind, mclust_num = 0, prefix, clust_type = "DAPC", plot_title = "output", cols = cols, pop_order = pop_order) {
    out_list <- list()
    out_list$csv <- if (clust_type == "DAPC") { run_DAPC(genind) } else { run_PCA(genind) } %>%
      mclust_wrapper(num_clust = args$mclust_num)
    out_list$plot <- plot_cluster_analysis(out_list$csv, plot_title = plot_title, cols = cols)
    return(out_list)
  }
  get_csv_plot_pair <- function(genind, mclust_num = 0, prefix, clust_type, plot_title = "output", cols = cols, pop_order = pop_order) {
    out_list <- list()
    out_list$csv <- if (clust_type == "DAPC") { run_DAPC(genind) } else { run_PCA(genind) } %>%
      mclust_wrapper(num_clust = args$mclust_num)
    out_list$plot <- plot_cluster_analysis(out_list$csv, plot_title = plot_title, cols = cols)
  
    
    return(out_list)
  }

  # FIXME END; BEGIN VERIFIED CODE AGAIN
  
  # Initialize output list
  out_singular <- list()
  
  out_singular$DAPC$all_inds$csv <- run_DAPC(genind) %>%
    mclust_wrapper(num_clust = args$mclust_num)
  plot_title <- get_plot_title(args$prefix, "DAPC (all individuals)")
  out_singular$DAPC$all_inds$plot <- plot_cluster_analysis(out_singular$DAPC$all_inds$csv, plot_title = plot_title, cols)
  
  # Populate DAPC content for only reference individuals
  classify_idx <- which(pop(genind) == args$classified_pop)
    # recycled code from above
  ref_genind <- subsample_genind(classify_idx, genind)
  out_singular$DAPC$ref_inds$csv <- run_DAPC(ref_genind) %>%
    mclust_wrapper(num_clust = args$mclust_num)
  plot_title <- get_plot_title(args$prefix, "DAPC (reference individuals)")
  out_singular$DAPC$ref_inds$plot <- plot_cluster_analysis(out_singular$DAPC$ref_inds$csv, plot_title = plot_title, cols)
  
  # If --pca flag set, populate PCA content
  if (args$pca) {
    message("\t--pca flag set: Running pcas")
    # PCA for all_inds
    out_singular$PCA$all_inds$csv <- run_PCA(genind) %>%
      mclust_wrapper(num_clust = args$mclust_num)
    plot_title <- get_plot_title(args$prefix, "PCA (all individuals)")
    out_singular$PCA$all_inds$plot <- plot_cluster_analysis(out_singular$PCA$all_inds$csv, plot_title = plot_title, cols)
    
    #PCA for ref_only
    out_singular$PCA$all_inds$csv <- run_PCA(ref_genind) %>%
      mclust_wrapper(num_clust = args$mclust_num)
    plot_title <- get_plot_title(args$prefix, "PCA (all individuals)")
    out_singular$PCA$ref_inds$plot <- plot_cluster_analysis(out_singular$PCA$all_inds$csv, plot_title = plot_title, cols)
  }
  message("End singular mode.\n")
}

##############################
### CLASSIFY MODE ###
##############################
if (args$classify) {
  message("Starting classify mode:")
  ### HELPER FUNCTIONS
  subsample_unknowns <- function(classified_pop, genind) {
    # We need to produce the subsampled sets of the genind object
    genind_list <- list() # initialize list
    # create subsets of of all (n-1) combinations of all n 'unknown' inds
    classify_idx <- which(pop(genind) == classified_pop)
    for (i in 1:length(classify_idx)) {
      subsamp <- classify_idx[-i]
      genind_list <- append(genind_list, subsample_genind(subsamp, genind))
    }
    return(genind_list)
  }
  
  # Initialize output list
  out_classify <- list()
  classify_genind <- subsample_unknowns(args$classified_pop, genind)
  csv_list <- as.list(1:length(classify_genind)) # this list will be directly appended to 'out_classify$DAPC$csv'

  ### Run DAPCs over all subsets - gets csvs
  classify_idx <- which(pop(genind) == args$classified_pop) # indices of all the inds in 'unknown'
  indnames_unknown <- rownames(genind[classify_idx]$tab) # names of inds at the indices above. Used to rename the output list
  
  for (i in 1:length(classify_genind)) {
    message(str_c("\tAnalyzing ",indnames_unknown[i]), " (",i,"/",length(classify_genind),")")
    iteration <- classify_genind[[i]]
    item <- as_tibble(run_DAPC(iteration, pop_order))
    csv_list[[i]] <- item
  }
  message("")
  # FIXME - maybe rework this section, b/c redundancy between 'csv_list' and 'out_classify$DAPC$csv'
  names(csv_list) <- indnames_unknown
  out_classify$DAPC$csv <- csv_list

  ### Plot all subsets (all items in 'out_classify$DAPC$csv' and/or 'csv_list')
  plot_list <- as.list(1:length(classify_genind))
  
  for (i in 1:length(classify_genind)) {
    message(str_c("\tPlotting ",indnames_unknown[i]), " (",i,"/",length(classify_genind),")")
    plot_title <- get_plot_title(args$prefix, plot_type = str_c(indnames_unknown[i]," DAPC"))
    plot <- plot_cluster_analysis(csv_list[[i]], plot_title = plot_title, cols = cols)
    plot_list[[i]] <- plot
  }
  message("")
  names(plot_list) <- indnames_unknown
  out_classify$DAPC$plot <- plot_list
  
  ### IF PCA FLAG IS SET, RUN PCAs
  if (args$pca) {
    message("--pca flag set: Running pcas")
    
    # Run PCAs
    csv_list <- as.list(1:length(classify_genind)) # re-initialize object to clear it
    for (i in 1:length(classify_genind)) {
      message(str_c("\tAnalyzing ",indnames_unknown[i]), " (",i,"/",length(classify_genind),")")
      iteration <- classify_genind[[i]]
      item <- as_tibble(run_PCA(iteration, pop_order))
      csv_list[[i]] <- item
    }
    message("")
    names(csv_list) <- indnames_unknown
    out_classify$PCA$csv <- csv_list
    
    ### Plot all subsets (all items in 'out_classify$PCA$csv' and/or 'csv_list')
    plot_list <- as.list(1:length(classify_genind))
    
    for (i in 1:length(classify_genind)) {
      message(str_c("\tPlotting ",indnames_unknown[i]), " (",i,"/",length(classify_genind),")")
      plot_title <- get_plot_title(args$prefix, plot_type = str_c(indnames_unknown[i]," PCA"))
      plot <- plot_cluster_analysis(csv_list[[i]], plot_title = plot_title, cols = cols)
      plot_list[[i]] <- plot
    }
    message("")
    names(plot_list) <- indnames_unknown
    out_classify$PCA$plot <- plot_list
  }

  message("End classify mode.\n")
}


##############################
### SAVE OUTPUTS ###
##############################
### HELPER FUNCTIONS
define_output_dir <- function(prefix = args$prefix, suffix, basedir = "./") {
  # Function to check for the existence of an output directory and create one in WD if it does not exist.
  # check to see if basedir has a backslash at the end of it
  if (!endsWith(basedir, "/")) { # add one if it's not there
    basedir <- str_c(basedir, "/") 
  }
  out_dir <- str_c(basedir, prefix, "_", suffix, "/")
  
  # check for the existence of the directory
  if (dir.exists(out_dir) == FALSE) {
    message(str_c("'", out_dir, "' created"))
    dir.create(out_dir)
  }
  else { message(str_c("'", out_dir,"' already created")) }
  return(out_dir)
}
define_output_dir_rec <- function(name = "out", wd = getwd()) {
  # if 'name' has a '/' in it, break the string down
  if (str_detect(name, "/")) {
    str_array <- strsplit(name,"/")[[1]]
    last_dir <- str_array[length(str_array)]
    
    prior_dirs <- str_replace(name, str_c("/",last_dir), "")
    # convert_prior_dirs to a string. Make recursive call to 'define_output_dir_rec'
    define_output_dir_rec(prior_dirs, wd)
  }
  
  # base case: create directory and message name
  out_dir <- str_c(getwd(), "/", name)
  if (dir.exists(out_dir) == FALSE) {
    dir.create(out_dir)
    message(str_c("'",out_dir, "' created"))
  }
  else { message(str_c("'", out_dir,"' already exists"))}
  # return out_dir on end
  return(out_dir)
}

save_plot <- function(plot, filename, filetype, out_dir) {
  filetype = filetype
  ggsave(
    filename = file.path(out_dir,filename),
    plot = plot,
    device = filetype,
    width = 500,
    height = 250,
    units = "mm",
    dpi = 300
  )
}

### CREATING DIRECTORY STRUCTURE
### FIXME - make this such that you use your fancy new recursive call structure
out_dir_base <- define_output_dir(suffix = "clustering_algs_out")
  out_dir_singular <- define_output_dir(args$prefix, suffix = "singular", out_dir_base)
    out_dir_singular_plot <- define_output_dir(args$prefix, suffix = "singular_plots", out_dir_singular)
    out_dir_singular_csv <- define_output_dir(args$prefix, suffix = "singular_csv", out_dir_singular)
 
  if (args$classify) {
    out_dir_classify <- define_output_dir(args$prefix, suffix = "classify", out_dir_base)
      out_dir_classify_plot <- define_output_dir(args$prefix, suffix = "classify_plots", out_dir_classify)
      out_dir_classify_csv <- define_output_dir(args$prefix, suffix = "classify_csv", out_dir_classify)
  }


### SAVE CSV AND PLOT OUTPUT TO THE APPROPRIATE DIRECTORIES
  # FIXME - same here; code can iterate through a list of filenames/plots
  
  ### Plot singular items
  clust_type <- "DAPC"
  message("Saving singular plots and coordinates:")
  # CSV
  coordname <- str_c(args$prefix, "_all_inds_", clust_type, ".csv")
  filestring <- str_c(out_dir_singular_csv, coordname)
  write.csv(out_singular$DAPC$all_inds$csv, filestring)
  #PLOT
  plotname <- str_c(args$prefix, "_all_inds_", clust_type, ".", args$filetype)
  save_plot(out_singular$DAPC$all_inds$plot, plotname, args$filetype, out_dir_singular_plot)

  # CSV
  coordname <- str_c(args$prefix, "_ref_inds_", clust_type, ".csv")
  filestring <- str_c(out_dir_singular_csv, coordname)
  write.csv(out_singular$DAPC$ref_inds$csv, filestring)
  #PLOT
  plotname <- str_c(args$prefix, "_ref_inds_", clust_type, ".", args$filetype)
  save_plot(out_singular$DAPC$ref_inds$plot, plotname, args$filetype, out_dir_singular_plot)
  
  #PCA
  if (args$pca) {
    clust_type <- "PCA"
    #CSV
    coordname <- str_c(args$prefix, "_all_inds_", clust_type, ".csv")
    filestring <- str_c(out_dir_singular_csv, coordname)
    write.csv(out_singular$PCA$all_inds$csv, filestring)
    #PLOT
      plotname <- str_c(args$prefix, "_all_inds_", clust_type, ".", args$filetype)
      save_plot(out_singular$PCA$all_inds$plot, plotname, args$filetype, out_dir_singular_plot)
      
      #CSV
      coordname <- str_c(args$prefix, "_ref_inds_", clust_type, ".csv")
      filestring <- str_c(out_dir_singular_csv, coordname)
      write.csv(out_singular$PCA$ref_inds$csv, filestring)
      #PLOT
      plotname <- str_c(args$prefix, "_ref_inds_", clust_type, ".", args$filetype)
      save_plot(out_singular$PCA$ref_inds$plot, plotname, args$filetype, out_dir_singular_plot)
    }
    
  ### Plot classify items
  if (args$classify) {
    message("Saving classify plots and coordinates:")
    
    # HELPER FUNCTIONS FOR SAVING CLASSIFY OBJECTS
    save_classify_csvs <- function(coord, indname, prefix, clust_type, out_dir) {
      filename_coord <- str_c(prefix, "_", indname, "_", clust_type, ".csv")
      write.csv(coord, file.path(out_dir, filename_coord))
    }
    save_classify_plot <- function(plot, filename, out_dir) {
      save_plot(plot, filestring, out_dir = out_dir, filetype = args$filetype)
    }
    get_classify_plotname <- function(indname, prefix, clust_type, filetype) {
      str_c(args$prefix, "_", indname, "_", clust_type, ".", args$filetype)
    }

    # SAVE DAPC COORDS
    clust_type <- "DAPC"
    map2(.x = out_classify$DAPC$csv, .y = indnames_unknown, .f = save_classify_csvs, prefix = args$prefix, clust_type = clust_type, out_dir = out_dir_classify_csv)
    
    # SAVE DAPC PLOTS
    ### FIND A WAY TO SAVE THE DAPC PLOTS
      # Get list of filenames for each plot
    classify_plot_title_DAPC <- purrr::map(.x = indnames_unknown, .f = get_classify_plotname, clust_type = clust_type) 
    i = 0
    for (i in 1:length(classify_genind)) {
      filename <- classify_plot_title_DAPC[[i]]
      plot <- out_classify$DAPC$plot[[i]]
      save_plot(plot = plot, filename = filename, out_dir = out_dir_classify_plot, filetype = args$filetype)
    }
    
    
    
    # SAVE PCA OBJECTS
    if(args$pca) {
      # Save coord objects
      clust_type <- "PCA"
      map2(.x = out_classify$PCA$csv, .y = indnames_unknown, .f = save_classify_csvs, prefix = args$prefix, clust_type = clust_type, out_dir = out_dir_classify_csv)
      
      # Save plot objects
      classify_plot_title_PCA <- purrr::map(.x = indnames_unknown, .f = get_classify_plotname, clust_type = clust_type) 
      i = 0
      for (i in 1:length(classify_genind)) {
        filename <- classify_plot_title_PCA[[i]]
        plot <- out_classify$PCAC$plot[[i]]
        save_plot(plot = plot, filename = filename, out_dir = out_dir_classify_plot, filetype = args$filetype)
      }
    }
    message(str_c("Output written to '", out_dir_classify, "'"))
  }
  
  message("Program end.")
