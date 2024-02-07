# load libraries
library(tidyverse)

# Helper
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

### read in all four CSVs
DAPC_all <- read.csv("AllInds_DAPC_ind_coords.csv")
PCA_all <- read.csv("AllInds_PCA_ind_coords.csv")
DAPC_noYemen <- read.csv("NoYemen_DAPC_ind_coords.csv")
PCA_noYemen <- read.csv("NoYemen_PCA_ind_coords.csv")
DAPC_noYemen_noMongolia <- read.csv("csv_output/noMongolia_DAPC_coords.csv")
PCA_noYemen_noMongolia <- read.csv("csv_output/noMongolia_PCA_coords.csv")
DAPC_noYemen_noMongolia_noEKazakh <- read.csv("csv_output/noEKazakh_DAPC_coords.csv")
PCA_noYemen_noMongolia_noEKazakh  <- read.csv("csv_output/noEKazakh_PCA_coords.csv")

### Plot theme:
label_theme <- theme(axis.title.x=element_text(size=14),
                     axis.title.y=element_text(size=14), 
                     plot.title=element_text(hjust=0.5, 
                                             size=16, 
                                             face="bold"))

#### SPLIT 0 - ALL INDS ####
# AllInds levels + labels:
pop_order_allInds <- c("Yemen",  "Israel","S_Iran",    "N_Iran",    "W_Uzbekistan",   "W_Kazakhstan",   "C_Uzbekistan",     "C_Kazakhstan",       "E_Kazakhstan",   "Mongolia")
label_order_allInds <- c("Yemen","Israel","South Iran","North Iran","West Uzbekistan","West Kazakhstan","Central Uzbekistan","Central Kazakhstan","East Kazakhstan","Mongolia")
cols_allInds <- c("#FFFF99","#B15928","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A","#A6CEE3","#1F78B4","#B2DF8A","#33A02C")

### Plot DAPC_all
title_string = "Asian Houbara DAPC (All Populations)"
DAPC_all$Pop <-factor(DAPC_all$Pop,
                      levels = pop_order_allInds,
                      labels = label_order_allInds)
DAPC_all_plot <- ggplot(DAPC_all, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
# Add colors:
  scale_fill_manual(values = cols_allInds) +
  scale_color_manual(values = cols_allInds) +
# Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
DAPC_all_plot

### Plot PCA_all
title_string = "Asian Houbara PCA (All Populations)"
PCA_all$Pop <-factor(PCA_all$Pop,
                      levels = pop_order_allInds,
                      labels = label_order_allInds)
PCA_all_plot <- ggplot(PCA_all, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_allInds) +
  scale_color_manual(values = cols_allInds) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
PCA_all_plot


#### SPLIT 1 - NO YEMEN ####
# NoYemen levels + labels:
pop_order_noYemen <- c("Israel",  "S_Iran",    "N_Iran",    "W_Uzbekistan",   "W_Kazakhstan",   "C_Uzbekistan",     "C_Kazakhstan",       "E_Kazakhstan",   "Mongolia")
label_order_noYemen <- c("Israel","South Iran","North Iran","West Uzbekistan","West Kazakhstan","Central Uzbekistan","Central Kazakhstan","East Kazakhstan","Mongolia")
cols_noYemen <- c("#B15928","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A","#A6CEE3","#1F78B4","#B2DF8A","#33A02C")

### Plot DAPC_NoYemen
title_string = "Asian Houbara DAPC (No Yemen)"
DAPC_noYemen$Pop <-factor(DAPC_noYemen$Pop,
                     levels = pop_order_noYemen,
                     labels = label_order_noYemen)
levels(DAPC_noYemen$Pop)

DAPC_noYemen_plot <- ggplot(DAPC_noYemen, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen) +
  scale_color_manual(values = cols_noYemen) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
DAPC_noYemen_plot


### Plot PCA_NoYemen
title_string = "Asian Houbara PCA (No Yemen)"
PCA_noYemen$Pop <-factor(PCA_noYemen$Pop,
                          levels = pop_order_noYemen,
                          labels = label_order_noYemen)
levels(PCA_noYemen$Pop)
levels(PCA_noYemen$Pop)
PCA_noYemen_plot <- ggplot(PCA_noYemen, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen) +
  scale_color_manual(values = cols_noYemen) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
PCA_noYemen_plot

#### SPLIT 2 - NO MONGOLIA ####
# NoYemen_NoMongolia levels + labels:
pop_order_noYemen_NoMongolia <- c("Israel",  "S_Iran",    "N_Iran",    "W_Uzbekistan",   "W_Kazakhstan",   "C_Uzbekistan",     "C_Kazakhstan",       "E_Kazakhstan")
label_order_noYemen_NoMongolia <- c("Israel","South Iran","North Iran","West Uzbekistan","West Kazakhstan","Central Uzbekistan","Central Kazakhstan","East Kazakhstan")
cols_noYemen_NoMongolia <- c("#B15928","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A","#A6CEE3","#1F78B4","#B2DF8A")

### Plot DAPC_NoYemen_NoMongolia
title_string = "Asian Houbara DAPC (No Yemen, No Mongolia)"
DAPC_noYemen_noMongolia$Pop <-factor(DAPC_noYemen_noMongolia$Pop,
                          levels = pop_order_noYemen_NoMongolia,
                          labels = label_order_noYemen_NoMongolia)
levels(DAPC_noYemen_noMongolia$Pop)

DAPC_noYemen_noMongolia_plot <- ggplot(DAPC_noYemen_noMongolia, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen_NoMongolia) +
  scale_color_manual(values = cols_noYemen_NoMongolia) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
DAPC_noYemen_noMongolia_plot


### Plot PCA_NoYemen_noMongolia
title_string = "Asian Houbara PCA (No Yemen, No Mongolia)"
PCA_noYemen_noMongolia$Pop <-factor(PCA_noYemen_noMongolia$Pop,
                         levels = pop_order_noYemen_NoMongolia,
                         labels = label_order_noYemen_NoMongolia)
levels(PCA_noYemen_noMongolia$Pop)
levels(PCA_noYemen_noMongolia$Pop)
PCA_noYemen_noMongolia_plot <- ggplot(PCA_noYemen_noMongolia, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen_NoMongolia) +
  scale_color_manual(values = cols_noYemen_NoMongolia) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
PCA_noYemen_noMongolia_plot

#### SPLIT 3 - NO E KAZAKH ####
# NoYemen_NoMongolia levels + labels:
pop_order_noYemen_NoMongolia_noEKazakh <- c("Israel",  "S_Iran",    "N_Iran",    "W_Uzbekistan",   "W_Kazakhstan",   "C_Uzbekistan",     "C_Kazakhstan")
label_order_noYemen_NoMongolia_noEKazakh <- c("Israel","South Iran","North Iran","West Uzbekistan","West Kazakhstan","Central Uzbekistan","Central Kazakhstan")
cols_noYemen_NoMongolia_noEKazakh <- c("#B15928","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A","#A6CEE3","#1F78B4")

### Plot DAPC_NoYemen_NoMongolia_noEKazakh
title_string = "Asian Houbara DAPC (No Yemen, No Mongolia, no EKazakh)"


DAPC_noYemen_noMongolia_noEKazakh_plot <- ggplot(DAPC_noYemen_noMongolia_noEKazakh, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen_NoMongolia_noEKazakh) +
  scale_color_manual(values = cols_noYemen_NoMongolia_noEKazakh) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
DAPC_noYemen_noMongolia_noEKazakh_plot


### Plot PCA_NoYemen_noMongolia_noEKazakh
title_string = "Asian Houbara PCA (No Yemen, No Mongolia, no EKazakh)"

PCA_noYemen_noMongolia_noEKazakh_plot <- ggplot(PCA_noYemen_noMongolia_noEKazakh, aes(x = Axis1, y = Axis2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = Pop), shape = 21, size = 3)+
  # Add colors:
  scale_fill_manual(values = cols_noYemen_NoMongolia_noEKazakh) +
  scale_color_manual(values = cols_noYemen_NoMongolia_noEKazakh) +
  # Add custom labels
  label_theme +
  xlab("Axis 1") +
  ylab("Axis 2") +
  ggtitle(title_string)
PCA_noYemen_noMongolia_noEKazakh_plot 


### Print plots to file system
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
save_plot <- function(plot, filename, out_dir = "./") {
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

save_plot(DAPC_all_plot, "DAPC_all_plot.pdf")
save_plot(PCA_all_plot, "PCA_all_plot.pdf")
save_plot(DAPC_noYemen_plot, "DAPC_noYemen_plot.pdf")
save_plot(PCA_noYemen_plot, "PCA_noYemen_plot.pdf")
save_plot(DAPC_noYemen_noMongolia_plot, "DAPC_noYemen_noMongolia_plot.pdf")
save_plot(PCA_noYemen_noMongolia_plot, "PCA_noYemen_noMongolia_plot.pdf")
save_plot(DAPC_noYemen_noMongolia_noEKazakh_plot, "DAPC_noYemen_noMongolia_noEKazakh_plot.pdf")
save_plot(PCA_noYemen_noMongolia_noEKazakh_plot, "PCA_noYemen_noMongolia_noEKazakh_plot.pdf")
