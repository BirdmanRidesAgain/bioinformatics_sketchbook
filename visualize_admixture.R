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

message("All necessary packages loaded.")
message("")

message("Parameters interpreted as:")
message(paste(" Input Q-file:", ARGS[1]))
message(paste(" Input .fam file:", ARGS[2]))
message(paste(" Output prefix:", ARGS[3]))
message("")


###################################################
### CODE BELOW THIS POINT IS NONFUNCTIONAL ###
###################################################

###################################################
### READ IN Q-FILE ###
###################################################
message("Reading in Q-file:")
message("")

contiglen <- readr::read_table(file = ARGS[1], col_names = c("","length"))
prefix <- ARGS[2]
prefix_plotname <- str_replace(prefix, "_", " ")
prefix_filename <- str_replace(prefix, " ", "_")
big_contig_threshold = 1000000

contiglen <- contiglen %>% 
  select(length) %>%
  arrange(desc(length)) %>%
  mutate(name = str_c(row_number())) %>%
  mutate(name = factor(name, levels = name)) %>%
  mutate(big_contig = 
           ifelse(length >= big_contig_threshold, yes = TRUE, no = FALSE))

###################################################
### READ IN FAM FILE ###
###################################################

###################################################
### PLOT CHROMOSOMER OUTPUT ###
###################################################
message("Plotting contig lengths:")
message("")

# Define ggplot2 theme:
assembly_theme <- theme(
  title = element_text(size = 18, face = "bold"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  legend.position = "none"
  #scale_x_continuous(n.breaks = 10)
  # scale_y_continuous()
)

# Plot the item
plot1 <- ggplot(contiglen, aes(x = name, y = length, fill = big_contig)) +
  geom_col(color = "black") +
  geom_hline(yintercept = big_contig_threshold, color = "red", linewidth = 1) +
  labs(title = str_c(prefix_plotname," contig lengths in descending order")) +
  xlab("Contig name") +
  ylab("Contig lenth (bp)") +
  theme(legend.position = "none") +
  theme_classic() +
  assembly_theme
plot1
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
