#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# Evanno_admixture.R
# KCollier 22 Jan 2024

# This is a command-line-usable script that takes a list of log likelihood values from admixture.
# It summarizes them and finds the best K value with the Evanno et al. method.

###################################################
### PARSE ARGUMENTS ###
###################################################
# Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  stop("Usage: ./Evanno_admixture.R <input_loglikelihoods.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  message("No prefix was provided. 'output' used for all outfiles.")
  ARGS[2] <- "output"
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
  library(patchwork)

  message("All necessary packages loaded.")
  message("")
  
  message("Parameters interpreted as:")
  message(paste(" Input loglikelihoods.txt:", ARGS[1]))
  message(paste(" Output prefix:", ARGS[2]))
  message("")
  
  table_loc <- ARGS[1]
  output_prefix <- ARGS[2]
  
  ###################################################
  ### HELPER FUNCTIONS ###
  ###################################################
  get_1st_deriv <- function(index) {
    mean_loglikes$mean_ln[index+1] - mean_loglikes$mean_ln[index]
  }
  
  ###################################################
  ### READ IN LOGLIKELIHOOD FILE ###
  ###################################################
  message("Reading in log likelihood file:")
  message("")
  
  table_loc <- "loglike.K1-K9.txt"
  #q_cols <- c("K1","K2")
  loglikes <- readr::read_table(file = table_loc) %>%
    mutate(ln = log(abs(log_likelihood)) * -1)

# Get summary table
  mean_loglikes <- loglikes %>%
    dplyr::group_by(num_K) %>%
    dplyr::summarize(mean_ln = mean(ln), std_dev = sd(ln))
 
# Get table of derivatives
  num_K = nrow(mean_loglikes) # used to set the number of loops to run
  # Create the 'derivs' tibble. We have to create this one separately, because applying the functions across rows of K values results in NA for K=1
  derivs <- tibble(ln_1st_deriv = NA, absolute_val_ln_2nd_deriv = NA, delta_K = NA) # We manually assign 'NAs' for K=1
  
  # Calculate and add rows for K2-K(Max_K-1)
    # The loop relies on adding values to i in order to add/subtract derivative values. See get_1st_deriv() for details.
  for (i in 1:(num_K-1)) {
    new_1st_deriv <- get_1st_deriv(i)
    new_2nd_deriv <- abs(get_1st_deriv(i + 1) - new_1st_deriv)
    new_delta_K <- (new_2nd_deriv / mean_loglikes$std_dev[i+1])
    
    derivs <- add_row(.data = derivs,
            ln_1st_deriv = new_1st_deriv,
            absolute_val_ln_2nd_deriv = new_2nd_deriv,
            delta_K = new_delta_K)
  }

# The first deriv. for K='num_K' is technically calculated, but is an invalid result.
# This is because there's no larger mean_loglike to subtract K='num_K''s mean_loglike from. Therefore, we set it back to 'NA'.
  derivs$ln_1st_deriv[num_K] <- NA 
  
##########################################################
### CREATE UNIFIED DATA FRAME AND SAVE TO LOCAL DISK ###
##########################################################
mean_loglikes <- bind_cols(mean_loglikes, derivs)  # binds columns and writes back to "mean_loglikes"
write_csv(mean_loglikes, file = stringr::str_c(output_prefix,"_maxK", num_K,".csv"))

  
##############################################
### PLOT THE EVANNO OUTPUT ###
##############################################
### OVERLAID CHART - BAD IDEA
#lines_overlaid <- ggplot(data = mean_loglikes) + 
#    geom_line(aes(x = num_K, y = mean_ln), color = 'red') +
#    geom_point(aes(x = num_K, y = mean_ln), color = 'red') + 
#    geom_line(aes(x = num_K, y = delta_K), color = 'blue') +
#    geom_point(aes(x = num_K, y = delta_K), color = 'blue') +
#    theme_classic()
  
mean_ln_only <- ggplot(data = mean_loglikes, aes(x = num_K, y = mean_ln)) + 
  xlab(label = "Value of K" ) +
  ylab(label = "Mean LnP(K)") +
  scale_x_continuous(breaks = seq.int(from = 0, to = 9, by = 1)) +
  scale_y_continuous(breaks = seq(from = 0, to = min(mean_loglikes$mean_ln), by = -0.01)) +
  geom_line(color = 'red') +
  geom_point(color = 'red') + 
  geom_errorbar(aes(ymin = (mean_ln - std_dev), 
                    ymax = (mean_ln + std_dev))) +
  theme_classic()
mean_ln_only
  
delta_K_only <- ggplot(data = mean_loglikes, aes(x = num_K, y = delta_K)) +
  #title(main = "") +
  xlab(label = "Value of K") +
  ylab(label = "Delta K") +
  scale_x_continuous(breaks = seq.int(from = 0, to = 9, by = 1)) +
  geom_line(color = 'blue') +
  geom_point(color = 'blue') +
  theme_classic()

# Use patchwork to display them side-by-side
mean_ln_deltaK_plot <- mean_ln_only + delta_K_only

###################################################
### SAVE PLOT TO LOCAL FILESYSTEM ###
###################################################
ggsave(
  stringr::str_c(output_prefix,"_maxK", num_K, ".pdf"),
  plot = mean_ln_deltaK_plot,
  device = "pdf",
  width = 500,
  height = 250,
  units = "mm",
  dpi = 300,
)
message("Program end.")