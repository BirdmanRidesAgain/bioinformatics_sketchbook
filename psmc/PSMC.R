#!/usr/bin/env Rscript 

### PARSE ARGUMENTS
  # Positional arguments
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  stop("Usage: ./PSMC.R <input.fa> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  ARGS[2] = "output"
}

### CHECKS FOR ALL NECESSARY PACKAGES
# Check if psmcr is installed. Install if it is not present.
if (!requireNamespace("psmcr", quietly = FALSE)) 
{
  install.packages("devtools")
  devtools::install_github("emmanuelparadis/psmcr/psmcr")
}
# Check if ape is installed. Install if it is not present.
if (!requireNamespace("ape", quietly = FALSE))
{
  install.packages("ape")
}
message("All necessary packages loaded.")
message("")

message("Parameters interpreted as:")
message(paste(" Input fasta:", ARGS[1]))
message(paste(" Output prefix:", ARGS[2]))
message("")

### PSMC PROCESSING
  # Import data
INPUT_FA <- ape::read.FASTA(ARGS[1])
OUT_PREFIX <- ARGS[2]

  # Run PSMC simulations
NUM_ITERATIONS <- 30
NUM_BOOTSTRAPS <- 100
NUM_CORES <- parallel::detectCores()

message("Beginning PSMC iterations...")
message(paste("# Iterations:", NUM_ITERATIONS))
message(paste("# Bootstrap replicates:", NUM_BOOTSTRAPS))
message(paste("# Available cores:", NUM_CORES))

message("")
PSMC_OUT <- psmcr::psmc(INPUT_FA, niters = NUM_ITERATIONS, B = NUM_BOOTSTRAPS, mc.cores = NUM_CORES)
saveRDS(PSMC_OUT, sub(" ","",paste(toString(OUT_PREFIX),".rds")))

q()
