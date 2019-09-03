
# Packages ----------------------------------------------------------------

library(pacman)
p_load(foreach, data.table, brainGraph, doParallel, tidyverse, abind)

# Number of cores for parallel
registerDoParallel(cores = (detectCores(logical = F)/2))

# Functions ---------------------------------------------------------------

#source("")

# Day and Directories -----------------------------------------------------

today <- format(Sys.Date(), "%Y%m%d")

savedir <- "outData/graphs/"
savedirDay <- paste0(savedir, today, "_")
