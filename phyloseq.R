#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}

library(phyloseq)
library(ggplot2)
library(dplyr)

set.seed(123)