#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------

R.version.string
getwd()

# Todo lo que descarguemos quedará dentro de la carpeta "data".
dir.create("/home/jmangino/metabarcoding/data", showWarnings = FALSE)

setwd("/home/jmangino/metabarcoding/data")

if (!requireNamespace("googledrive", quietly = TRUE)) {
  install.packages("googledrive")
}

library(googledrive)

drive_deauth() # No hace falta iniciar sesión porque la carpeta es pública.

carpeta_taller <- as_id(
  "https://drive.google.com/drive/folders/1SOhCrJXgrtucWavxsxizLY4FErL3Almk"
)
