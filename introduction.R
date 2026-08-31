#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------

R.version.string
getwd()

# Todo lo que descarguemos quedará dentro de la carpeta "data".
#dir.create("/home/jmangino/metabarcoding/data", showWarnings = FALSE)

setwd("/home/jmangino/metabarcoding/data")

if (!requireNamespace("googledrive", quietly = TRUE)) {
  install.packages("googledrive")
}

library(googledrive)

drive_deauth() # No hace falta iniciar sesión porque la carpeta es pública.

carpeta_taller <- as_id(
  "https://drive.google.com/drive/folders/1SOhCrJXgrtucWavxsxizLY4FErL3Almk"
)

archivos_principales <- drive_ls(carpeta_taller)
archivos_principales

metadata_local <- file.path("/home/jmangino/metabarcoding/data", "sampleData_alimentos.csv")

# Descargamos la metadata solamente si todavía no está en la computadora.
if (!file.exists(metadata_local)) {
  metadata_drive <- archivos_principales[
    archivos_principales$name == "sampleData_alimentos.csv",
  ]

  if (nrow(metadata_drive) == 0) {
    stop("No se encontró sampleData_alimentos.csv en Google Drive.")
  }

  drive_download(
    metadata_drive,
    path = metadata_local,
    overwrite = FALSE
  )
}

# Leemos la tabla.
metadata <- read.csv(
  metadata_local,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
