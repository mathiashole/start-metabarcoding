#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------

if (!requireNamespace("mia", quietly = TRUE)) {
  BiocManager::install("mia")
}

library(mia)
library(vegan)
library(ggplot2)

setwd("/home/jmangino/metabarcoding/data")

# Cargar el objeto tse de la carpeta de drive Día 4.
tse <- readRDS("tse_dia4_base.rds")

# Ver qué contiene el objeto.
tse

dim(assay(tse, "counts"))  
# Cantidad de ASVs (filas) y cantidad de muestras (columnas)

rownames(assay(tse, "counts"))[1:6]  
# Muestra los nombres de las primeras seis ASVs

colnames(assay(tse, "counts"))[1:6]  
# Muestra los identificadores de las primeras seis muestras

assayNames(tse)  # Muestra los nombres de las distintas tablas almacenadas en el objeto tse

assay(tse, "counts")[1:6, 1:6]  
# Conteos originales: número de lecturas de cada ASV en cada muestra

assay(tse, "relabundance")[1:6, 1:6]  
# Abundancias relativas: proporción que representa cada ASV dentro de cada muestra

assay(tse, "subsampled")[1:6, 1:6]  
# Conteos luego de la rarefacción: todas las muestras fueron igualadas a la misma profundidad

assay(tse, "rarefiedRelabundance")[1:6, 1:6]  
# Abundancias relativas calculadas a partir de los conteos rarefactados

assay(tse, "pa")[1:6, 1:6]  
# Presencia/ausencia: 1 indica que la ASV está presente y 0 que está ausente

# Un pequeño chequeo

colSums(assay(tse, "counts"))[1:6]  
# Cantidad total de lecturas originales en cada una de las primeras seis muestras

colSums(assay(tse, "relabundance"))[1:6]  
# Cada columna debería sumar 1 porque son proporciones

colSums(assay(tse, "subsampled"))[1:6]  
# Todas deberían sumar 6513: profundidad utilizada para rarefaccionar

colSums(assay(tse, "rarefiedRelabundance"))[1:6]  
# Cada muestra debería sumar 1 porque son abundancias relativas

table(assay(tse, "pa"))  

