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

# 2. SELECCIONAR LAS MUESTRAS -------------------------------------------------

# Para este ejercicio trabajaremos solamente con muestras de Sticking y de
# After singeing
tse_d4 <- tse[
  ,
  tse$Sampling.position %in% c("Sticking", "After singeing")
]
# Selecciona muestras de la superficie antes y después del chamuscado

# Eliminar categorías que ya no están presentes después de seleccionar muestras.
tse_d4$Sampling.position <- droplevels(tse_d4$Sampling.position)
tse_d4$Farmer <- droplevels(tse_d4$Farmer)

# Comprobar cuántas muestras quedaron en cada grupo.
table(tse_d4$Sampling.position)

# 3. DIVERSIDAD ALFA ----------------------------------------------------------

# La diversidad alfa describe la diversidad dentro de cada muestra.

# Riqueza observada: número de ASV detectadas.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "observed_richness",
  name = "observed"
)

# Shannon: combina riqueza y distribución de abundancias.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "shannon_diversity",
  name = "shannon"
)

# Uniformidad de Simpson: indica cuán uniformes son las abundancias.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "simpson_evenness",
  name = "simpson_evenness"
)

# Chao1: estima la riqueza considerando las ASV poco frecuentes.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "counts",
  index = "chao1_richness",
  name = "chao1"
)

# Construir una tabla con los resultados.
alpha <- as.data.frame(colData(tse_d4))[
  ,
  c(
    "Sample", "Sampling.position", "observed", "shannon",
    "simpson_evenness", "chao1"
  )
]

alpha

# 4. COMPARAR LA DIVERSIDAD ALFA ---------------------------------------------

# Las muestras de Sticking y After singeing corresponden a los mismos
# cuatro animales. Por lo tanto, las observaciones están pareadas.

# Extraer el número que identifica al animal desde el nombre de la muestra.

alpha$animal <- gsub("\\D", "", alpha$Sample)


# Ordenar las etapas para mostrarlas en su orden temporal.

alpha$Sampling.position <- factor(
  alpha$Sampling.position,
  levels = c("Sticking", "After singeing")
)


# Comprobar que los animales estén en el mismo orden en ambas etapas.

animales_sticking <- alpha$animal[
  alpha$Sampling.position == "Sticking"
]

animales_singeing <- alpha$animal[
  alpha$Sampling.position == "After singeing"
]

pares_animales <- data.frame(
  Sticking = animales_sticking,
  After_singeing = animales_singeing
)

pares_animales