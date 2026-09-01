#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("Biostrings")
}

library(dada2)
library(Biostrings)

path <- setwd("/home/jmangino/metabarcoding/data")

list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.metadata <- read.csv("sampleData_alimentos.csv")

library(ggplot2)

p <- plotQualityProfile(fnFs[1:2])
ggsave("perfil_calidad.png", plot = p, width = 7, height = 5, dpi = 300)

p2 <- plotQualityProfile(fnRs[1:2])
ggsave("perfil_calidad_R.png", plot = p2, width = 7, height = 5, dpi = 300, device = "png")

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(290, 200),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

save(out, file = "out")

load("errF")

load("errR")

png("perfil_errorF.png", width = 2100, height = 1500, res = 300)
# 2. Dibuja el gráfico dentro del archivo
plotErrors(errF, nominalQ = TRUE)
# 3. Cierra y guarda el archivo de forma segura
dev.off()

png("perfil_errorR.png", width = 2100, height = 1500, res = 300)

plotErrors(errR, nominalQ = TRUE)

dev.off()

# --------------------------------------------------------------

derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names
save(derepFs, file = "derepFs")

derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.names
save(derepRs, file = "derepRs")

# --------------------------------------------------------------

dadaFs <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
save(dadaFs, file = "dadaFs")

dadaRs <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
save(dadaRs, file = "dadaRs")

dadaFs[[1]]
dadaRs[[1]]

mean(dadaFs[[1]]$quality)
mean(dadaRs[[1]]$quality)

# --------------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
save(mergers, file = "mergers")

head(mergers[[1]])

# --------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)
save(seqtab, file = "seqtab")

dim(seqtab)

table(nchar(getSequences(seqtab)))

seq_dis <- table(nchar(getSequences(seqtab)))

png("barplot.png", width = 2100, height = 1500, res = 300)

barplot(
  seq_dis,
  col = "lightblue",
  border = "darkblue",
  main = "Distribución de longitudes de ASVs",
  xlab = "Longitud de secuencia (bp)",
  ylab = "Número de ASVs"
)

dev.off()


seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% seq(440, 466)]
save(seqtab2, file = "seqtab2")

table(nchar(getSequences(seqtab2)))

seq_dis2 <- table(nchar(getSequences(seqtab2)))

png("barplot2.png", width = 2100, height = 1500, res = 300)

barplot(
  seq_dis2,
  col = "lightblue",
  border = "darkblue",
  main = "Distribución de longitudes de ASVs filtradas",
  xlab = "Longitud de secuencia (bp)",
  ylab = "Número de ASVs"
)
dev.off()

seqtab.nochim <- removeBimeraDenovo(
  seqtab2,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
save(seqtab.nochim, file = "seqtab.nochim")
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab2)

# --------------------------------------------------------------

dna <- DNAStringSet(getSequences(seqtab.nochim))
meansequencelength <- as.data.frame(table(nchar(getSequences(dna))))
meansequencelength$Var1 <- as.numeric(as.character(meansequencelength$Var1))
meansequencelength$total <- meansequencelength$Var1 * meansequencelength$Freq
sum(meansequencelength$total) / sum(meansequencelength$Freq)

getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# --------------------------------------------------------------

taxa <- assignTaxonomy(
  seqtab.nochim,
  "/home/jmangino/metabarcoding/data/silva_nr_v132_train_set.fa.gz",
  multithread = TRUE
)

taxa <- addSpecies(
  taxa,
  "/home/jmangino/metabarcoding/data/silva_species_assignment_v132.fa.gz"
)

save(taxa, file = "taxa")