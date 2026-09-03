#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------

library(Biostrings)
library(DECIPHER)

asv_map <- readRDS("asv_map.rds")

# --- Hacer un objeto DNAStringSet nombrado por ASV no por secuencia ---
seqs <- DNAStringSet(asv_map$Sequence)
names(seqs) <- asv_map$ASV

# --- Alineamiento ---
alignment <- AlignSeqs(seqs, anchor = NA, verbose = FALSE)

# Exportar el alineamiento para FastTree
writeXStringSet(alignment, filepath = "alignment_ps.fasta")


## Bash #########################################################################################
# --- Hacer árbol de ML con FastTree usando el alineamiento (etiquetas serán ASV1, ASV2, ...) ---
FastTree -gtr -nt alignment_ps.fasta > ml-tree_ps.nwk
## Bash #########################################################################################


## R
# --- Cargar el árbol ---
tree <- ape::read.tree("ml-tree_ps.nwk")
str(tree)
plot(tree)
treer <- phangorn::midpoint(tree)
ape::write.tree(treer, file = "ml-tree_ps_midpoint.nwk")

# Confirmar que las etiquetas de los tips coinciden con ps
setequal(tree$tip.label, taxa_names(ps))   # TRUE
sum(tree$tip.label %in% taxa_names(ps))    # = ntaxa(ps)

# --- Incluir el árbol y las secuenicas al objeto phyloseq ---
ps_final <- merge_phyloseq(ps, phy_tree(tree), refseq(seqs))

ps_final
refseq(ps_final)

saveRDS(ps_final, "ps_final.rds")


##################################################################
##################################################################
library(mia)
library(miaViz)
library(phyloseq)

# -- Convertir phyloseq a TreeSummarizedExperiment
tse <- convertFromPhyloseq(ps_final)
tse


# -- Plot del árbol con la anotación (taxonomia)
pdf("arbol-tse.pdf")
plotRowTree(
  tse,
  tip_colour_by = "Phylum",   # or Genus, etc. - any column in rowData(tse)
  tip_size_by   = "Genus"     # optional
) +
  theme_bw()
dev.off()


pdf("arbol-tse_r.pdf")
plotRowTree(
  tse_r,
  tip_colour_by = "Phylum",   # or Genus, etc. - any column in rowData(tse)
  tip_size_by   = "Genus"     # optional
) +
  theme_bw()
dev.off()


# -- Relacionar el árbol con los metadatos
# aggregate/transform abundances, then combine tree + heatmap of samples.

# Option A: plot tree with a heatmap of abundances per sample using miaViz
pdf("arbol-tse_muestras.pdf")
plotRowTree(
  tser,
  tip_colour_by = "Phylum"
) 

# Combine with an abundance heatmap aligned to tree tips
plotRowTree(
  tser,
  tip_colour_by = "Phylum",
  add_legend = TRUE
)
dev.off()

# -- Arbol + heatmap de las abundancias
top_taxa <- getTop(tse, top = 30, method = "mean")
tse_top <- tse[top_taxa, ]

pdf("arbol-tse_heatmap")
plotRowTree(
  tse_top,
  tip_colour_by = "Phylum"
) +
  theme_bw() +
  labs(title = "Phylogenetic tree - top 30 ASVs")
dev.off()

# -- If you want sample-level grouping/metadata reflected (e.g. tree + boxplots by group)
# This is usually done by combining a tree plot with a separate ggplot of sample-level data,
# aligned using patchwork, since plotRowTree focuses on row (taxa) annotations, not columns (samples).

library(patchwork)

tree_plot <- plotRowTree(tse_top, tip_colour_by = "Phylum") + theme_bw()

# Ejemplo: abundancia por Sampling_position para los top taxa
abund_df <- as.data.frame(t(assay(tse_top, "counts")))
abund_df$Sample <- rownames(abund_df)
abund_long <- reshape2::melt(abund_df, id.vars = "Sample", variable.name = "ASV", value.name = "Count")
abund_long <- merge(abund_long, as.data.frame(colData(tse_top)), by.x = "Sample", by.y = "row.names")

abund_plot <- ggplot(abund_long, aes(x = ASV, y = Count, fill = Sampling_position)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw()

pdf("arbol-tse_abund.pdf")
tree_plot + abund_plot
dev.off()

##################################################################
library(mia)
library(ComplexHeatmap)
library(circlize)
library(ape)
library(phangorn)

# -- Subset a un numero de taxa manejable
top_taxa <- getTop(tse, top = 100, method = "mean")
tse_top  <- tse[top_taxa, ]

# -- Recortar el árbol para que coincida con el subset previo
tree_full <- rowTree(tse_top)
tree <- keep.tip(tree_full, rownames(tse_top))
length(tree$tip.label) == nrow(tse_top)   # should be TRUE

# -- Poner una raíz al árbol si no la tiene
if (!is.rooted(tree)) {
  tree <- phangorn::midpoint(tree)
}

# -- Constuir una matriz de abundancia, ordenada para ajustarse a tree$tip.label
mat <- as.matrix(assay(tse_top, "counts"))
mat <- mat[tree$tip.label, ]

mat_rel <- prop.table(mat, margin = 2) * 100   # relative abundance %

# -- Convertir el objeto phylo para Heatmap()
tree_ultra <- chronos(tree)
hc   <- as.hclust(tree_ultra)
dend <- as.dendrogram(hc)

# -- Metadatos
sample_meta <- as.data.frame(colData(tse_top))
sample_meta <- sample_meta[colnames(mat_rel), , drop = FALSE]

col <- HeatmapAnnotation(
  Sampling_position = sample_meta$Sampling_position,
  col = list(
    Sampling_position = c(
      "Sticking"  = "#1b9e77",
      "Ripening"  = "#d95f02",
      "Packaging" = "#7570b3"
      # adjust to match levels(factor(sample_meta$Sampling_position))
    )
  ),
  annotation_name_side = "left"
)

# -- Heatmap + Arbol + Metadatos
ht <- Heatmap(
  mat_rel,
  name = "Rel. Abundance (%)",
  col = colorRamp2(c(0, max(mat_rel)), c("white", "darkred")),
  cluster_rows = dend,        # draws the phylogenetic tree as the row dendrogram
  cluster_columns = TRUE,
  top_annotation = col,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4)
)


##################################################################
##################################################################
# Hacer el árbol en R

# R
library(dada2)
library(phangorn)
library(DECIPHER)

mat <- readRDS(file="seqtab.nochim.rds")
secs <- getSequences(mat)
names(secs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))

