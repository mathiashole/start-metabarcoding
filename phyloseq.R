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

# --------------------------------------------------------------------------------

data_dir <- "/home/jmangino/metabarcoding/data"

seqtab_file <- file.path(data_dir, "seqtab.nochim.rds")
taxa_file   <- file.path(data_dir, "taxa_gtdb")
meta_file   <- file.path(data_dir, "sampleData_alimentos.csv")

# ---------------------------------------------------------------------------------

seqtab.nochim <- readRDS(seqtab_file)

class(seqtab.nochim)
dim(seqtab.nochim)

head(rownames(seqtab.nochim))
substr(colnames(seqtab.nochim)[1:3], 1, 80)

# --------------------------------------------------------------------------------

metadata <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(metadata)
dim(metadata)

colnames(metadata)

rownames(metadata) <- metadata$Run

# Creamos un nombre sin espacios para usarlo más fácilmente en el práctico.
metadata$Sampling_position <- metadata[["Sampling position"]]

# --------------------------------------------------------------------------------

load(taxa_file)

class(taxa_gtdb)
dim(taxa_gtdb)
head(taxa_gtdb)

colnames(taxa_gtdb)

# --------------------------------------------------------------------------------

samples_in_counts <- rownames(seqtab.nochim)
samples_in_meta   <- rownames(metadata)

setdiff(samples_in_counts, samples_in_meta)
setdiff(samples_in_meta, samples_in_counts)

# Reordenar metadata de acuerdo con la tabla de ASVs
metadata_ps <- metadata[samples_in_counts, , drop = FALSE]

identical(
  rownames(seqtab.nochim),
  rownames(metadata_ps)
)

# --------------------------------------------------------------------------------

asvs_in_counts <- colnames(seqtab.nochim)
asvs_in_tax    <- rownames(taxa_gtdb)

length(asvs_in_counts)
length(asvs_in_tax)

common_asvs <- intersect(
  asvs_in_counts,
  asvs_in_tax
)

length(common_asvs)

seqtab_ps <- seqtab.nochim[, common_asvs, drop = FALSE]
taxa_ps   <- taxa_gtdb[common_asvs, , drop = FALSE]

# --------------------------------------------------------------------------------
# renombrar la secuencias de ASVs para que sean más fáciles de manejar

# Guardamos primero la correspondencia entre el nuevo ID y la secuencia original.
asv_sequences <- colnames(seqtab_ps)
asv_ids <- paste0("ASV", seq_along(asv_sequences))

asv_map <- data.frame(
  ASV = asv_ids,
  Sequence = asv_sequences,
  stringsAsFactors = FALSE
)

# Renombramos de la misma forma la tabla de abundancias y la tabla taxonómica.
colnames(seqtab_ps) <- asv_ids
rownames(taxa_ps) <- asv_ids

head(asv_map)
head(colnames(seqtab_ps))
head(rownames(taxa_ps))

# --------------------------------------------------------------------------------
# Create phyloseq object

OTU <- otu_table(
  seqtab_ps,
  taxa_are_rows = FALSE
)

OTU
taxa_are_rows(OTU)

# --------------------------------------------------------------------------------
# create taxa table
TAX <- tax_table(
  as.matrix(taxa_ps)
)

TAX

# --------------------------------------------------------------------------------
# create sample data
SAM <- sample_data(
  metadata_ps
)

SAM

# --------------------------------------------------------------------------------
# create phyloseq object
ps <- phyloseq(
  otu_table = OTU,
  tax_table = TAX,
  sample_data = SAM
)

ps

nsamples(ps)
ntaxa(ps)

# --------------------------------------------------------------------------------
# explore the phyloseq object
head(sample_names(ps))
head(taxa_names(ps))

rank_names(ps)
sample_variables(ps)

otu_table(ps)[1:5, 1:5]
head(as(tax_table(ps), "matrix"))
head(data.frame(sample_data(ps)))

#--------------------------------------------------------------------------------
# depth of sequencing
depth <- sample_sums(ps)

head(depth)
summary(depth)

min(depth)
median(depth)
max(depth)

max(depth) / min(depth)

# plots of depth of sequencing

depth_df <- data.frame(
  Reads = sample_sums(ps),
  sample_data(ps),
  check.names = FALSE
)

png("depth_plot.png", width = 2100, height = 1500, res = 300)

ggplot(
  depth_df,
  aes(
    x = reorder(Sample, Reads),
    y = Reads,
    fill = Sampling_position
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Muestra",
    y = "Número de reads",
    fill = "Posición de muestreo",
    title = "Profundidad de secuenciación"
  ) +
  theme_bw()

dev.off()

# --------------------------------------------------------------------------------
# explore the phyloseq object
head(
  sort(
    taxa_sums(ps),
    decreasing = TRUE
  ),
  20
)

# --------------------------------------------------------------------------------
# Subsetting and filtering
# filtering samples
ps_sticking <- subset_samples(
  ps,
  Sampling_position == "Sticking"
)

ps_sticking

# ---------------------------------------------------------------------------------
# filtering taxa
ps_firmicutes <- subset_taxa(
  ps,
  Phylum == "Firmicutes"
)

ps_firmicutes

# --------------------------------------------------------------------------------
# Keep ASV with at least 100 reads across all samples

ps_abundant <- prune_taxa(
  taxa_sums(ps) >= 100,
  ps
)

ps_abundant

# --------------------------------------------------------------------------------
# Keep ASV with at least 5 sample
ps_prevalent <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= 5,
  prune = TRUE
)

ps_prevalent

# ---------------------------------------------------------------------------------
# Clustering taxonomically
ps_phylum <- tax_glom(
  ps,
  taxrank = "Phylum",
  NArm = TRUE
)

ps_genus <- tax_glom(
  ps,
  taxrank = "Genus",
  NArm = TRUE
)

c(
  ASVs = ntaxa(ps),
  Phyla = ntaxa(ps_phylum),
  Genera = ntaxa(ps_genus)
)

# --------------------------------------------------------------------------------
# Taxonomic composition

ps_phylum_rel <- transform_sample_counts(
  ps_phylum,
  function(x) x / sum(x)
)

plot_bar(
  ps_phylum_rel,
  x = "Sample",
  fill = "Phylum"
) +
  labs(
    x = "Muestra",
    y = "Abundancia relativa",
    title = "Composición taxonómica a nivel de Phylum"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )