# Script to perform weighted correlation gene co-expression network analysis (WGCNA)
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

# Load libraries
library(tidyverse)
library(DESeq2)

## Set directories
root_dir <- "/Users/shehbeel/Documents/OpenPBTA-miRNA-Analysis"
data_dir <- file.path(root_dir, "data")
#wgcna_dir <- file.path(root_dir, "scratch", "wgcna-mirna")
output_dir <- file.path(root_dir, "scratch", "output")

## output files for final plots
wgcna_pdf <- file.path(output_dir, "wgcna.pdf")
wgcna_tiff <- file.path(output_dir, "wgcna.tiff")

# Load datasets
mirna_clusters <- paste0(data_dir, "/", "pbta_clusters_mirna.csv")
#mirna_meta <- paste0(data_dir, "/", "pbta_meta_mirna.rds")
mirna_counts <- paste0(data_dir, "/", "pbta_raw_counts_mirna.rds")
clusters <- read.csv(mirna_clusters)
#meta <- readRDS(mirna_meta)
counts <- readRDS(mirna_counts)


## DATA PREPROCESSING

## Fix column names
# Remove first four characters in sample_ids
names(counts) <- substring(names(counts),5)
# Remove last two characters in sample_ids
names(counts) <- substr(names(counts), 1, nchar(names(counts))-2)

# Rearrange counts column names by cluster data
counts <- counts[,clusters$sample_id]

# Convert the assigned cluster values into factors
clusters$cluster <- factor(clusters$cluster)

# Normalize counts with DESeq
dds <- DESeqDataSetFromMatrix(counts,
                              clusters,
                              design = ~cluster)
# Remove lowly expressed genes
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# Main DESeq
ddsDE <- DESeq(dds)

# Normalized counts
#expr_normalized <- counts(ddsDE, normalized=TRUE)

# ## VST normalization
vsd <- varianceStabilizingTransformation(ddsDE)
# library(genefilter)
wpn_vsd <- getVarianceStabilizedData(ddsDE)
# rv_wpn <- rowVars(wpn_vsd)
# summary(rv_wpn)
# q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
# q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
# expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
expr_normalized <- wpn_vsd

# Plot the normalized expressions
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )


## WGCNA ANALYSIS ##
# Transpose normalized counts matrix
input_mat <- t(expr_normalized)

# BiocManager::install("WGCNA")
library(WGCNA)
allowWGCNAThreads() # allow multi-threading (optional)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

## Make plots of the soft-thresholds
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

## Pick a soft threshold power based on the plots, preferably near the curve of the plot


picked_power = 6
temp_cor <- cor       
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

# Create network using the blockwiseModules command
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,  # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

# Return cor function to original namespace
cor <- temp_cor

## Look at the modules
# Convert labels to colors for plotting
mergedColors <- labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)

## Relate module (cluster) assignments to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Export list of genes and their modules
write_delim(module_df,
            file = "gene_modules_no_vst.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME","", .)

# Add treatment names
#MEs0$treatment <- row.names(MEs0)
MEs0$sample_id <- row.names(MEs0)

# Add clusters
MEs0 <- merge(MEs0, clusters, by="sample_id")
MEs0 <- MEs0 %>% select("MEbrown":"cluster")

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-cluster) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=cluster, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Cluster Relationships", y = "Modules", x="Cluster", fill="corr")


## Examine Expression Profiles
# pick out a few modules of interest here
modules_of_interest = c("grey", "turquoise", "blue", "yellow", "brown")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Sample_ID",
       y = "Normalized Expression")


## Generate and Export Networks









