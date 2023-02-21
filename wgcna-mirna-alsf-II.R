# Script to perform weighted correlation gene co-expression network analysis (WGCNA)
# Author: Shehbeel Arif
# Code adapted from ALSF (https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#4_Identifying_co-expression_gene_modules_with_WGCNA_-_RNA-seq)
# Children's Hospital of Philadelphia

# Load libraries
library(tidyverse)
#library(ggplot2)
library(DESeq2)
library(WGCNA)
library(ggpubr)
allowWGCNAThreads() # allow multi-threading (optional)

################################################################################
# SETTING UP DIRECTORIES & LOADING DATA
# Set up analysis directories
# Create the data folder if it doesn't exist
# if (!dir.exists("data")) {
#   dir.create("data")
# }

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}




## Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, "scratch", "wgcna-mirna", "results")
plots_dir <- file.path(root_dir, "scratch", "wgcna-mirna", "plots")

## output files for final plots
# wgcna_pdf <- file.path(output_dir, "wgcna.pdf")
# wgcna_tiff <- file.path(output_dir, "wgcna.tiff")

# Load datasets
mirna_clusters <- paste0(data_dir, "/", "pbta_clusters_mirna.csv")
mirna_counts <- paste0(data_dir, "/", "pbta_raw_counts_mirna.rds")
clusters <- read.csv(mirna_clusters)
counts <- readRDS(mirna_counts)

################################################################################
## DATA PREPROCESSING

## Fix column names
# Remove first four characters in sample_ids
names(counts) <- substring(names(counts),5)
# Remove last two characters in sample_ids
names(counts) <- substr(names(counts), 1, nchar(names(counts))-2)

# Rearrange counts column names by cluster data
counts <- counts[,clusters$Sample_ID]

# Convert the assigned cluster values into factors
clusters$Cluster <- factor(clusters$Cluster)

################################################################################
# QUALITY CONTROL - OUTLIER DETECTION
# Detect outlier genes
gsg <- goodSamplesGenes(t(counts))
summary(gsg)
gsg$allOK

table(gsg$goodGenes) # No outlier genes detected
table(gsg$goodSamples) # No outlier samples detected

# Remove genes that are detectd as outliers
#counts <- counts[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering - method 1
htree <- hclust(dist(t(counts)), method = "average")
plot(htree)
# Sample outliers: 7316-496 (Ependymoma), 7316-230 (Ganglioglioma), 7316-88 (Ependymoma)

# Detect outlier samples using PCA - method 2

pca <- prcomp(t(counts))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
# Sample outliers: 7316-496 (Ependymoma)

# Exclude outlier samples
samples.to.be.excluded <- c("7316-496")
counts.subset <- counts[, !colnames(counts) %in% samples.to.be.excluded]

# Remove the outlier sample from colData as well
clusters <- clusters[clusters$Sample_ID != "7316-496",]

################################################################################
# NORMALIZATION

# Normalize counts with DESeq
dds <- DESeqDataSetFromMatrix(countData = counts.subset,
                              colData = clusters,
                              design = ~1) # not specifying the model

# Remove all genes with counts <15 in more than 75% of samples (0.75 * 237 samples = 177.5)
## Suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) > 178]
nrow(dds75) # 602 genes left

# Perform variance stabilization
#dds_norm <- vst(dds75) # Did not work b/c it is less than 'nsub' rows
dds_norm <- varianceStabilizingTransformation(dds75)

# Transpose the counts
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

### Plot the normalized expressions
norm_counts_df <- data.frame(t(normalized_counts)) %>%
  mutate(
    Gene_id = row.names(t(normalized_counts))
  ) %>%
  pivot_longer(-Gene_id)

norm_counts_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "VST Normalized miRNA Expression",
    x = "Sample ID",
    y = "VST Normalized Expression"
  )

################################################################################
### WGCNA ANALYSIS ###
################################################################################
# Determine soft-threshold for scale-free network

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(normalized_counts,
                         #blockSize = 30,
                         powerVector = powers,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed", # genes are similar if they are strongly correlated; "signed" applied only to genes that are positively correlated
                         verbose = 5)

## Make plots of the soft-thresholds
par(mfrow = c(1,2));
cex1 = 0.8;

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
abline(h = 0.80, col = "red")
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


################################################################################
# Run WGCNA
chosen_power <- 10
temp_cor <- cor       
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = chosen_power, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234 # there's some randomness associated with this calculation
                          # so we should set a seed
)


# Return cor function to original namespace
cor <- temp_cor

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

# Look at number of genes per Modules
table(bwnet$colors)
# 0   1   2   3 
# 330 185  55  32 

################################################################################
# WGCNA MODULE DENDROGRAM

# Convert labels to colors for plotting
mergedColors <- labels2colors(bwnet$colors)
unmergedColors <- labels2colors(bwnet$unmergedColors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(mergedColors[bwnet$blockGenes[[1]]], unmergedColors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05 )


################################################################################
## WHICH MODULES HAVE BIGGEST DIFFERENCES ACROSS CLUSTER GROUPS?

# Check if SAMPLE_ID orders are the same as module_eigengenes
all.equal(clusters$Sample_ID, rownames(module_eigengenes))

# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ clusters$Cluster)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

# Save the Limma stats data to a CSV file
readr::write_csv(stats_df,
                 file = file.path("results", "wgcna-mirna-cluster-limma.csv")
)

################################################################################
# CORRELATION BETWEEN MODULES AND CLUSTERS

# As a sanity check, let’s use ggplot to see what module 3’s eigengene looks like between treatment groups.
module_df <- module_eigengenes %>%
  tibble::rownames_to_column("Sample_ID") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(clusters %>%
                      dplyr::select(Sample_ID, Cluster),
                    by = c("Sample_ID" = "Sample_ID")
  )

# Save the module eigengenes Sample_ID dataframe
readr::write_csv(stats_df,
                 file = file.path("results", "wgcna-mirna-sample_id-eigengenes.csv")
)

# Reorder modules so similar modules are next to each other
module_order <- names(module_df)

# Add clusters
temp_module_df <- module_df %>% select("ME2":"Cluster")

# tidy & plot data
temp_module_df <- temp_module_df %>%
  pivot_longer(-Cluster) %>%
  mutate(
    module = factor(name, levels = module_order)
  )

temp_module_df %>% ggplot(., aes(x=Cluster, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "#f7f7f7",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Cluster Relationships", y = "Modules", x="Cluster", fill="corr")


################################################################################
# Make boxplot of Module Eigengene expression by clusters

# Create list of pairwise vectors to statistically compare clusters
my_comparisons <- list(c("1","2"),
                       c("3","2"),
                       c("4","2")
                       )

# Make plot
ggplot(
  module_df,
  aes(
    x = Cluster,
    y = ME3,
    color = Cluster
  )
) +
  labs(title = "ME3 Expression in Clusters") +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  # Add global p-value
  stat_compare_means(label.y = 0.35) + # ME0=0.15, ME1=0.25, ME2=0.15, ME3=0.25
  # Add pairwise comparisons p-values
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') +
  theme_classic()

################################################################################
## What genes are part of Module 3?
# Add Modules to genes
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
# Filter for genes part of only ME3
gene_module_key %>%
  dplyr::filter(module == "ME3")

# Save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path("results", "PBTA_miRNA_wgcna_gene_to_module.tsv")
)


################################################################################
# MODULE-SPECIFIC HEATMAPS

# Make custom heatmap function
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = clusters,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its Sample_ID
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("Sample_ID")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(Sample_ID, Cluster) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "Sample_ID") %>%
    # Arrange by patient and time point
    dplyr::arrange(Cluster) %>%
    # Store sample
    tibble::column_to_rownames("Sample_ID")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    cluster = col_annot_df$Cluster,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each clusters in Cluster
    col = list(cluster = c("1" = "#f1a340", "2" = "#998ec3", "3" = "#0E86D4", "4" = "#E7625F"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = TRUE,
                                     show_column_names = FALSE,
                                     column_title = "ME3 Gene Expression"
  )
  
  # Return heatmap
  return(heatmap)
}

# Make Module 3 heatmap
mod_3_heatmap <- make_module_heatmap(module_name = "ME3")

# Print out the plot
mod_3_heatmap


# Save this plot to PNG
png(file.path("results", "PBTA_mirna_module_3_heatmap.png"))
mod_3_heatmap
#dev.off()

## For comparison look at other modules
#mod_2_heatmap <- make_module_heatmap(module_name = "ME2")
# Print out the plot
#mod_2_heatmap


################################################################################
## FINDING DRIVER GENES
chooseTopHubInEachModule(normalized_counts,
                         gene_module_key$module, 
                         power = chosen_power, 
                         type = "signed")

# ME0           ME1           ME2           ME3 
# "miR-32-5p" "miR-4695-5p"  "miR-154-5p"   "miR-17-5p" 


################################################################################
## EXAMINE EXPRESSION PROFILES

# pick out a few modules of interest here
modules_of_interest <- c("ME0", "ME1", "ME2", "ME3")

# Pull out list of genes in that module
submod <- gene_module_key %>%
  subset(module %in% modules_of_interest)

row.names(gene_module_key) <- gene_module_key$gene

# Get normalized expression for those genes
temp_normalized_counts <- t(normalized_counts)
subexpr <- temp_normalized_counts[submod$gene,]

submod_df <- data.frame(subexpr) %>%
  mutate(
    gene = row.names(.)
  ) %>%
  pivot_longer(-gene) %>%
  mutate(
    module = gene_module_key[gene,]$module
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Sample_ID",
       y = "normalized expression") +
  scale_color_manual(values=c("grey", "turquoise", "blue", "brown"))

################################################################################
## GENERATE NETWORKS
# The network file can be generated for Cytoscape or as an edge/vertices file.

genes_of_interest <- module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest <- norm_counts[genes_of_interest$gene_id,]


# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM <- TOMsimilarityFromExpr(t(expr_of_interest),
                             power = picked_power)

# Add gene names to row and columns
row.names(TOM) <- row.names(expr_of_interest)
colnames(TOM) <- row.names(expr_of_interest)

edge_list <- data.frame(TOM) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  mutate(gene2 = gsub("*\\.","-",as.character(gene2))) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )


# Select rows if gene1 and gen2 are from same module
edge_list[edge_list$module1 == edge_list$module2,]
# Select rows if gene1 and gen2 have a correlation >0.1
new_edge_list <- edge_list[edge_list$correlation > 0.5,]

a <- new_edge_list$gene1
b <- new_edge_list$gene2
idx <- order(c(seq_along(a), seq_along(b)))
edges <- c(a,b)[idx]

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")

write_delim(as.data.frame(TOM),
            file = "TOM_edgelist.tsv",
            delim = "\t")

read.csv("edgelist.tsv", sep="\t")

# Creating a Network graph using igraph
# Tutorial: https://kateto.net/networks-r-igraph
library(igraph) # Load the igraph package
rm(list = ls()) # Remove all the objects we created so far.

g1 <- graph( edges=edges, directed=F) 

plot(g1)


