# ALSF guide
# Script to perform weighted correlation gene co-expression network analysis (WGCNA)
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

# Load libraries
library(tidyverse)
#library(gridExtra)
#library(ggplot2)
library(DESeq2)
library(WGCNA)
allowWGCNAThreads() # allow multi-threading (optional)

## Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "scratch", "output")

## output files for final plots
wgcna_pdf <- file.path(output_dir, "wgcna.pdf")
wgcna_tiff <- file.path(output_dir, "wgcna.tiff")

# Load datasets
mirna_clusters <- paste0(data_dir, "/", "pbta_clusters_mirna.csv")
mirna_counts <- paste0(data_dir, "/", "pbta_raw_counts_mirna.rds")
clusters <- read.csv(mirna_clusters)
counts <- readRDS(mirna_counts)


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

# Remove all genes with counts <15 in more than 75% of samples (0.75 * 238 samples = 178.5)
## Suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) > 179]
nrow(dds75) # 598 genes left

# Perform variance stabilization
#dds_norm <- vst(dds75) # Did not work b/c it is less than 'nsub' rows
dds_norm <- varianceStabilizingTransformation(dds75)
#dds_norm <- varianceStabilizingTransformation(dds)


## WGCNA ANALYSIS ##
# Transpose the counts
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

## NETWORK CONSTRUCTION
# Determine parameters for WGCNA
sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

# Run WGCNA
bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 9, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
)

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)


## Which modules have biggest differences across treatment groups?
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

# As a sanity check, let’s use ggplot to see what module 18’s eigengene looks like between treatment groups.
module_df <- module_eigengenes %>%
  tibble::rownames_to_column("Sample_ID") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(clusters %>%
                      dplyr::select(Sample_ID, Cluster),
                    by = c("Sample_ID" = "Sample_ID")
  )

################################################################################
# Reorder modules so similar modules are next to each other

module_order <- names(module_df)

# Add clusters
temp_module_df <- module_df %>% select("ME2":"Cluster")

# tidy & plot data
mME <- temp_module_df %>%
  pivot_longer(-Cluster) %>%
  mutate(
    #name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=Cluster, y=name, fill=value)) +
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



# Make boxplot of Module Eigengene expression by clusters
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
  theme_classic()

## What genes are part of Module 3?
# Add Modules to genes
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
# Filter for genes part of only ME3
gene_module_key %>%
  dplyr::filter(module == "ME3")

# save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path("results", "PBTA_miRNA_wgcna_gene_to_module.tsv")
)



################################################################################
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
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
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
dev.off()

## For comparison look at other
mod_0_heatmap <- make_module_heatmap(module_name = "ME0")
# Print out the plot
mod_0_heatmap







