# Comparing Unsupervised Clustering Results using 5% and 10% genes
# Author: Shehbeel Arif
# Center for Data-Driven Discovery in Biomedicine (D3b)
# Children's Hospital of Philadelphia

clust10 <- read.csv('pam_pearson_10_k4/pp10k4_cluster_members.csv', row.names = 1)
clust5 <- read.csv('pam_pearson_5_k4/pp5k4_cluster_members.csv', row.names = 1)

library(mclust)
adjustedRandIndex(clust5$cluster, clust10$cluster)
# 0.8310456 agreement between the two classifications