# Creating a Network graph using igraph
# Tutorial: https://kateto.net/networks-r-igraph
library(igraph) # Load the igraph package
#rm(list = ls()) # Remove all the objects we created so far.

# Load in edge_list
edge_list <- read.csv("edgelist.tsv", sep="\t")

# Select rows if gene1 and gen2 are from same module
new_edge_list <- edge_list[edge_list$module1 == edge_list$module2,]
# Select rows with blue module
new_edge_list <- new_edge_list[new_edge_list$module1 == "yellow",] # Blue = 0.4, Brown = 0.2, Grey = 0.45, Turquoise = 0.5, Yellow = 0.25
# Select rows if gene1 and gen2 have a correlation >0.1
new_edge_list <- new_edge_list[new_edge_list$correlation > 0.25,]

a <- new_edge_list$gene1
b <- new_edge_list$gene2
idx <- order(c(seq_along(a), seq_along(b)))
edges <- c(a,b)[idx]

color_a <- new_edge_list$module1
color_b <- new_edge_list$module2
color_idx <- order(c(seq_along(color_a), seq_along(color_b)))
color <- c(color_a,color_b)[color_idx]

g1 <- graph( edges=edges, directed=F) 

plot(g1, vertex.color=color)
dev.off()

# Create data
#data <- matrix(sample(0:1, 400, replace=TRUE, prob=c(0.8,0.2)), nrow=20)
# TOM_edge_list <- read.csv("TOM_edgelist.tsv", sep="\t")
#network <- graph_from_adjacency_matrix(TOM_edge_list, mode='undirected', diag=F)
#plot(network)



