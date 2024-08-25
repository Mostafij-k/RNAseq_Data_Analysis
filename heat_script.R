library(ggplot2) 		#heatmap plotting
library(tidyverse)  		# data manipulation
library(cluster)  




library(readr)
dataset= read.csv('heatmap_fold.csv', header = TRUE, sep = ',')
row.names(dataset)= dataset$X
dataset= dataset[, -1]
dataset= data.matrix(dataset)

pheatmap(dataset)

y= dataset






dataset = log10(dataset+1)

dataset= scale(dataset, scale = TRUE)

y=dataset

mycol= colorpanel(75, 'red')


mycol <-  colorpanel(75, "red", "black", "green")

heatmap.2(
  y,
  Rowv = as.dendrogram(hr),
  Colv = as.dendrogram(hc),
  col = mycol,
  density.info = "none",
  trace = "none",
  dendrogram = "both",
  scale = "row",
  labRow = NULL,
  labCol = NULL,
  margins = c(5, 10),
  RowSideColors = mycolhc
)

library(pheatmap)

pheatmap(y) #default parameter

#parameters to modify. Try changing cluster_cols and clusters_rows to turn on and off hierarchical clustering 

pheatmap(y ,kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete")


