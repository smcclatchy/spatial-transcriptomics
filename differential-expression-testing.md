---
title: 'Differential Expression Testing'
teaching: 60
exercises: 20
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is the purpose of differential expression testing in bioinformatics?
- Can Moran's I algorithm independently identify region-specific differential expressions that align with results obtained from expert annotations?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Identify differentially expressed genes across different layers using expert annotations.
- Utilize Moran's I algorithm to find spatially variable genes.
- Explore the correlation between genes identified through expert annotations and those detected by Moran's I algorithm.

::::::::::::::::::::::::::::::::::::::::::::::::



## Introduction to Differential Expression Testing

Differential expression testing is crucial in bioinformatics for identifying genes that show significant differences in expression across different samples or groups. 
This method helps find genes that are upregulated or downregulated in specific contexts, providing insights into biological functions and disease mechanisms.

## Moran's I Statistic

Moran's I is a measure used to assess spatial autocorrelation in data, indicating whether similar values are clustered, dispersed, or random. 
In bioinformatics, it's applied to detect genes whose expression patterns exhibit clear spatial structure, aiding in understanding spatially localized biological processes.

## Data Preparation


``` r
plot_seurat_object <- seurat_object[,!is.na(seurat_object$layer_guess)]
```

``` error
Error in eval(expr, envir, enclos): object 'seurat_object' not found
```

``` r
plot_seurat_object$Layers <- plot_seurat_object$layer_guess
```

``` error
Error in eval(expr, envir, enclos): object 'plot_seurat_object' not found
```

``` r
unique_clusters <- unique(plot_seurat_object$Layers)
```

``` error
Error in eval(expr, envir, enclos): object 'plot_seurat_object' not found
```

``` r
num_clusters <- length(unique_clusters)
```

``` error
Error in eval(expr, envir, enclos): object 'unique_clusters' not found
```

``` r
palette <- carto_pal(num_clusters, "Safe")
```

``` error
Error in carto_pal(num_clusters, "Safe"): could not find function "carto_pal"
```

``` r
names(palette) <- unique_clusters
```

``` error
Error in eval(expr, envir, enclos): object 'unique_clusters' not found
```

## Differential Expression Analysis

### Differential Expression Using Expert's Annotation

Based on the exprerts brain layers annotation, as it is indicated here:


``` r
plot_seurat_object <- seurat_object[,!is.na(seurat_object$layer_guess)]
```

``` error
Error in eval(expr, envir, enclos): object 'seurat_object' not found
```

``` r
plot_seurat_object$Layers <- plot_seurat_object$layer_guess
```

``` error
Error in eval(expr, envir, enclos): object 'plot_seurat_object' not found
```

``` r
unique_clusters <- unique(plot_seurat_object$Layers)
```

``` error
Error in eval(expr, envir, enclos): object 'plot_seurat_object' not found
```

``` r
num_clusters <- length(unique_clusters)
```

``` error
Error in eval(expr, envir, enclos): object 'unique_clusters' not found
```

``` r
palette <- carto_pal(num_clusters, "Safe")
```

``` error
Error in carto_pal(num_clusters, "Safe"): could not find function "carto_pal"
```

``` r
names(palette) <- unique_clusters
```

``` error
Error in eval(expr, envir, enclos): object 'unique_clusters' not found
```

``` r
p <- SpatialDimPlot(plot_seurat_object, group.by = 'Layers', cols=palette) +
  theme(legend.position = "right")
```

``` error
Error in eval(expr, envir, enclos): object 'plot_seurat_object' not found
```

``` r
print(p)
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'p' not found
```

We identify genes that are upregulated in each brain region in comparison to other regions.


``` r
Idents(seurat_object) <- "layer_guess"
```

``` error
Error: object 'seurat_object' not found
```

``` r
brain2 <- FindAllMarkers(seurat_object, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

``` error
Error in eval(expr, envir, enclos): object 'seurat_object' not found
```

### Spatial Differential Expression Using Moran's I

We identify the genes whose expression patterns exhibit clear spatial structure using Moran's I algorithm.


``` r
brain <- FindSpatiallyVariableFeatures(seurat_object, assay = "SCT", features = VariableFeatures(seurat_object)[1:1000], selection.method = "moransi")
```

``` error
Error in eval(expr, envir, enclos): object 'seurat_object' not found
```

## Correlation of Differnetially Expressed Genes in each Brain Region and genes with highset Moran's I value.

### Heatmap of Differential Expression


``` r
library(pheatmap)
library(dplyr)

# Unique clusters
unique_clusters <- unique(brain2$cluster)
```

``` error
Error in eval(expr, envir, enclos): object 'brain2' not found
```

``` r
# Create list of data frames filtered by cluster
cluster_data_frames <- lapply(setNames(unique_clusters, unique_clusters), function(cluster) {
  brain2 %>% filter(cluster == !!as.character(cluster))
})
```

``` error
Error in eval(expr, envir, enclos): object 'unique_clusters' not found
```

``` r
# Get and sort 'MoransI_observed' values
wer_sorted <- brain@assays[["SCT"]]@meta.features %>%
  arrange(desc(MoransI_observed)) %>%
  slice_head(n = 100)
```

``` error
Error in eval(expr, envir, enclos): object 'brain' not found
```

``` r
# Initialize p-value adjustment matrix
p_val_adj_matrix <- matrix(1, nrow = nrow(wer_sorted), ncol = length(cluster_data_frames), 
                           dimnames = list(rownames(wer_sorted), names(cluster_data_frames)))
```

``` error
Error in eval(expr, envir, enclos): object 'wer_sorted' not found
```

``` r
# Fill the matrix with adjusted p-values
for (i in rownames(wer_sorted)) {
  for (j in seq_along(cluster_data_frames)) {
    df <- cluster_data_frames[[j]]
    if (i %in% df$gene) {
      p_val_adj_value <- df$p_val_adj[df$gene == i]
      p_val_adj_matrix[i, j] <- p_val_adj_value
    }
  }
}
```

``` error
Error in eval(expr, envir, enclos): object 'wer_sorted' not found
```

``` r
# Generate and save heatmap
g <- pheatmap(p_val_adj_matrix, 
              cluster_rows = TRUE, 
              cluster_cols = TRUE, 
              display_numbers = FALSE, 
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
              main = "Heatmap of DE p-values of spatially DE genes ")
```

``` error
Error in eval(expr, envir, enclos): object 'p_val_adj_matrix' not found
```

``` r
print(g)
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'g' not found
```

The heatmap visualization reveals a key finding of our analysis: genes displaying the highest Moran's I values show distinct expression patterns that align with specific brain regions identified through expert annotations. 
This observation underscores the spatial correlation of gene expression, highlighting its potential relevance in understanding regional brain functions and pathologies.

:::::::::::::::::::::::::::::::::: keypoints

- Differential expression testing pinpoints genes with significant expression variations across samples, helping to decode biological and disease mechanisms. 
- Moran's I statistic is applied to reveal spatial autocorrelation in gene expression, critical for examining spatially dependent biological activities.
- Moran's I algorithm effectively identifies genes expressed in anatomically distinct regions, as validated from the correlation analysis with the DE genes from the annotated regions.

:::::::::::::::::::::::::::::::::::::::::::::



