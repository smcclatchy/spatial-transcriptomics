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

## Differential Expression Analysis

### Differential Expression Using Expert's Annotation

We will begin by performing differential expression across the annotated layers. 
As a reminder, those look like:


``` r
SpatialDimPlotColorSafe(filter_st[, !is.na(filter_st[[]]$layer_guess)], "layer_guess") + labs(fill="Layer") 
```

``` warning
Warning: Not validating Centroids objects
Not validating Centroids objects
```

``` warning
Warning: Not validating FOV objects
Not validating FOV objects
Not validating FOV objects
Not validating FOV objects
Not validating FOV objects
Not validating FOV objects
```

``` warning
Warning: Not validating Seurat objects
```

``` output
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
```

<img src="fig/differential-expression-testing-rendered-layers-1.png" style="display: block; margin: auto;" />

We identify genes that are upregulated in each brain region in comparison to other regions.


``` r
Idents(seurat_object)  <- "layer_guess"
```

``` error
Error: object 'seurat_object' not found
```

``` r
brain2                  <- FindAllMarkers(filter_st, 
                                          assay = "SCT", 
                                          only.pos = TRUE, 
                                          min.pct = 0.25, 
                                          logfc.threshold = 0.25)
```

``` output
Calculating cluster 0
```

``` output
For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
(default method for FindMarkers) please install the presto package
--------------------------------------------
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
--------------------------------------------
After installation of presto, Seurat will automatically use the more 
efficient implementation (no further action necessary).
This message will be shown once per session
```

``` output
Calculating cluster 1
```

``` output
Calculating cluster 2
```

``` output
Calculating cluster 3
```

``` output
Calculating cluster 4
```

``` output
Calculating cluster 5
```

``` output
Calculating cluster 6
```

``` output
Calculating cluster 7
```

``` output
Calculating cluster 8
```

### Spatial Differential Expression Using Moran's I

We identify the genes whose expression patterns exhibit clear spatial structure using Moran's I algorithm.


``` r
brain <- FindSpatiallyVariableFeatures(filter_st, 
                                       assay            = "SCT", 
                                       features         = VariableFeatures(filter_st)[1:1000], 
                                       selection.method = "moransi")
```

``` output
Computing Moran's I
```

``` warning
Warning in dist(x = pos): NAs introduced by coercion
```

``` output
Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

``` output
Also defined by 'BiocGenerics'
```

``` output
Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

``` output
Also defined by 'BiocGenerics'
```

## Correlation of Differentially Expressed Genes in each Brain Region and genes with highset Moran's I value.

### Heatmap of Differential Expression


``` r
# Unique clusters
unique_clusters        <- unique(brain2$cluster)

# Create list of data frames filtered by cluster
cluster_data_frames    <- lapply(
                            setNames(unique_clusters, unique_clusters), 
                            function(cluster) {
                              brain2 %>% filter(cluster == !!as.character(cluster))
                            }
                          )

# Get and sort 'MoransI_observed' values
wer_sorted             <- brain@assays[["SCT"]]@meta.features %>%
                            arrange(desc(MoransI_observed)) %>%
                            slice_head(n = 100)

# Initialize p-value adjustment matrix
p_val_adj_matrix       <- matrix(
                            1, 
                            nrow = nrow(wer_sorted), 
                            ncol = length(cluster_data_frames), 
                            dimnames = list(
                              rownames(wer_sorted), 
                              names(cluster_data_frames)
                            )
                          )

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

# Define colors using a color ramp
color_palette <- colorRamp2(c(0, 0.5, 1), c("navy", "white", "firebrick3"))
```

``` error
Error in colorRamp2(c(0, 0.5, 1), c("navy", "white", "firebrick3")): could not find function "colorRamp2"
```

``` r
# Create the heatmap
heatmap <- Heatmap(p_val_adj_matrix,
                   name              = "DE p-values", # Title for the heatmap legend
                   cluster_rows      = TRUE, 
                   cluster_cols      = TRUE,
                   show_row_names    = FALSE, 
                   show_column_names = FALSE,
                   show_row_dend     = TRUE, 
                   show_column_dend  = TRUE,
                   cell_fun          = function(j, i, x, y, width, height, fill) {
                                       if (!is.na(p_val_adj_matrix[i, j])) {
                                           grid.text(sprintf("%.2f", p_val_adj_matrix[i, j]), x, y)
                                       }
                                   }, # Optionally display numbers in cells
                   col               = color_palette)
```

``` error
Error in Heatmap(p_val_adj_matrix, name = "DE p-values", cluster_rows = TRUE, : unused argument (cluster_cols = TRUE)
```

``` r
# Draw the heatmap
draw(heatmap, main = "Heatmap of DE p-values of spatially DE genes")
```

``` error
Error: unable to find an inherited method for function 'draw' for signature 'object = "function"'
```

The heatmap visualization reveals a key finding of our analysis: genes displaying the highest Moran's I values show distinct expression patterns that align with specific brain regions identified through expert annotations. 
This observation underscores the spatial correlation of gene expression, highlighting its potential relevance in understanding regional brain functions and pathologies.

:::::::::::::::::::::::::::::::::: keypoints

- Differential expression testing pinpoints genes with significant expression variations across samples, helping to decode biological and disease mechanisms. 
- Moran's I statistic is applied to reveal spatial autocorrelation in gene expression, critical for examining spatially dependent biological activities.
- Moran's I algorithm effectively identifies genes expressed in anatomically distinct regions, as validated from the correlation analysis with the DE genes from the annotated regions.

:::::::::::::::::::::::::::::::::::::::::::::



