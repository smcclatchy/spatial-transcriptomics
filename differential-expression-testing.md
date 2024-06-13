---
title: 'Differential Expression Testing'
teaching: 60
exercises: 20
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is the purpose of differential expression testing in bioinformatics?
- Can Moran's I algorithm independently identify region-specific differential
expressions that align with results obtained from expert annotations?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Identify differentially expressed genes across different layers using expert 
annotations.
- Utilize Moran's I algorithm to find spatially variable genes.
- Explore the correlation between genes identified through expert annotations 
and those detected by Moran's I algorithm.

::::::::::::::::::::::::::::::::::::::::::::::::



## Introduction to Differential Expression Testing

Differential expression testing is crucial in bioinformatics for identifying 
genes that show significant differences in expression across different samples 
or groups. 
This method helps find genes that are upregulated or downregulated in specific 
contexts, providing insights into biological functions and disease mechanisms.

## Moran's I Statistic

Moran's I is a measure used to assess spatial autocorrelation in data, 
indicating whether similar values are clustered, dispersed, or random. 
In bioinformatics, it's applied to detect genes whose expression patterns 
exhibit clear spatial structure, aiding in understanding spatially localized 
biological processes.

## Differential Expression Analysis

### Differential Expression Using Expert's Annotation

We will begin by performing differential expression across the annotated layers. 
As a reminder, those look like:


``` r
SpatialDimPlotColorSafe(filter_st[, !is.na(filter_st[[]]$layer_guess)], "layer_guess") + 
  labs(fill = "Layer") 
```

<img src="fig/differential-expression-testing-rendered-layers-1.png" style="display: block; margin: auto;" />

We identify genes that are upregulated in each brain region in comparison to 
other regions.


``` r
Idents(filter_st)  <- "layer_guess"
de_genes           <- FindAllMarkers(filter_st, 
                                     assay    = "SCT", 
                                     only.pos = TRUE, 
                                     min.pct  = 0.25, 
                                     logfc.threshold = 0.25)
```

``` output
Calculating cluster Layer3
```

``` output
Calculating cluster Layer1
```

``` output
Calculating cluster WM
```

``` output
Calculating cluster Layer5
```

``` output
Calculating cluster Layer6
```

``` output
Calculating cluster Layer2
```

``` output
Calculating cluster Layer4
```

### Spatial Differential Expression Using Moran's I

We identify the genes whose expression patterns exhibit clear spatial structure 
using Moran's I algorithm.

> DMG: Can one of you elaborate on what Moran's I is? And possibly add a reference?


``` r
svg <- 
  FindSpatiallyVariableFeatures(filter_st, 
                                assay            = "SCT", 
                                features         = VariableFeatures(filter_st)[1:1000], 
                                selection.method = "moransi")
```

## Correlation of Differentially Expressed Genes in each Brain Region and genes with highset Moran's I value.

### Heatmap of Differential Expression


``` r
# Get and sort 'MoransI_observed' values
morans_i_genes <- svg@assays[["SCT"]]@meta.features %>%
                    rownames_to_column("gene") %>%
                    arrange(desc(MoransI_observed)) %>%
                    slice_head(n = 100)

# Merge the Moran's I values with the DE genes
df <- merge(morans_i_genes, de_genes, all.x = TRUE, by = "gene")
df <- subset(df, !is.na(cluster))

# Create a matrix whose rows are the spatially variable genes (indicated by Moran's I),
# whose columns are the clusters, and whose entries are the adjusted DE pvalue for the
# corresponding gene and cluster.
p_val_adj_matrix <- df %>%
                       select(gene, cluster,p_val_adj) %>%
                       pivot_wider(names_from = cluster, values_from = p_val_adj, values_fill = 1.0) %>%
                       column_to_rownames("gene") %>%
                       as.matrix()

# Create a heatmap of the DE pvalues of spatially variable genes
Heatmap(p_val_adj_matrix2,
        column_title      = "Heatmap of DE p-values of spatially DE genes",
        name              = "DE p-values", # Title for the heatmap legend
        row_title         = "Spatially variable genes",
        cluster_rows      = TRUE, 
        cluster_columns   = TRUE,
        show_row_names    = FALSE, 
        show_column_names = TRUE,
        show_row_dend     = FALSE, 
        show_column_dend  = TRUE)
```

``` error
Error in eval(expr, envir, enclos): object 'p_val_adj_matrix2' not found
```

The heatmap visualization reveals a key finding of our analysis: genes displaying the highest Moran's I values show distinct expression patterns that align with specific brain regions identified through expert annotations. 
This observation underscores the spatial correlation of gene expression, highlighting its potential relevance in understanding regional brain functions and pathologies.

:::::::::::::::::::::::::::::::::: keypoints

- Differential expression testing pinpoints genes with significant expression variations across samples, helping to decode biological and disease mechanisms. 
- Moran's I statistic is applied to reveal spatial autocorrelation in gene expression, critical for examining spatially dependent biological activities.
- Moran's I algorithm effectively identifies genes expressed in anatomically distinct regions, as validated from the correlation analysis with the DE genes from the annotated regions.

:::::::::::::::::::::::::::::::::::::::::::::



