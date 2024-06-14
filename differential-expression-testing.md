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

## Differential Expression Analysis

### Differential Expression Using Expert's Annotation

We will begin by performing differential expression across the annotated layers. 
We will use the source publication's spot annotation, which is stored in the
"layer_guess" column of the metadata. As a reminder, those look like:


``` r
SpatialDimPlotColorSafe(filter_st[, !is.na(filter_st[[]]$layer_guess)], "layer_guess") + 
  labs(fill = "Layer") 
```

<img src="fig/differential-expression-testing-rendered-layers-1.png" style="display: block; margin: auto;" />

We identify genes that are upregulated in each annotated brain region using the 
[FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) function 
in Seurat. This performs a "one versus rest" comparison of a gene's expression 
in one region relative to that gene's expression in all other regions. The 
default test used here, the Wilcoxon Rank Sum test, 
will use an efficient implementation within the 
[presto](https://github.com/immunogenomics/presto) library, if installed. The 
speedup over the default implementation is substantial, and we highly recommend 
installing presto and using Seurat v5, which leverages it.


``` r
Idents(filter_st)  <- "layer_guess"
de_genes           <- FindAllMarkers(filter_st, 
                                     assay    = "SCT",
                                     verbose  = FALSE,
                                     only.pos = TRUE, 
                                     min.pct  = 0.25, 
                                     logfc.threshold = 0.25)
```

## Moran's I Statistic

Moran's I is a measure used to assess spatial autocorrelation in data, 
indicating whether similar values of a feature (e.g., expression levels of
a gene) are clustered, dispersed, or random 
([Jackson, et al. 2010](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/1476-072X-9-33)). 
These correspond to Moran's I values that are positive, negative, or near zero, 
respectively.

Here, we can apply it to detect genes whose expression patterns 
exhibit spatial structure, which may reflect region-specific, biological function.
That is, we anticipate that spatially variable genes will exhibit region-specific
expression. Let's check that hypothesis by first computing spatially variable genes
and then assessing whether they are differentially expressed across regions.

![Moran's I statistic quantifies spatial correlation. **Top Left:** Checkerboard pattern results in negative Moran's I, indicating anti-correlation. **Top Right:** Linear gradient shows a high positive Moran's I, reflecting a strong spatial gradient. **Bottom Left:** Random pattern leads to a Moran's I near zero, suggesting no significant spatial autocorrelation. **Bottom Right:** 'Ink blot' pattern demonstrates positive autocorrelation, indicative of a clustered or spreading pattern. Relationships are calculated using direct, equally weighted neighbors, normalized for each cell. ](https://upload.wikimedia.org/wikipedia/commons/f/f0/Moran%27s_I_example.png){alt="Moran's I statistic quantifies spatial correlation."}

Image by <a href="https://commons.wikimedia.org/wiki/File:Moran%27s_I_example.png">WikiNukalito</a>, <a href="https://creativecommons.org/licenses/by-sa/4.0">CC BY-SA 4.0</a>, via Wikimedia Commons.

### Spatial Differential Expression Using Moran's I

We identify the genes whose expression patterns exhibit clear spatial structure 
using Moran's I algorithm. We have selected the top 1,000 genes, which should
be sufficient to identify brain regions.


``` r
svg <- 
  FindSpatiallyVariableFeatures(filter_st, 
                                assay            = "SCT", 
                                features         = VariableFeatures(filter_st)[1:1000], 
                                selection.method = "moransi")
```

## Correlation of Differentially Expressed Genes in each Brain Region and genes with highset Moran's I value.

### Heatmap of Differential Expression

In this case, we have spot annotation for this tissue section from the author
publication. Let's check the degree of correlation between the author's 
brain region annotation and the spatially differentially expressed genes. To do
this, we will plot a heatmap of the p-values of the most spatially 
differentially expressed genes, organized by brain region.


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
Heatmap(p_val_adj_matrix,
        column_title      = "Heatmap of DE p-values of spatially DE genes",
        name              = "DE p-values", # Title for the heatmap legend
        row_title         = "Spatially variable genes",
        cluster_rows      = TRUE, 
        cluster_columns   = FALSE,
        show_row_names    = FALSE, 
        show_column_names = TRUE,
        show_row_dend     = FALSE, 
        show_column_dend  = TRUE)
```

<img src="fig/differential-expression-testing-rendered-heatmap-de-1.png" style="display: block; margin: auto;" />

The heatmap visualization reveals a key finding of our analysis: genes displaying 
the highest Moran's I values show distinct expression patterns that align with 
specific brain regions identified through expert annotations. This observation 
underscores the spatial correlation of gene expression, highlighting its 
potential relevance in understanding regional brain functions and pathologies.

:::::::::::::::::::::::::::::::::: keypoints

- Differential expression testing pinpoints genes with significant expression 
variations across samples, helping to decode biological and disease mechanisms. 
- Moran's I statistic is applied to reveal spatial autocorrelation in gene 
expression, critical for examining spatially dependent biological activities.
- Moran's I algorithm effectively identifies genes expressed in anatomically 
distinct regions, as validated from the correlation analysis with the DE genes 
from the annotated regions.

:::::::::::::::::::::::::::::::::::::::::::::



