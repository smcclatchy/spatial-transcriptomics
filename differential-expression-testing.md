---
title: 'Differential Expression Testing'
teaching: 60
exercises: 20
---

:::::::::::::::::::::::::::::::::::::: questions 

- How can we access region-specific gene expression using differential expression?
- How can we access spatially varying gene expression using spatial statistics?
::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Identify differentially expressed genes across different layers defined by expert 
annotations.
- Utilize the Moran's I statistic to find spatially variable genes.
- Explore the relationship between layer-specific genes identified by differential expression and spatially varying genes identified by Moran's I.

::::::::::::::::::::::::::::::::::::::::::::::::



## Spot-level Differential Expression Using Expert Annotation

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

The resulting table indicates DE p-values and adjusted p-values for each gene,
along with the percentage of spots in which the gene was detected (pct.1) 
in the corresponding cluster and the percentage of spots in which it was
detected in all <em>other</em> clusters (pct.2):


``` r
head(de_genes)
```

``` output
                p_val avg_log2FC pct.1 pct.2     p_val_adj cluster    gene
MT-CO1  7.105768e-163  0.3732786 1.000 1.000 1.278754e-158  Layer3  MT-CO1
ENC1    1.654660e-142  0.8693181 0.998 0.929 2.977725e-138  Layer3    ENC1
MT-ATP6 6.825005e-127  0.3333658 1.000 1.000 1.228228e-122  Layer3 MT-ATP6
MT-CYB  8.895560e-118  0.3290618 1.000 1.000 1.600845e-113  Layer3  MT-CYB
MT-CO2  3.349139e-116  0.2869993 1.000 1.000 6.027110e-112  Layer3  MT-CO2
MT-CO3  5.357942e-105  0.2766345 1.000 1.000 9.642152e-101  Layer3  MT-CO3
```

It is not uncommon to have pathologist annotations of regions. The Visium assay is performed
on a tissue that is also stained for hematoxylin and eosin (H&E) -- a routine practice in
pathology for diagnosing cancer, for example. A pathologist could manually annotate this H&E
image using a histology viewer, such as QuPath.

In cases where we do not have expert annotations, we could perform the analysis above across
clusters. Indeed, this is often the first step to interpreting a cluster -- comparing its marker
genes to those of regions expected within the tissue. It is important to remember, however,
that these markers are for clusters of spots -- aggregations of cells -- not individual cells, as
would be the case in scRNA-seq. Those aggregations may be of cells with similar types, in which
case analysis may be similar to that of scRNA-seq, or of different cell types, in which case the
interpretation would be quite different than with scRNA-seq.

## Moran's I Statistic

We next consider an alternative to the above approach that relies on defined regions -- either
through expert annotation or via clustering. Instead, we will apply the Moran's I statistic.
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
using Moran's I algorithm, as implemented in [FindSpatiallyVariableFeatures](https://satijalab.org/seurat/reference/findspatiallyvariablefeatures). 
We have selected the top 1,000 genes, which should
be sufficient to identify brain regions. The following will take several minutes to run.


``` r
svg <- 
  FindSpatiallyVariableFeatures(filter_st, 
                                assay            = "SCT", 
                                features         = VariableFeatures(filter_st)[1:1000], 
                                selection.method = "moransi")
```

## Correlation of Region-specific Differentially Expressed Genes and Spatially Variable Genes

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



