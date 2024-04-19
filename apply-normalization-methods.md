---
title: 'Normalization in Spatial Transcriptomics'
teaching: 15
exercises: 5
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do we determine the necessity of normalization in spatial transcriptomics?
- What insights do additional modalities like H&E staining provide in assessing normalization needs?
- How do specific normalization techniques like SCTransform and log scaling work and when should they be applied?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Assess the need for normalization using both spatial transcriptomics data and ancillary modalities like H&E staining.
- Understand the specific applications and mechanisms of normalization techniques such as SCTransform and log scaling.
- Implement adaptive normalization strategies that accurately reflect both absolute and relative cellular information.

::::::::::::::::::::::::::::::::::::::::::::::::

## Understanding Normalization in Spatial Transcriptomics

Normalization in spatial transcriptomics must be carefully tailored to each dataset, balancing the technical corrections with the preservation of biologically meaningful signals.

### Assessing Normalization Needs

#### Using H&E Staining to Guide Normalization
Hematoxylin and Eosin (H&E) staining is critical for preliminary assessments of tissue sections. It highlights structural and pathological features, guiding the interpretation of transcriptomic data. For example, high RNA counts in a necrotic region, typically characterized by reduced cellular material, might suggest technical artifacts, indicating a need for normalization.

#### Mean-Variance Plotting
Mean-variance plots are an essential tool for assessing gene expression variability relative to the mean expression levels across different spots. By plotting the variance of gene expression against the mean expression level, researchers can identify genes with variance that deviates significantly from what would be expected under normal biological conditions. This can be particularly useful for spotting genes that are overly influenced by technical artifacts or biological outliers, suggesting the corresponding normalization choices.

### Normalization Techniques and Their Implications

#### SCTransform
SCTransform is a normalization method for single-cell and spatial transcriptomics that uses a regularized negative binomial regression to stabilize variance across expression levels. It selects highly variable genes and corrects for technical noise by modeling gene expression counts with Pearson residuals. This approach effectively adjusts for confounding factors such as sequencing depth, facilitating more accurate downstream analyses like clustering.

#### LogNormalize
LogNormalize is a specific normalization method that scales gene expression data to account for differences in cell-specific total RNA counts. This process involves dividing the raw gene expression counts in each cell by the total counts in that cell, multiplying by a scale factor, and then applying a natural logarithm transformation using `log1p` (log(1+x)). This method helps in reducing the skewness caused by highly expressed genes and stabilizes the variance across the dataset, making it more suitable for downstream analytical comparisons.

### No One-Size-Fits-All Approach
Raw read counts provide essential insights into the absolute cell type densities within a sample, which are crucial for mapping cellular distribution. In contrast, normalized data adjusts for technical variations like sequencing depth and RNA capture efficiency, thus revealing the relative proportions of cell types and identifying specific tissue structures, such as epithelium or fibrosis. Saiselet et al. demonstrated that while normalized data effectively identify distinct morphologies, raw counts are vital for detecting areas with unusual cell-type concentrations, such as high epithelial regions expressing vimentin (VIM) (Saiselet, M., et al., Journal of Molecular Cell Biology, 2020).
Hence, choosing the right normalization method depends on the specific characteristics of each dataset and the biological questions at hand. Researchers must understand the impact of each normalization strategy on both the biological and technical aspects of their data.



## Why Normalize?

### Description of normalization techniques

### Data Interpretation

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Analyze Normalized vs. Raw Data

Consider two datasets from the same spatial transcriptomics experiment: one raw and one normalized. What differences would you expect in the analysis results between these two datasets?

:::::::::::::::::::::::: solution 

Normalized data should show reduced variance due to technical factors and enhanced clarity in detecting biological differences. Raw data may exhibit skewed expression profiles influenced by technical artifacts.


::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::

## Normalization Methods

::::::::::::::::::::::::::::::::::::: keypoints 

- Normalization is essential but must be selectively applied based on the unique characteristics of each dataset and the specific biological questions at hand.
- Techniques like SCTransform and log scaling offer ways to balance technical correction with biological integrity.
- Examining both raw and normalized data can provide comprehensive insights into the absolute and relative characteristics of cellular components in spatial transcriptomics.

::::::::::::::::::::::::::::::::::::::::::::::::


