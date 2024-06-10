---
title: Setup
---

Spatially resolved transcriptomics is a cutting-edge molecular profiling 
technology that measures transcript expression and position within a tissue.
This course will teach you how to analyze spatial transcriptomics data generated 
by the [Visium 10X platform][visium]. By the end of this course, you will be 
able to:

- Describe different spatial transcriptomics technologies and the best uses of 
each.
- Identify important elements for good experimental design.
- Analyze spatial transcriptomics data, including quality control, transcript quantification, and cell-type assignment within regions.

## Prerequisites
To succeed in this course, you need to have

- proficiency in the R programming language;
- knowledge of bulk RNA and single-cell sequence analysis.

For this lesson, you will be working in the R programming language and the 
RStudio development environment. You will be installing this software on your 
laptop and downloading the data set. Installing software and downloading data 
may take 2-3 hours. Please take care of this before the workshop so that you are 
able to participate at the start of the course.

## Software Setup
Please download and install R version 4.4.0 (Puppy Cup). To interact with R, 
we use RStudio. You can also download the latest stable version of RStudio, 
although this is not as critical as the latest R version is. If you don't have
administrative rights to your laptop, please ask the IT help desk to install 
software for you. Once you have installed R and RStudio, open RStudio to verify 
that the installation was successful.

R Library Installation

Next, we will install the required packages for this lesson. Note that the 
[spacexr](https://github.com/dmcable/spacexr) package takes a long time to 
download.

In RStudio, copy and paste the following commands into the Console:

```r
pkgs <- c("BiocManager", "data.table",  "ggExtra", "hdf5r",
          "here",        "igraph",      "leiden", "Matrix", "matrixStats", 
          "plyr",        "rcartocolor", "remotes",
          "Rfast2",      "Seurat",      "tidyverse", "R.utils")
for(pkg in pkgs) {
  if(!require(pkg, character.only=TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
BiocManager::install(c("glmGamPoi", "rhdf5", "ComplexHeatmap"))

options(timeout = 1e6)
remotes::install_github("immunogenomics/presto",  build_vignettes = FALSE)
remotes::install_github("dmcable/spacexr",        build_vignettes = FALSE)
```

Once the installation has finished, copy and paste the following commands into 
the console to verify that both packages installed correctly.

```r
library(BiocManager)
library(data.table)
library(R.utils)
library(hdf5r)
library(here)
library(presto)
library(rcartocolor)
library(remotes)
library(spacexr)
library(sparseMatrixStats)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
```

## Project Setup

1. Create a new project called `spatialRNA`.
  Click the File menu button, then New Project.
  Click New Directory.
  Click New Project.
  Type spatialRNA as the directory name. Create the project anywhere you like, 
  but don't forget where you put it!
  Click the Create Project button. This will create a file called 
  `spatialRNA.Rproj` in the directory you just created. In the future you can 
  double-click on this file to open RStudio in this directory. This will be the 
  easiest way to interact with the files/code you produce in this workshop.

2. Use the Files tab to create a `data` folder to hold the data, a `scripts` 
folder to house your scripts, and a `results` folder to hold results. 
Alternatively, you can copy and paste the following commands into the R console 
for step 2 only. You still need to create a project with step 1.

```r
dir.create("data")
dir.create("scripts")
dir.create("results")
```

## Data Set Download

We will be working with brain data from 
[Maynard et al., Nature Neuroscience, 2021](https://www.nature.com/articles/s41593-020-00787-0). 
We have created links on Box from which you can download the data. Once you have 
opened your `spatialRNA` project in RStudio, run the code below to download the 
data into your `data` directory.

```r
dir.create("data/151508/spatial", recursive = TRUE)
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/vfbgloxx9ciu04hj9i9f5jjzglfrx0by.h5",
              destfile = "data/151508/151508_raw_feature_bc_matrix.h5",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/puetvwocuf14kzyogtds7y1o7cme8dvp.h5",
              destfile = "data/151508/151508_filtered_feature_bc_matrix.h5",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/zxulmiupeurvrmwe0tz1y0hslzrsj4wb.json",
              destfile = "data/151508/spatial/scalefactors_json.json",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/w1bjdguo3emfb5uvwlys26yrudw08dm0.png",
              destfile = "data/151508/spatial/tissue_hires_image.png",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/xgmx8tqdfjndejr1hp3r29r4531rd6s2.png",
              destfile = "data/151508/spatial/tissue_lowres_image.png",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/gw2e7d47tihg25df8hahelrhx42y345o.csv",
              destfile = "data/151508/spatial/tissue_positions_list.csv",
              mode     = "wb")

dir.create("data/151673/spatial", recursive = TRUE)
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/ge38lg6u1i45n3grusrdyd3nccukl489.h5",
              destfile = "data/151673/151673_raw_feature_bc_matrix.h5", mode = "wb")
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/m8btvh1y9tjszfal99k2cvr1f32lh1si.h5",
              destfile = "data/151673/151673_filtered_feature_bc_matrix.h5", mode = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/feb8gnawor51ojh2ci4mlshhxobpnaji.json",
              destfile = "data/151673/spatial/scalefactors_json.json",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/ejyx4qkv62p5t0njwf5z8px2mqcjxnd7.png",
              destfile = "data/151673/spatial/tissue_hires_image.png",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/8tgbt4654zbxwsqr3vwlk64dyzzulhlb.png",
              destfile = "data/151673/spatial/tissue_lowres_image.png",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/drlayml5otq7n2xedndm0qsqly58g306.csv",
              destfile = "data/151673/spatial/tissue_positions_list.csv",
              mode     = "wb")


download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/ny1wokl6sz1xjzz68aftbk209se5nvws.tsv",
              destfile = "data/spot-meta.tsv",
              mode     = "wb")

dir.create("data/scRNA-seq", recursive = TRUE)
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/ydu9rbdhum5qrvuijze23qwz7dlztefo.tsv",
              destfile = "data/scRNA-seq/sc_cell_types.tsv",
              mode     = "wb")
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/lasxuiq5wi3ms1jnokzmm7pp4hptr8ma.gz",
              destfile = "data/scRNA-seq/sc_counts.tsv.gz",
              mode     = "wb")

download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/dt2chlmxtjajxfnlpfzolz1tvb6kic7p.rds",
              destfile = "data/rctd-sample-1.rds",
              mode     = "wb")
```

## Session Info

```r
sessionInfo()
```

<!-- Globus link:  http://research.libd.org/globus/jhpce_HumanPilot10x/index.html -->

Development of this lesson was funded by a [Jackson Laboratory](https://www.jax.org/) 
Director's Innovation Fund award to Dr. Gary Churchill. Lesson authors are 
grateful for this support.
