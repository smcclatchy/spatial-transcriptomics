---
title: Setup
---

For this lesson, you will be working in the R programming language and the RStudio
develpment environment. You will be installing this software on your laptop and
downloading the data set.

## Software Setup

Installing the software may take up to 30 minutes. You may also need to contact your local 
Information Technology Help Desk to get permission or assistance installing the software. 
You may be able to install the applications using JAX Self Service software. Please do this 
before the workshop. We will not delay the start of the course while you install software. 
We will help you in advance to make sure that you have everything that you need.

If you do not already have R and RStudio installed, download and install the following software:

    R/4.3.3: Select the installation for your operating system (Windows, Mac, or Linux).
    RStudio: Download the free Rstudio Desktop.

You do not need to install this exact version of R, but it would be good to make sure your R 
is relatively recent (say, updated within the past year).

Once you have installed R and RStudio, open RStudio to verify that the installation was 
successful.
R Library Installation

Next, we will install the required packages for this lesson. Note that the [spacexr](https://github.com/dmcable/spacexr) 
package takes a long time to download. It may take up to 30 minutes for it to install.

In RStudio, copy and paste the following commands into the Console:

```r
install.packages(c("BiocManager", "data.table",  "doMC",  "ggExtra", "hdf5r",
                   "here",        "igraph",      "leiden", "Matrix", "matrixStats", 
                   "plyr",        "rcartocolor", "remotes",
                   "Rfast2",      "Seurat",      "tidyverse"), dependencies = TRUE)
BiocManager::install(c("glmGamPoi", "rhdf5"))

options(timeout = 1e6)
remotes::install_github("immunogenomics/harmony", build_vignettes = FALSE)
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

Once the installation has finished, copy and paste the following commands into the 
console to verify that both packages installed correctly.

```r
library(BiocManager)
library(data.table)
#library(doMC)
#library(ggExtra)
library(hdf5r)
library(here)
#library(igraph)
#library(leiden)
#library(Matrix)
#library(matrixStats)
#library(plyr)
library(rcartocolor)
library(remotes)
#library(Rfast2)
library(spacexr)
library(Seurat)
library(tidyverse)
```

## Project Setup

    Create a new project called "spatialRNA".
        Click the File menu button, then New Project.
        Click New Directory.
        Click New Project.
        Type spatialRNA as the directory name. Create the project anywhere you like, but don't 
        forget where you put it!
        Click the Create Project button. This will create a file called "spatialRNA.Rproj" in the 
        directory you just created. In the future you can double-click on this file to open 
        RStudio in this directory. This will be the easiest way to interact with the files/code 
        you produce in this workshop.

    Use the Files tab to create a data folder to hold the data, a scripts folder to house 
    your scripts, and a results folder to hold results. Alternatively, you can copy and paste 
    the following commands into the R console for step 2 only. You still need to create a 
    project with step 1.

```r
dir.create("data")
dir.create("scripts")
dir.create("results")
```

## Data Set Download

We will be working with brain data from 
[Maynard et al., Nature Neuroscience, 2021](https://www.nature.com/articles/s41593-020-00787-0). 
We have created links on Box from which you can download the data. Once you have opened your
"spatialRNA" project in RStudio, run the code below to download the data into your "data" 
directory.

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
```

## Session Info

```r
sessionInfo()
```

<!-- Globus link:  http://research.libd.org/globus/jhpce_HumanPilot10x/index.html -->
