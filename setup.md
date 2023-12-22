---
title: Setup
---

For this lesson, you will be working in the R programming language and the RStudio
develpment environment. You will be installing this software on your laptop and
downloading the data set.

## Software Setup
Packages: Seurat, spacexr

Installing the software may take up to 30 minutes. You may also need to contact your local 
Information Technology Help Desk to get permission or assistance installing the software. 
You may be able to install the applications using JAX Self Service software. Please do this 
before the workshop. We will not delay the start of the course while you install software. 
We will help you in advance to make sure that you have everything that you need.

If you do not already have R and RStudio installed, download and install the following software:

    R/4.3.1: Select the installation for your operating system (Windows, Mac, or Linux).
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
install.packages(c("BiocManager", "hdf5r", "Matrix", "Seurat", "remotes", "tidyverse"), dependencies = TRUE)
options(timeout = 1e6)
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

Once the installation has finished, copy and paste the following commands into the 
console to verify that both packages installed correctly.

```r
library(tidyverse)
library(Matrix)
library(hdf5r)
library(Seurat)
library(spacexr)
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
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/f9e5nshrfzk0k5bdr11h7ocjp5prrzox.h5",
              destfile = "data/151508_raw_feature_bc_matrix.h5")
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/sddrhl3ronu8nk94ja2gcifnte6lk9lt.h5",
              destfile = "data/151508_filtered_feature_bc_matrix.h5")
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/xycr1otk4hhgcbsec6vu45k2s9sisnbt.h5",
              destfile = "data/151675_raw_feature_bc_matrix.h5")
download.file(url = "https://thejacksonlaboratory.box.com/shared/static/4xoq4xcbt74zld76ifqlsagc7e36l9qs.h5",
              destfile = "data/151675_filtered_feature_bc_matrix.h5")
```

<!-- Globus link:  http://research.libd.org/globus/jhpce_HumanPilot10x/index.html -->
