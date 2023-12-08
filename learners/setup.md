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

In RStudio, copy and paste the following commands into the Console:

```r
install.packages(c("BiocManager", "Matrix", "Seurat", "spacexr", "tidyverse"), dependencies = TRUE)
BiocManager::install("rhdf5")
```

Once the installation has finished, copy and paste the following commands into the 
console to verify that both packages installed correctly.

```r
library(tidyverse)
library(Matrix)
library(rhdf5)
library(Seurat)
library(spacexr)
```

## Project Setup

    Create a new project called scRNA.
        Click the File menu button, then New Project.
        Click New Directory.
        Click New Project.
        Type scRNA as the directory name. Create the project anywhere you like, but don't 
        forget where you put it!
        Click the Create Project button. This will create a file called scRNA.Rproj in the 
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
:::::::::::::::::::::::::

## Data Set Download

<!--
FIXME: place any data you want learners to use in `episodes/data` and then use
       a relative link ( [data zip file](data/lesson-data.zip) ) to provide a
       link to it, replacing the example.com link.
-->

Original data is on a Globus endpoint:
http://research.libd.org/globus/jhpce_HumanPilot10x/index.html

But Antonios has it hosted on Box: <https://thejacksonlaboratory.ent.box.com/s/kqo4d25nba067qcvf908t5obqiym3zak/folder/238879909912>

TBD: Create static links that don't expire for each directory.

```r
download.file(url = "???", destfile = "data")
```
