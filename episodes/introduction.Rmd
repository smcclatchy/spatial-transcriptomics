---
title: Spatially Resolved Transcriptomics in Life Sciences Research
teaching:
exercises:
---

:::::: questions
 - What is spatial transcriptomics? 
 - What research questions or problems can spatial transcriptomics address?
 - How do the technologies work?
 - Which technology will we learn about in this lesson?
::::::

:::::: objectives
 - Describe why and how spatial transcriptomics can be used in research.
 - Describe how spatial transcriptomics technology works. 
 - Describe how spatial transcriptomics addresses the limitations of single-cell or bulk RNA sequencing technologies. 
::::::

## Spatial transcriptomics in biomedical research
Investigating the organization of cells and tissues is fundamental to life 
sciences research. Cells and tissues situated in different regions of an organ 
can possess diverse functions and cell types. These cells in turn are influenced 
by varying tissue microenvironments, receiving and processing distinct 
information from that microenvironment. Co-located cells can communicate 
directly with one another through chemical and mechanical signals, responding to 
these signals with changes in their own state. Thus, knowing the spatial 
organization of cells in a tissue can reveal cell and tissue function.

![A cross-section of skeletal muscle tissue showing muscle cells and a small nerve.](https://upload.wikimedia.org/wikipedia/commons/c/c3/Skeletal_muscle_-_cross_section%2C_nerve_bundle.jpg){alt='A cross-section of human skeletal muscle showing muscle cells and a nerve nearby. Stained with hematoxylin and eosin.'}

<a href="https://commons.wikimedia.org/wiki/File:Skeletal_muscle_-_cross_section,_nerve_bundle.jpg">Department of Histology, Jagiellonian University Medical College</a>
<a href="https://creativecommons.org/licenses/by-sa/3.0/deed.en" rel="license">CC BY-SA 3.0 DEED</a>

Describing spatial organization and cell signals, specifically gene expression 
signals, is the focus of spatially resolved transcriptomics. Spatial patterns of 
gene expression determine how genes are regulated within a tissue system and how 
those tissues and their component cells function. Spatial transcriptomic (ST) 
methods map cell position in a tissue, clarifying the physical relationships 
between individual cells and cellular structures. ST simultaneously measures
gene expression, delivering valuable information about cell phenotype, state, 
and cell and tissue organization and function. The combination of cellular 
expression and position sheds light on signals a cell sends or receives through 
cell-to-cell interactions or from soluble signaling molecules nearby. Spatial 
information localizes cell signaling while delivering comprehensive gene 
expression profiling within tissues. 

![Signaling between adjacent cells. The Notch protein functions as a receptor for ligands that activate or inhibit such receptors. Receptor-ligand interactions ground cell signaling and communication, often requiring close proximity between cells. ](https://upload.wikimedia.org/wikipedia/commons/0/04/Notchccr.svg){alt='alt text for accessibility purposes'}

<a href="https://commons.wikimedia.org/wiki/File:Notchccr.svg">Fred the Oyster</a> Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Spatial transcriptomics addresses a key obstacle in single-cell and bulk RNA
sequencing studies: their loss of spatial information. Spatial organization and structure determine function in most tissues and organs, so capturing both 
spatial and expression information is critical for understanding tissue function 
in neuroscience, cancer biology, developmental biology, and most other fields. Immunology is the only field that doesn't depend on spatial information to 
determine structure and function because immune cells often circulate in liquid 
media like blood and lymph. 

## Spatial transcriptomics technologies
Spatial transcriptomics technologies broadly fall within two groups: 
imaging-based and sequencing-based methods. Both imaging- and sequence-based 
datasets are available through the [The BRAIN Initiative - Cell Census Network](https://biccn.org/). Sequencing-based datasets are featured in the 
[Human Cell Atlas](https://data.humancellatlas.org/). These technologies vary in 
ability to profile entire transcriptomes, deliver single-cell resolution, and 
detect genes efficiently. 

### Imaging-based technologies
Imaging-based technologies read transcriptomes in situ 
using microscopy and feature single-cell or even single-molecule resolution. 
They identify messenger RNA (mRNA) species through hybridization with fluorescent probes. These 
probes are gene-specific in fluorescence in situ hybridization (FISH). 

![Overview of fluorescence in situ hybridization (FISH). ](fig/FISH_Overview.png){alt='a general schematic showing fluorescence in situ hybridization'}

<a href="https://commons.wikimedia.org/wiki/File:Spatial_Transcriptomics_Overview.png">
Adapted from Spatial Transcriptomics Overview by SlifertheRyeDragon.</a> Image 
created with Biorender.com. Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Current FISH methods employ multiple hybridization rounds, with risk of error 
for each transcript growing exponentially with each round. FISH methods are 
limited in the size of tissue that they can profile and most accept fresh-frozen 
(FF) tissue only. They can also be time-consuming and expensive due to 
microscopic imaging they require. Since they target specific genes, they can 
only detect genes that are in the probe set employed. They have high spatial resolution though, even delivering single-molecule resolution in single-molecule 
FISH (smFISH).

Conventional FISH methods have few distinct color channels that limit the number 
of genes that can be simultaneously analyzed. Multiplexed error-robust FISH 
(MERFISH) overcomes this problem, greatly increasing the number of RNA species
that can be simultaneously imaged in single cells by using binary code gene 
labeling in multiple rounds of hybridization.

![Schematic representation of multiplexed error-robust FISH (MERFISH). Binary codes assigned to mRNA species of interest, where “1” represents a short fluorescent DNA probe. b, Consecutive hybridization rounds, bleaching in between is implied, but not shown for clarity. At the end of the sixth round, it is possible to tell different mRNAs apart due to the decoded combinations of “1” and “0”. ](https://upload.wikimedia.org/wikipedia/commons/thumb/a/ad/MERFISH_Diagram.png/1024px-MERFISH_Diagram.png){alt='alt text for accessibility purposes'}

<a href="https://commons.wikimedia.org/wiki/File:MERFISH_Diagram.png">SlifertheRyeDragon</a>, <a href="https://creativecommons.org/licenses/by-sa/4.0">CC BY-SA 4.0</a>, via Wikimedia Commons

In situ sequencing amplifies and sequences mRNAs directly within a block 
or section of fresh-frozen (FF) or formalin-fixed paraffin embedded (FFPE) 
tissue. 

![Overview of in situ sequencing.](fig/ISS_Overview.png){alt='a general schematic showing in situ sequencing'}
<a href="https://commons.wikimedia.org/wiki/File:Spatial_Transcriptomics_Overview.png">
Adapted from Spatial Transcriptomics Overview by SlifertheRyeDragon.</a> Image 
created with Biorender.com. Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Messenger RNA (mRNA) is reverse transcribed to complementary DNA (cDNA) within 
tissue sections. A "padlock" probe binds to the cDNA, which is then 
circularized. Following circularization, the cDNA is amplified (copied) by 
rolling-circle amplification (RCA), then sequenced by ligation (joining). Probes 
profile one or two bases at a time using different fluorophores, eventually 
revealing the identity of the cDNA through imaging. Since it requires imaging,
in situ sequencing is an imaging-based method even though it involves 
sequencing. In situ sequencing can accommodate larger tissue sections than can 
FISH, though FISH methods are more efficient at detecting mRNA of genes in the 
probe set. Like FISH, in situ sequencing requires considerable imaging time on a 
microscope but delivers high spatial resolution. 

### Sequencing-based technologies
Sequencing-based methods capture, sequence, and count mRNA in situ using 
next-generation sequencing while retaining positional information. This is 
distinct from in situ sequencing because next-generation sequencing is employed, 
not imaging. The main advantage of sequenced-based methods is that they are 
unbiased in capturing mRNA because they don't rely on predefined probesets 
targeting specific genes. Sequencing-based methods retain spatial information 
through laser-capture microdissection (LCM), microfluidics, or through 
ligation of mRNAs to arrays of barcoded probes that record position. Since 
sequencing-based technologies use probes that record position instead of 
gene-specific probes, they can deliver transcriptome-wide profiling. 

LCM-based methods employ lasers to cut a tissue slice or fuse tissue to a 
membrane followed by sequencing of individual cells. 

![Overview of laser-capture microdissection .](fig/LCM_Overview.png){alt='a general schematic showing laser-capture microdissection '}
<a href="https://commons.wikimedia.org/wiki/File:Spatial_Transcriptomics_Overview.png">
Adapted from Spatial Transcriptomics Overview by SlifertheRyeDragon.</a> Image 
created with Biorender.com. Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

LCM techniques process tissue sections and individual cells for transcriptomic 
profiling by isolating regions of interest. They are useful for profiling 
transcriptomes as a first pass and for identifying RNA isoforms, but their blunt 
approach to capturing spatial expression data limits spatial resolution and 
requires many samples for sequencing. Since they focus on regions of interest, 
it is often not possible to obtain a picture of spatial expression across a 
whole tissue. LCM is an older technology that has long been used with FFPE 
tissues. Modern LCM-based approaches include Nanostring's GeoMx DSP and STRP-seq. 

Microfluidics places a chip with multiple barcode-containing channels onto a 
tissue section followed by a second chip with channels perpendicular to the 
first. The barcodes are then ligated to each other to create an array of 
unique barcodes on the tissue. This "deterministic barcoding" is employed in
DBiT-seq. DBiT-seq can be used with FFPE tissues. This approach is helpful to
avoid diffusion of mRNA away from capture areas, though a disadvantage is that
cells often sit astride multiple capture areas.

![Workflow schematic of DBiT-seq on FFPE samples.](fig/Workflow-of-DBiT-seq-on-FFPE-samples-a-Scheme-of-DBiT-seq-on-FFPE-samples.png){alt='a general schematic showing a microfluidics workflow with DBit-seq on formalin-fixed paraffin embedded (FFPE) tissue'}

<a href="https://www.researchgate.net/figure/Workflow-of-DBiT-seq-on-FFPE-samples-a-Scheme-of-DBiT-seq-on-FFPE-samples_fig2_346261659"><img src="https://www.researchgate.net/publication/346261659/figure/fig2/AS:961767006040117@1606314531966/Workflow-of-DBiT-seq-on-FFPE-samples-a-Scheme-of-DBiT-seq-on-FFPE-samples.png"/></a>
<a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Other array-based methods capture mRNA with spatially-barcoded probes and 
sequence them. They can profile larger tissue sections than can FISH or in situ 
sequencing and they don't rely on time-consuming microscopic imaging. Spatial resolution is lower, however.

In this lesson we will use data from positionally barcoded arrays.  

![A sequencing-based spatial transcriptomics method using printed spots on a slide. ](https://upload.wikimedia.org/wikipedia/commons/1/14/Spatial_transcriptomics_ii.png){alt='A graphic showing printed spots on a glass slide that are identified by a barcode and that contain primers to capture mRNA from the tissue laid on top of them'}


<a href="https://commons.wikimedia.org/wiki/User:Jasquatch">James Chell</a>, <a href="https://commons.wikimedia.org/wiki/File:Spatial_transcriptomics_ii.png">Spatial transcriptomics ii</a>, <a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" rel="license">CC BY-SA 4.0</a>

| Technology | Gene detection efficiency | Transcriptome-wide profiling | Spatial resolution | Tissue area |
| ------------ | :------: | :------: | :------: | :------: |
| FISH                | +        | -        | +        | -        |
| In situ sequencing  | -        | -        | +        | +        |
| LCM-based           |          | +        | -        | -        |
| Microfluidics       | -        | +        | -        | +        |
| Array-based         | -        | +        | -        | +        |
Table 1. Relative strengths and weaknesses of spatial transcriptomics
technologies by general category.

::::::::::::::::::::::::::::: callout

This introduction to the technologies is intended to help you 
navigate a complex technological landscape so that you can learn more about
existing technologies. It is not intended to be comprehensive.

:::::::::::::::::::::::::::::::::::::

The diversity in spatial transcriptomics technologies is enormous and rapidly
developing.  If you would like to learn more about spatial transcriptomics technologies, please see the
<a href="./reference.html"><strong>list of references located here</strong></a>.

::::::::::::::::::::::::::::::::::::: challenge 

## Discussion: Which technology is right for your research?

1. Would an imaging-based or a sequencing-based solution be preferable for your
research? Why?  

2. From the descriptions above, which technology do you think would best suit
your research? Would you use fresh-frozen (FF) or formalin-fixed tissues 
embedded in paraffin (FFPE)?Even if your institution does not offer service 
using a specific technology, describe which best suits your research and why you
think it's best suited.  

:::::::::::::::::::::::: solution 

varying - refer to strengths and weaknesses of each category, best uses of 
different technologies and how they have been used in research to date

::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

## 10X Genomics Visium technology
In this lesson we will use data from an array-based method called Visium that is
offered by 10X Genomics. Visium is an upgrade and commercialization of the 
spatial transcriptomics method described in 2016 in 
[Science, 353(6294)](https://doi.org/10.1126/science.aaf2403) 
and illustrated in general in the figure above. In brief, thin tissue sections
are placed atop spots printed with spatial barcodes. When the tissue is 
permeabilized, mRNA is released from the cells and hybridized to the spatial
barcodes. Hybridized mRNA is reverse transcribed to cDNA and then sequenced.
Spatial transcriptomics combines two key modes: histological imaging and gene 
expression profiling. Histological imaging captures tissue morphology with
standard staining protocols while expression profiling is captured by sequencing
spatially barcoded cDNA.

Sequencing-based datasets have grown faster than have imaging-based datasets, 
with Visium dominating published datasets. Unlike most other sequencing-based 
technologies, Visium accommodates both FF or FFPE tissue. Each spot provides 
average gene expression measurements from between one to ten cells, approaching
single-cell resolution. Average gene expression measurements are combined with
histological images that couple molecular detail and tissue morphology and 
structure. 

:::::: keypoints
 - Spatial transcriptomics provides the location of individual cells relative to
 neighboring cells and cell structures.
 - A cell's location is useful data for describing its phenotype, state, and
 cell and tissue function.
 - Spatial transcriptomics addresses a key obstacle in single-cell studies: the
 loss of spatial information when tissues are dissociated.
 - The main goal of spatial transcriptomics studies is to integrate expression
 with spatial information.

::::::
