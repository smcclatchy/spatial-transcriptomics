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
sciences research. Tissues in different regions of an organ can possess diverse
functions and cell types and are influenced by varying tissue microenvironments. 
Spatial patterns of gene expression determine how genes are regulated within a 
tissue system and how those tissues and their component cells function. 
Spatially resolved transcriptomic methods map positions of individual cells in a 
tissue while simultaneously measuring their expression, delivering valuable 
information about cell phenotype, state and cell and tissue organization and 
function. Cell position  clarifies the physical relationships between individual 
cells and cellular structures. When combined with cellular expression, cell
position sheds light on signals a cell sends or receives through cell-to-cell 
interactions or from soluble signaling molecules nearby. Spatial information 
localizes cell signaling while delivering comprehensive gene expression 
profiling within tissues. 

![Signaling between adjacent cells. The Notch protein functions as a receptor for ligands that activate or inhibit such receptors. Receptor-ligand interactions ground cell signaling and communication, often requiring close proximity between cells. ](https://upload.wikimedia.org/wikipedia/commons/0/04/Notchccr.svg){alt='alt text for accessibility purposes'}

<a href="https://commons.wikimedia.org/wiki/File:Notchccr.svg">Fred the Oyster</a> Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Spatial transcriptomics addresses a key obstacle in single-cell and bulk RNA
sequencing studies: the loss of spatial information when tissues are 
dissociated. Spatial organization and structure determine function in most
tissues and organs, so capturing both spatial and expression information is 
critical for understanding tissue function in neuroscience, cancer biology, 
developmental biology, and most other fields. Immunology is the only field that 
doesn't depend on spatial information to determine structure and function 
because immune cells often circulate in liquid media like blood and lymph. 

## Spatial transcriptomics technologies
Spatial transcriptomics technologies broadly fall within two groups: 
imaging-based and sequencing-based methods. These technologies vary in ability 
to profile entire transcriptomes, deliver single-cell resolution, and detect
genes efficiently. Imaging-based technologies read transcriptomes in situ using 
microscopy and feature single-cell or even single-molecule resolution. They 
identify mRNA species through hybridization with fluorescent probes. These 
probes are gene-specific in fluorescence in situ hybridization (FISH). Current 
FISH methods employ multiple hybridization rounds, with risk of error for each 
transcript growing exponentially with each round. FISH methods are limited in 
the size of tissue that they can profile and can be time-consuming due to 
microscopic imaging they require. They also can only detect genes that are in
the probe set employed. They have high spatial resolution though, even
delivering single-molecule resolution in single-molecule FISH (smFISH).

In situ sequencing amplifies and sequences mRNAs directly within a block 
or section of tissue. Probes profile one or two bases at a time using different
fluorophores, eventually revealing the identity of the mRNA through imaging. In 
situ sequencing can accommodate larger tissue sections than can FISH, though 
FISH methods are more efficient at detecting mRNA of genes in the probe set. 
Like FISH, in situ sequencing requires considerable imaging time on a microscope 
but delivers high spatial resolution. 

Sequencing-based methods capture, sequence, and count mRNA in situ using 
next-generation sequencing while retaining positional information. This is 
distinct from in situ sequencing because next-generation sequencing is employed, 
not imaging. Sequencing-based methods retain spatial information through laser-
capture microdissection (LCM) and microfluidics, or through ligation of mRNAs to 
arrays of barcoded probes that record position. LCM employs lasers to cut a 
tissue slice or fuse tissue to a membrane followed by sequencing of individual
cells. Microfluidics places a chip with multiple barcode-containing channels 
onto a tissue section followed by a second chip with perpendicular channels and 
a new set of barcodes. Since sequencing-based technologies use probes that 
record position instead of using gene-specific probes, they can deliver 
transcriptome-wide profiling. 

LCM techniques process tissue sections and individual cells for transcriptomic 
profiling by isolating regions of interest. They are useful for profiling 
transcriptomes as a first pass and for identifying RNA isoforms, but their blunt 
approach to capturing spatial expression data limits spatial resolution and 
requires many samples for sequencing. Since they focus on regions of interest,
it is often not possible to obtain a picture of spatial expression across a 
whole tissue. Array-based methods capture mRNA with spatially-barcoded probes 
and sequence them, or alternatively print the array onto tissue using 
microfluidic channels to produce a unique barcode at every position in the 
tissue. Array-based methods can profile larger tissue sections than can FISH or 
in situ sequencing and they don't rely on microscopic imaging, which can be 
quite time consuming. Spatial resolution is lower, however.

| Technology | Gene detection efficiency | Transcriptome-wide profiling | 
Spatial resolution | Tissue area |
| ------------ | :------: | :------: | :------: | :------: |
| FISH                | +        | -        | +        | -        |
| In situ sequencing  | -        | -        | +        | +        |
| LCM                 |          | +        | -        | -        |
| Microfluidics       |          | +        | -        | +        |
| Array-based         | -        | +        | -        | +        |
Table 1. Relative strengths and weaknesses of spatial transcriptomics
technologies by general category.

Array-based methods use tissue frozen below the temperature at which RNA 
degrades (fresh-frozen; FF), though some array-based methods like ... can use 
formalin-fixed tissues embedded in paraffin (FFPE). 

::::::::::::::::::::::::::::: callout

The diversity in spatial transcriptomics technologies is enormous and rapidly
developing. This introduction to the technologies is intended to help you 
navigate a complex technological landscape so that you can learn more about
existing technologies. It is not intended to be comprehensive. If you would
like to learn more about spatial transcriptomics technologies, please see the
list of references located <a href="./reference.html"><strong>here</strong></a>.

:::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge 

## Discussion: Which technology is right for your research?

1. Would an imaging-based or a sequencing-based solution be preferable for your
research? Why?

:::::::::::::::::::::::: solution 

varying - refer to strengths and weaknesses of each category

::::::::::::::::::::::::

2. From the descriptions above, which technology do you think would best suit
your research? Even if your institution does not offer service using a 
specific technology, describe which best suits your research and why you think
it's best suited. Would you use fresh-frozen (FF) or formalin-fixed tissues 
embedded in paraffin (FFPE)? Why is the technology best suited to your research?

:::::::::::::::::::::::: solution 

varying - describe best uses of different technologies and how they have been
used in research to date

::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

In this lesson we will use data from positionally barcoded arrays.  

![A sequencing-based spatial transcriptomics method using printed spots on a slide. ](https://upload.wikimedia.org/wikipedia/commons/1/14/Spatial_transcriptomics_ii.png){alt='alt text for
accessibility purposes'}
<a href="https://commons.wikimedia.org/wiki/User:Jasquatch">James Chell</a>, <a href="https://commons.wikimedia.org/wiki/File:Spatial_transcriptomics_ii. png">Spatial transcriptomics ii</a>, <a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" rel="license">CC BY-SA 4.0</a>

## 10X Genomics Visium technology

:::::: keypoints
 - Spatial transcriptomics provides the location of individual cells relative to neighboring cells and cell structures.
 - A cell's location is useful data for describing its phenotype, state, and cell and tissue function.
 - Spatial transcriptomics addresses a key obstacle in single-cell studies: the loss of spatial information when tissues are dissociated.
::::::
