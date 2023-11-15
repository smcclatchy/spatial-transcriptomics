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
sciences research. Spatially resolved transcriptomic methods map positions of individual cells in a tissue and measure their expression, delivering valuable
information about cell and tissue organization and function. Cell position and
expression can describe a cell's phenotype, state, and cell and tissue function
while clarifying the physical relationships between individual cells and 
cellular structures. A cell's position sheds light on signals it receives 
through cell-to-cell interactions or from dissolved signaling molecules nearby. 
Spatial information locates cell signaling while delivering comprehensive gene expression profiling within tissues. 

![Signaling between adjacent cells. The Notch protein functions as a receptor for ligands that activate or inhibit such receptors. Receptor-ligand interactions ground cell signaling and communication, often requiring close proximity between cells. ](https://upload.wikimedia.org/wikipedia/commons/0/04/Notchccr.svg){alt='alt text for accessibility purposes'}

<a href="https://commons.wikimedia.org/wiki/File:Notchccr.svg">Fred the Oyster</a> Public domain, via Wikimedia Commons <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="license">CC BY-SA 4.0 DEED</a>

Spatial transcriptomics addresses a key obstacle in single-cell studies: the 
loss of spatial information when tissues are dissociated to release viable,
intact cells. Spatial organization and structure determine function in most
tissues and organs, so capturing spatial information is critical for 
understanding tissue function in neuroscience, cancer biology, developmental 
biology, and most other fields. Immunology is the only field that doesn't depend 
on spatial information to determine structure and function because immune cells 
often circulate in liquid media like blood and lymph. 

## Spatial transcriptomics technologies
Spatial transcriptomics technologies broadly fall within two groups: 
imaging-based and sequencing-based methods. Imaging-based technologies read 
transcriptomes in situ using microscopy. They identify mRNA species through
hybridization with flourescent gene-specific probes as in fluorescence in 
situ hybridization (FISH). Current FISH methods employ multiple hybridization
rounds, with risk of error for each transcript growing exponentially with each
round. These methods are limited in the size of tissue that they can profile and
can be time-consuming due to microscopic imaging. 

In situ sequencing amplifies and sequences mRNAs directly within a block 
or section of tissue. Probes profile one or two bases at a time using different
fluorophores, eventually revealing the identity of the mRNA through imaging. In 
situ sequencing can accommodate larger tissue sections than can FISH, though 
FISH methods are more efficient at detecting mRNA. Like FISH, in situ sequencing requires considerable imaging time on a microscope. 

Sequencing-based methods capture, sequence, and use next-generation sequencing 
to count mRNA in situ while retaining positional information. This is distinct
from in situ sequencing because next-generation sequencing is employed, not 
imaging. Sequencing-based methods retain spatial information through 
microdissection and microfluidics, or through ligation of mRNAs to arrays of 
barcoded probes. Microdissection techniques process tissue sections for 
transcriptomic profiling. The are useful for profiling transcriptomes as a first 
pass, but have limited spatial resolution and require many samples for 
sequencing. Array-based methods capture mRNA with spatially-barcoded probes and
sequence them, or alternatively print the array onto tissue using microfluidic
channels to produce a unique barcode at every position in the tissue. 

In this lesson we will use data from positionally barcoded arrays. The probes are not gene-specific; rather, they record position.
Array-based methods use tissue frozen below the temperature at which RNA 
degrades, though some methods can use formalin-fixed tissues embedded in 
paraffin. 

## 10X Visium technology

![A sequencing-based spatial transcriptomics method using printed spots on a slide. ](https://upload.wikimedia.org/wikipedia/commons/1/14/Spatial_transcriptomics_ii.png){alt='alt text for
accessibility purposes'}
<a href="https://commons.wikimedia.org/wiki/User:Jasquatch">James Chell</a>, <a href="https://commons.wikimedia.org/wiki/File:Spatial_transcriptomics_ii. png">Spatial transcriptomics ii</a>, <a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" rel="license">CC BY-SA 4.0</a>

:::::: keypoints
 - Spatial transcriptomics provides the location of individual cells relative to neighboring cells and cell structures.
 - A cell's location is useful data for describing its phenotype, state, and cell and tissue function.
 - Spatial transcriptomics addresses a key obstacle in single-cell studies: the loss of spatial information when tissues are dissociated.
::::::
