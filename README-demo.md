# ELASTomics analysis
Code from Shiomi et al. (2024) - "High-throughput mechanical phenotyping and transcriptomics of single cells."


## System requirements
All code has been tested on R version 4.2.2.


## Installation guide
Rstudio is downloaded from https://www.rstudio.com/.

The libraries required for the analysis are listed at the beginning of the R scripts.

The typical install time is within 30 min.


## Demo

###  Data analysis
To download the data analysis code, git clone using the following command:

    git clone https://github.com/RIKEN-Microfluidics-Lab/ELASTomics.git


### Scripts

DEMO＿ELASTomics.R outputs figures corresconding to Figs. 3b-g, and Supplementary Fig. 12 in the manuscript of Shiomi et al. (2024).

Expected outputs are commented as corresponding figure numbers at the end of each analysis code.

Execution of all analyses takes less than a few days.
SCENIC analysis takes the most time.

The code is released under the GNU Public License (GPL 3.0).

### 10x format data

10x format to be loaded in Seurat is available from https://riken-share.ent.box.com/s/r7noivs3qn33y1etalbr4lx9v0mltga8.
The directories correspond to the respective samples as follows:


- TIG1-C: 10x data of none-electroporated TIG-1 cells (PDL = 42, 59)

- TIG1-C: 10x data of electroporated TIG-1 cells (PDL = 42, 59)


## Raw sequence data
 All demultiplexed sequencing data have been deposited on the Sequencing Read Archive and are available for download under accession PRJNA841462.
