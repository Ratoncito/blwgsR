# 📦 blwgsR

## BootLeg Whole Genome Sequencing R package #TODO change this name maybe, something that incorprates ATAC and RNAseq, pacbio, and 10X
*an R package for parsing and combining pacbio long read sequencing with 10X multiomic data*

*R*
---

### 📦 R Installation

```R
devtools::install_github("ratoncito/blwgsR")


### test run the function

#we provide a example dataset from the Fitzwalter 2025 paper, this is a subset of 100,000 reads from this paper

#this function has many optional inputs, for looking at the HTT_exon1 gene these options work well, for other genes this will require different inputs, for investigating all genes see the _ function (in progress)
read_info <- parse_reads()

#now that we have the reads we can count the CAG repeats like so

read_info$cag_repeats <- count_cags()

#now we can demultiplex the reads into cells
cell_info <- demultiplex()
