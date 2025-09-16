# 📦 blwgsR

## BootLeg Whole Genome Sequencing R package #TODO change this name maybe, something that incorprates ATAC and RNAseq, pacbio, and 10X
*an R package for parsing and combining pacbio long read sequencing with 10X multiomic data*

*R*
---

### 📦 R Installation

```R
devtools::install_github("ratoncito/blwgsR", force = TRUE)
library(blwgsR)

#parse the read into inserts barcodes sampleIDs 
read_info <- parse_reads(blwgsR::Fitzwalter2025, out.file = "test.csv")

#count CAG repeats and demultiplex the reads into cells 
cell_info <- demultiplex_reads(reads = read_info, out.file = "None")


![chromium_barcode_info.png]
