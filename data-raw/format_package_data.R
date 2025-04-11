#load in important metadata 
#TODO how to have this within the R package??
#TODO need to figure out how to attach as a r data object 
#these data can be found within a download of the cellranger software from 10
#specifically here -> cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz
#and here -> cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
Barcode10X_WhiteList <- read.csv("data-raw/atac_737K-arc-v1.txt.gz", header = FALSE)$V1
Multiome_key = data.frame(rna = read.csv("data-raw/rna_737K-arc-v1.txt.gz", header = FALSE)$V1,
                          atac = read.csv("data-raw/atac_737K-arc-v1.txt.gz", header = FALSE)$V1)

#TODO where did I get this from?
#the 10X website I believe 
SampleIndex_WhiteList <- read.csv("data-raw/Single_Index_Kit_N_Set_A (2).csv", header = FALSE)

#format metadata
#TODO preformat maybe so they are automatically in the R package and we can save some compute 
colnames(SampleIndex_WhiteList) <- c("sample", "ref1", "ref2", "ref3", "ref4")
SampleIndex_WhiteList$sample <- gsub("SI-NA-", "", SampleIndex_WhiteList$sample)
SampleIndex_WhiteList <- melt(SampleIndex_WhiteList, id = "sample")
colnames(SampleIndex_WhiteList) <- c("sample", "index_type", "index")

#save all of of necessary data for inclusion in the published R package
use_data(SampleIndex_WhiteList)
use_data(Barcode10X_WhiteList)
use_data(Multiome_key)
