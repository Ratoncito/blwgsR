#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
demultiplex_reads <- function(){
  info <- read.csv("results/all_data_htt_hits_info3.csv")
  message("\n\n\n================================#| STEP 5: DEMULTIPLEXING |#================================")
  # Create a function for getting the mode
  getmode <- function(v) {
    uniqv <- BiocGenerics::unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  #merge cell information with sample key
  info <- merge(info, SampleIndex_WhiteList[,c(1,3)], by.x = "sampleid", by.y = "index")
  #make a column that indicates a unique cell
  info$barcode_sample <- paste(info$barcode, info$sample, sep = "_")
  message(paste("unique cells recovered: ", length(BiocGenerics::unique(info$barcode_sampleid))))
  
  message("cell per sample recovered: ")
  message(BiocGenerics::table(cell_info$sample))
  
  
  cell_info <- parallel::mclapply(BiocGenerics::unique(info$barcode_sample), function(cell){
    sub <- info[info$barcode_sampleid == cell,]
    
    res <- data.frame(atac = sub$barcode[1], 
                      cag_repeats_mean = mean(sub$cag_repeats), 
                      cag_repeats_median = median(sub$cag_repeats), 
                      cag_repeats_mode = getmode(sub$cag_repeats), 
                      cag_repeats_var = var(sub$cag_repeats),
                      cag_repeats_sd = sd(sub$cag_repeats),
                      number_of_reads = nrow(sub),
                      sample = sub$sample[1], 
                      barcode_sample = cell)
    
    return(res)
  }, mc.cores = detectCores() - 1)
  cell_info <- do.call(rbind, cell_info)
  
  #
  #message(length(BiocGenerics::unique(cell_info$barcode_sample)))
  
  #merge cell information with the atac to rna key
  cell_info <- merge(cell_info, Multiome_key, by = "atac")
  
  #check for errors 
  #TODO how do I address this?
  colSums(is.na(cell_info))  # Count missing values per column
  
  message(paste("Number of unique cells: ", nrow(cell_info)))
  message(paste("Min number of reads per cell: ", min(cell_info$number_of_reads)))
  message(paste("Max number of reads per cell: ", max(cell_info$number_of_reads)))
  message(paste("Mean number of reads per cell: ", mean(cell_info$number_of_reads)))
  
  write.csv(cell_info, file = "results/all_data_htt_hits_cell_info3.csv")
  message("cell level info saved")
  
  if(plot){
    ggplot2::ggplot(data = cell_info, aes(x = cag_repeats_mode, y = cag_repeats_sd)) + geom_point()
    
    ggplot2::ggplot(data = cell_info, aes(x = log10(cag_repeats_mode), y = log10(cag_repeats_var + 1))) + 
      geom_point() + 
      geom_smooth(method = "lm", color = "red") + 
      geom_smooth() + 
      theme_bw()
  }
  
}