# This program takes the allelic imbalance files and split it into two groups 
# of files, REF and ALT. Each k-mer is written to its own separate file in the
# provided output directory.
# Author: Ricardo D'O. Albanus

options(stringsAsFactors = F)

library(dplyr)
library(tidyr)

# Parse commands
args <- commandArgs(T)
infile <- args[1]
outdir <- args[2]

# # Test
# infile <- file.path("work/kmers/6-mers/allelic_imbalance",
#                     "buenrostro_rep1.significant.kmers.txt")
# outdir <- "test_split"


# Helper functions
write_output <- function(split_element, allele){
    stopifnot(allele %in% c("REF", "ALT"))
    if(allele == "REF"){
        kmer <- unique(split_element$ref_kmer)
    } else{
        kmer <- unique(split_element$alt_kmer)
    }
    if(length(kmer) != 1){
        stop("Must be given unique a k-mer!")
    }
    outhandle <- sprintf("%s__%s__%s.bed", sample, allele, kmer)
    outfile <- file.path(outdir, outhandle)
    write.table(split_element, outfile, sep = "\t", quote = F, col.names = F, 
                row.names = F)
}


# Make sure we have the output directory
if(!dir.exists(outdir)){
    dir.create(outdir)
}

# Read file and extract sample name
d <- read.table(infile, header = T, sep = "\t")
sample <- gsub(".*\\/(.*)\\.signif.*", "\\1", infile)

# Split by REF k-mer
d_ref <- split(d, d$ref_kmer)
lapply(d_ref, write_output, "REF")

# Split by ALT k-mer
d_alt <- split(d, d$alt_kmer)
lapply(d_alt, write_output, "ALT")
