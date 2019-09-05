options(stringsAsFactors = F)
library(optparse)

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="input file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output folder", metavar="character"),
    make_option(c("-c", "--colnum"), type="integer", default=NULL, 
                help="starting column", metavar="integer"),
    make_option(c("-n", "--name"), type="character", default=NULL, 
                help="motif name", metavar="character"),
    make_option(c("-b", "--bins"), type="character", default=NULL, 
                help="bins file", metavar="character"),
    make_option(c("-t", "--tp"), type="character", default=NULL, 
                help="true positives file", metavar="character"),
    make_option(c("-p", "--peaks"), type="character", default=NULL, 
                help="ATAC-seq peaks file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
args <- parse_args(opt_parser)

## Step 1: Define variables and load files

infile   <- args$file
colstart <- args$colnum
motif    <- args$name
chipfile <- args$tp
peakfile <- args$peaks 
binfile  <- args$bins
outpath  <- args$output
fdrfile  <- args$fdr

# # Test
# infile   <- 'work/tf-binding/concatenated_files/buenrostro_rep1/CTCF_known2.80.bed'
# colstart <- 7
# motif    <- 'CTCF_known2.80'
# chipfile <- file.path("/lab/work/albanus/gm12878/motif_intersects_encode2017",
#                       "intersections/CTCF__CTCF_known2__ENCFF096AKZ-CTCF.bed")
# peakfile <- 'work/tf-binding/motifs_in_atac-seq_peaks/buenrostro_rep1/CTCF_known2.80.bed'
# binfile  <- 'work/tf-binding/bin_names.txt'
# outpath  <- 'f1_gm12878'

# Read files
bins <- as.character(read.table(file = binfile, sep ="\t"))
n <- length(bins)
columns <- (colstart -1) + (1:n)

dir.create(outpath, showWarnings = F)
handle <- motif
outfile <- paste0(outpath, "/", motif, ".bed")

d <- read.table(file = infile)
tp <- read.table(file = chipfile)
peaks <- read.table(file = peakfile)
peaks <- peaks[peaks$V7 == 1,]

# Get binarized ATAC-seq peaks and TP calls
d$tp <- paste(d$V1,d$V2,d$V3,sep="_") %in% paste(tp$V1,tp$V2,tp$V3,sep="_")
d$peaks <- paste(d$V1,d$V2,d$V3,sep="_") %in% paste(peaks$V1,peaks$V2,peaks$V3,sep="_")


# Define thresholds
thresholds <- c(
    -log10(0.05),  # BMO
    .99,  # CENTIPEDE
    -log10(0.05),  # DNase2TF
    0,  # HINT
    0.99, # sbCENTIPEDE
    0, 0 # ATAC-peaks (twice)
)

ntp <- sum(d$tp)
ntn <- nrow(d) - ntp

methods <- c(bins, 'WithinPeaks')
metrics <- list(f1 = NULL, precision = NULL, recall = NULL)
for(i in 1:n){
    d <- d[order(d[,columns[i]], decreasing = T),]
    
    d$index <- 1:nrow(d)
    d$recall <- cumsum(d$tp)/ntp
    # Calculate Precision for each of the positions (TP / (TP + FP))
    d$ntp <- cumsum(d$tp)
    d$nfp <- cumsum(d$tp == F)
    d$precision <- d$ntp/(d$nfp + d$ntp)
    
    if(thresholds[i] != 0){
        idxFdr <- sum(d[,columns[i]] >= thresholds[i])
    } else {
        idxFdr <- sum(d[,columns[i]] > 0)
    }
    
    precision = d$precision[idxFdr]
    recall = d$recall[idxFdr]
    
    metrics$precision[i] <- precision
    metrics$recall[i] <- recall
    metrics$f1[i] <- (precision * recall) / (precision + recall) * 2
}

# Add within peaks F1
d <- d[order(d$peaks, decreasing = T),]

d$index <- 1:nrow(d)
d$recall <- cumsum(d$tp)/ntp
# Calculate Precision for each of the positions (TP / (TP + FP))
d$ntp <- cumsum(d$tp)
d$nfp <- cumsum(d$tp == F)
d$precision <- d$ntp/(d$nfp + d$ntp)
idxFdr <- sum(d$peaks)

precision = d$precision[idxFdr]
recall = d$recall[idxFdr]

metrics$precision[n+1] <- precision
metrics$recall[n+1] <- recall
metrics$f1[n+1] <- (precision * recall) / (precision + recall) * 2

out <- data.frame(precision = metrics$precision,
                  recall = metrics$recall,
                  f1 = metrics$f1,
                  method = c('BMO', 'CENTIPEDE', "DNase2TF", "HINT", 
                             "CENTIPEDE_sb", "Peaks_log10", 'WithinPeaks'),
                  experiment = handle,
                  total_true_positives = ntp)

write.table(out, file = paste(outfile, '.out', sep =""), quote = F, 
            row.names = F, col.names = F, sep = "\t")
