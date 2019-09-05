#####
# get_ROC_all.R
#   This script will take a list of named files and ChIP-seq associations, and generate one overlayed
#  ROC plot for comparing two methods. In addition, I'm including an extra output with all the reported
#  AUCs plus the size of the ChIP-seq experiments for downstream filtering purposes.

# Usage: Rscript scripts/get_ROC_hint.R file scoreCol motif truePositives binNames output

# 1) file:          Bed file with scores to be ranked
# 2) scoreCol:      Which column the scores start
# 3) motif:         Name of the motif. Used for plot title
# 4) truePositives: Bed file containing a subset of coordinates to be used as true positives
# 5) binNames:      One-lined, tab-separated file with the names of the categories present in the input file
# 6) outputPath:    Path of the output files (no '/' at the end)

# Author: Ricardo D'O. Albanus
#####


options(stringsAsFactors = F)
suppressWarnings(library(ROCR))
suppressWarnings(library(PRROC))
suppressWarnings(library(RColorBrewer))
args<-commandArgs(TRUE)

if(length(args) != 6)
{
  stop("\n\nWarning: Wrong number of arguments!\nPlease supply: file scoreCol motif truePositives binNames output.\n\n")
}


## Step 1: Define variables and load files

infile   <- args[1]
colstart <- as.numeric(args[2])
motif    <- args[3]
chipfile <- args[4]
binfile  <- args[5]
outpath  <- args[6]

bins <- as.character(read.table(file = binfile, sep ="\t"))
n <- length(bins)
columns <- (colstart -1) + (1:n)
outfile <- paste(outpath, motif, sep ="/")

dat <- read.table(file = infile)
chip <- read.table(file = chipfile)


## Step 2: Define prediction function

makePred <- function(x, TP, column){
  # Make comparison table
  data <- data.frame(ref = coords, 
                     score = dat[,column])  # dat[,column] is the posteriors
  data$tp <- as.numeric(data$ref %in% TP)   # 0/1 vector to see if coordinate is in TP

  # Make performance objects and get area under curve (AUC)
  pred <- prediction( predictions = data$score, labels = data$tp)         
  perf <- performance(pred, measure = "tpr", x.measure = "fpr") # Syntax for making ROC curves (as in ?performance) 
  auc <- round(as.numeric(performance(pred, "auc")@y.values),4) # http://stackoverflow.com/questions/4903092/calculate-auc-in-r
  pauc <- round(as.numeric(performance(pred, "auc", fpr.stop = 0.05)@y.values)/0.05, 4) # Calculate partial AUC at .05 threshold
  # Return list with performance object and AUC
  RP.perf <- performance(pred, "prec", "rec")
  # Calculate PR curve AUC with PRROC
  pr_auc <- pr.curve(scores.class0 = data$score[data$tp == 1 & !is.na(data$score)],
                     scores.class1 = data$score[data$tp == 0 & !is.na(data$score)],
                     sorted = F, curve = F)
  return(list(perf = perf, auc = auc, pauc = pauc, pr = RP.perf, prauc = pr_auc$auc.integral))
}


## Step 3: Run analysis

# Make true positives (TP) vector:
# concatenate coords and remove whitespace for matching
TP <- gsub(" ", "", paste(chip[,1], chip[,2], chip[,3], sep = "_"), fixed = T) 
coords <- gsub(" ", "", paste(dat[,1], dat[,2], dat[,3], sep = "_"), fixed = T)

# Make performance objects using the make_pred function
for (i in 1:n){
  assign(paste('perf', i, sep = ''), makePred(dat, TP, columns[i]))
}; rm(i)
# Prepare and write output file
auc_all <- NULL
pauc_all <- NULL
prauc_all <- NULL
for(i in 1:n){
  pauc_all[i] <- get(paste('perf', i, sep=''))$pauc
  auc_all[i] <- get(paste('perf', i, sep=''))$auc
  prauc_all[i] <- get(paste('perf', i, sep=''))$prauc
}; rm(i)

output1 <- rbind(bins, auc_all)
output1 <- cbind(output1, c('ChIPseq_size', length(TP)))
write.table(output1, file = paste(outfile, '.1.0.out', sep = ''), quote = F, row.names = F, sep = '\t', col.names = F)
output2 <- rbind(bins, pauc_all)
output2 <- cbind(output2, c('ChIPseq_size', length(TP)))
write.table(output2, file = paste(outfile, '.0.1.out', sep = ''), quote = F, row.names = F, sep = '\t', col.names = F)
output3 <- rbind(bins, prauc_all)
output3 <- cbind(output3, c('ChIPseq_size', length(TP)))
write.table(output3, file = paste(outfile, '.prauc.out', sep = ''), quote = F, row.names = F, sep = '\t', col.names = F)

## Step 4: Plot

if(n <= 12){
  cols <- brewer.pal(8, 'Dark2')
} else{
  cols <- rainbow(n)
}

orderedROC <- order(auc_all, decreasing = T)
orderedPR <- order(prauc_all, decreasing = T)

# Plot ROC curves
cairo_pdf(file = paste(outfile, '.pdf', sep = ''), height=6, width=6)
for(i in 1:n){
    plot(get(paste('perf',i, sep =''))$perf, col = cols[i], lwd = 2, add = i>=2, 
         xlim = c(0,1), ylim = c(0,1), xaxs='i', yaxs='i', main = motif)
}
abline(a = 0, b = 1, lwd = 1.5, lty = 2)
vecLegend <- NULL
for(i in 1:n){
  vecLegend[i] <- paste(bins[orderedROC[i]], paste('(AUC = ', round(auc_all[orderedROC[i]], 3), ')', sep = ''),
                        sep = '\t')
}
legend(x = 'bottomright', vecLegend, col = cols[orderedROC], pch = 15, y.intersp = 1, bg = 'white')
dev.off()

# Plot Precision-Recall curves
cairo_pdf(file = paste(outfile, '_pr.pdf', sep = ''), height=6, width=6)
for(i in 1:n){
    plot(get(paste('perf',i, sep =''))$pr, col = cols[i], lwd = 2, add = i>=2, 
         xlim = c(0,1), ylim = c(0,1), xaxs='i', yaxs='i', main = motif)
}
vecLegend <- NULL
for(i in 1:n){
  vecLegend[i] <- paste(bins[orderedPR[i]], paste('(prAUC = ', round(prauc_all[orderedPR[i]], 3),')',sep = ''),
                        sep = '\t')
}
legend(x = 'topright', vecLegend, col = cols[orderedPR], pch = 15, y.intersp = 1, bg = 'white')
dev.off()


## Step 5: Manual troubleshooting

# myplot <- function(n){
#   d <- dat[,c(columns[n],17)]
#   d <- d[order(d[,1], decreasing = T),]
#   plot(d$tp, type = 'l', main = bins[n])
# }
# 
# for(i in 1:n){
#   myplot(i)
#   rm(i)
# }