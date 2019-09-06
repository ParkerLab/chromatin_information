options(stringsAsFactors = F)

args <- commandArgs(T)

print(args)

if(length(args) < 3){
    stop("User should provide 3-4 arguments: {input} {output} {FDR} {fraction}")
}
infile <- args[1]
outfile <- args[2]

if (args[3] == "nominal"){
  fdr_thresh <- args[3]
} else {
  fdr_thresh <- as.numeric(args[3]) 
}
if(is.na(args[4])) {
    min_fraction <- 0.1
} else {
    min_fraction <- as.numeric(args[4])
}

d <- read.table(infile, header = T, sep = "\t")
print(dim(d))

# Keep significant regions
if (is.numeric(fdr_thresh)) {
  d$padj <- p.adjust(d$p, method = "BH")
  d <- d[d$padj <= fdr_thresh,]  
} else {
  d$padj <- d$p
  d <- d[d$padj <= 0.05,]  
}


print(dim(d))

# Keep regions where there's at least x% of each allele
d <- d[abs(d$fraction_ref - .5) <= .5 - min_fraction,]

print(dim(d))

# Write output
write.table(d, outfile, sep = "\t", col.names = T, row.names = F, quote = F)