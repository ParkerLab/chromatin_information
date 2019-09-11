# mat2jaspar.R
# Converts pwms into jaspar format

## Input
# ID (not used)
# 455 46 146 353
# 174 25 759 42
# 0 941 0 59
# 5 0 988 7

## Output
# > ID
# A [455 174 0 5 ]
# C [46 25 941 0 ]
# G [146 759 0 988 ]
# T [353 42 59 7 ]

options(stringsAsFactors = F)

args <- commandArgs(T)
if(length(args) ==2){
    infile <- args[1]
    id     <- args[2]
} else(stop("Must provide input file and motif ID!"))

# infile = "/lab/data/motifs/pwm/ENCODE2013/hg18_matrixes/CTCF_known2.mat"
# id = "CTCF_known2"

mat <- read.table(infile, skip = 1, col.names = c("A", "C", "G", "T"))
mat <- t(mat)

# Convert to JASPAR format
write(paste0(">", id), file = stdout())
for(i in 1:4){
    dat <- paste(mat[i,], '', collapse = '')
    out <- paste0(rownames(mat)[i], " [", dat, "]", collapse = '')
    write(out, file = stdout())
}
