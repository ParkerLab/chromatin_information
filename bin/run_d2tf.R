library(dnase2tf)
args <- commandArgs(TRUE)

# Input data

mapfiledir <- '/lab/data/reference/human/hg19/mappability'
assemseqdir <- '/lab/data/reference/human/hg19'

# User input
datafilepath          <- args[1]
outputfilepath        <- args[2]
paired                <- args[3]
hotspotfilename       <- args[4]
dftfilename_for_human <- args[5]

# Run code
dnase2tf(datafilepath, hotspotfilename, mapfiledir, outputfilepath,
    assemseqdir = assemseqdir,
    dftfilename = dftfilename_for_human,
    maxw = 30,          # Maximum motif width
    minw = 6,           # Minumum motif width
    z_threshold = -2,   # Z-score threshold
    numworker = 1,      # Number of cores
    paired = paired)    # Paired-end or not
