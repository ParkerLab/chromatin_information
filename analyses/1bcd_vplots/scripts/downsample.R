options(stringsAsFactors = F)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

seed <- 17071234

args <- commandArgs(T)
if (length(args) != 4) {
  stop("Must provide 4 arguments: infile, n_frags, n_motifs, output")
}

# setwd("/home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN")
# args[1] <- "work/vplots/lab_gm12878/AP1__AP1_known5__ENCFF784PEF-JUNB.vsignal.filtered.gz"
# args[2] <- "work/min_size.lab_gm12878.txt"
# args[3] <- "work/min_motifs.lab_gm12878.txt"
# args[4] <- "work/vplots/lab_gm12878/AP1__AP1_known5__ENCFF784PEF-JUNB.downsampled.vsignal.gz"

infile <- args[1]
frags <- args[2]
nmotifs <- args[3]
outfile <- args[4]

d <- data.table::fread(infile, data.table = F)
frags_n <- read.table(frags)[,1]
motif_n <- read.table(nmotifs)[,1]

idx_rank <- d %>% 
  group_by(FeatureNumber) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  mutate(rank = 1:n())

set.seed(seed)
left_join(d, idx_rank, by = "FeatureNumber") %>% 
  arrange(rank) %>% 
  filter(rank <= motif_n) %>%
  sample_n(frags_n) %>% 
  select(-n, -rank, -rank) %>% 
  write.table(gzfile(outfile), col.names = T, row.names = F, 
              quote = F, sep = "\t")

