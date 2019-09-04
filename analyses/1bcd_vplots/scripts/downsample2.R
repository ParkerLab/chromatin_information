options(stringsAsFactors = F)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

seed <- 17071234

args <- commandArgs(T)
if (length(args) != 3) {
  stop("Must provide 3 arguments: infile, n_motifs, output")
}

# setwd("/home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN")
# args[1] <- "work/vplots/lab_gm12878/AP1__AP1_known5__ENCFF784PEF-JUNB.vsignal.filtered.gz"
# args[2] <- "work/min_motifs.lab_gm12878.txt"
# args[3] <- "work/vplots_same_n_motif/lab_gm12878/AP1__AP1_known5__ENCFF784PEF-JUNB.downsampled.vsignal.gz"

infile <- args[1]
nmotifs <- args[2]
outfile <- args[3]

d <- data.table::fread(infile, data.table = F)
motif_n <- read.table(nmotifs)[,1]

set.seed(seed)
idx_rank <- d %>% 
  group_by(FeatureNumber) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  sample_n(length(.$FeatureNumber)) %>% 
  mutate(rank = 1:n())

left_join(d, idx_rank, by = "FeatureNumber") %>% 
  arrange(rank) %>% 
  filter(rank <= motif_n) %>%
  select(-n, -rank, -rank) %>% 
  write.table(gzfile(outfile), col.names = T, row.names = F, 
              quote = F, sep = "\t")

