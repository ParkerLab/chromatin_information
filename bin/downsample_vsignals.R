options(stringsAsFactors = F)
setwd("/lab/work/albanus/2018_cohesin_and_tfs")

library(dplyr)
library(tidyr)

# args <- commandArgs(T)
# exp <- args[1]
# coh <- args[2]
# outdir <- args[3]

exp <- "CTCF__CTCF_known2__ENCFF963PJY-CTCF"
coh <- "ENCFF002CPK"
outdir <- "work/vplots_sonicated"

read_vsignal <- function(sample, exp, coh, type) {
  handle <- sprintf("%s__%s.%s.vsignal.gz", exp, coh, type)
  f <- file.path("work/vsignal", sample, handle)
  read.table(f, header = T) %>% 
    filter(FragmentSize > 40) %>% 
    # Sort by number of reads per motif instance
    group_by(FeatureNumber) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    arrange(-n) %>% 
    # Recalculate indices using sorted motif instances
    mutate(FeatureNumber = !duplicated(FeatureNumber)) %>% 
    mutate(FeatureNumber = cumsum(FeatureNumber))
}

sub_p <- read_vsignal("lab_gm12878_subsampled_3", exp, coh, "plus")
sub_m <- read_vsignal("lab_gm12878_subsampled_3", exp, coh, "minus")

son_p <- read_vsignal("lab_gm12878_sonicated_3", exp, coh, "plus")
son_m <- read_vsignal("lab_gm12878_sonicated_3", exp, coh, "minus")

# Downsample by fragments (prioritizing low signal motifs)
downsample_frags <- function(vsignal, max_frags, inv = T) {
  if(inv) {
    ww <- 1 / vsignal$n
  } else {
    ww = vsignal$n
  }
  vsignal %>% 
    sample_n(max_frags, weight = ww) %>% 
    # Sort by number of reads per motif instance
    group_by(FeatureNumber) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    arrange(-n)
}

nfrags <- c(nrow(sub_m), nrow(sub_p), nrow(son_m), nrow(son_p)) %>% 
  min()

d_sub_p <- downsample_frags(sub_p, nfrags)
d_sub_m <- downsample_frags(sub_m, nfrags)

d_son_p <- downsample_frags(son_p, nfrags)
d_son_m <- downsample_frags(son_m, nfrags)

# Downsample by motifs
downsample_motifs <- function(vsignal, max_motifs) {
  vsignal %>% 
    filter(FeatureNumber <= max_motifs)
}

nmotifs <- c(
  length(unique(d_sub_m$FeatureNumber)),
  length(unique(d_sub_p$FeatureNumber)),
  length(unique(d_son_m$FeatureNumber)),
  length(unique(d_son_p$FeatureNumber))
) %>% 
  min()

d_sub_p <- downsample_motifs(d_sub_p, nmotifs)
d_sub_m <- downsample_motifs(d_sub_m, nmotifs)

d_son_p <- downsample_motifs(d_son_p, nmotifs)
d_son_m <- downsample_motifs(d_son_m, nmotifs)


nfrags <- c(nrow(d_sub_m), nrow(d_sub_p), nrow(d_son_m), nrow(d_son_p)) %>% 
  min()

d_sub_p <- downsample_frags(d_sub_p, nfrags) %>% select(-n)
d_sub_m <- downsample_frags(d_sub_m, nfrags) %>% select(-n)

d_son_p <- downsample_frags(d_son_p, nfrags) %>% select(-n)
d_son_m <- downsample_frags(d_son_m, nfrags) %>% select(-n)


c(nrow(d_sub_m), nrow(d_sub_p), nrow(d_son_m), nrow(d_son_p))

c(
  length(unique(d_son_p$FeatureNumber)), length(unique(d_son_p$FeatureNumber)),
  length(unique(d_son_p$FeatureNumber)), length(unique(d_son_p$FeatureNumber))
)

# Output data
write_out <- function(vsignal_sub, sample, exp, coh, type) {
  outdir_with_sample <- file.path(outdir, sample)
  dir.create(outdir_with_sample, showWarnings = F, recursive = T)
  
  handle <- sprintf("%s__%s.%s.vsignal.gz", exp, coh, type)
  f <- file.path(outdir_with_sample, handle)
  write.table(vsignal_sub, f, quote = F, sep = "\t", 
              col.names = T, row.names = F)
}

write_out(d_sub_p, "lab_gm12878_subsampled_3", exp, coh, "plus")
write_out(d_sub_m, "lab_gm12878_subsampled_3", exp, coh, "minus")
write_out(d_son_p, "lab_gm12878_sonicated_3", exp, coh, "plus")
write_out(d_son_m, "lab_gm12878_sonicated_3", exp, coh, "minus")
