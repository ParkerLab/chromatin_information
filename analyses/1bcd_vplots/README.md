# Make ChIP-seq V-plots and calculate f-VICEs

## Usage
Make V-plots and calculate f-VICEs for the TFs with ChIP-seq data. f-VICEs correspond to the ones reported in Figure 1D.
```sh
snakemake --configfile config.yml
```
Make V-plots downsampling to same number of motifs and fragments, used for Figures 1B and 1C.
```sh
snakemake -s rules/same_n.smk --configfile rules/same_n.yml
```