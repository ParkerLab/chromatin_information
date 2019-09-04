# Make ChIP-seq V-plots and calculate f-VICEs

## Usage
Make V-plots and calculate f-VICEs for the TFs with ChIP-seq data
```sh
snakemake --configfile config.yml
```
Make V-plots downsampling to same number of motifs and fragments
```sh
snakemake -s rules/same_n.smk --configfile rules/same_n.yml
```