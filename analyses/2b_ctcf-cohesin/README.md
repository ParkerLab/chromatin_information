# CTCF-cohesin +/- regions and sonicated ATAC-seq data

## Usage
### Generate V-plots for CTCF-cohesin +/- regions
```sh
snakemake -npr --configfile config.yml -q
```

### Generate V-plots directly comparing the sonicated and non-sonicated sample
```sh
sh sonicated.sh
```