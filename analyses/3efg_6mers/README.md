# Perform 6-mer analyses
## Usage
### Calculate 6-mer f-VICEs
```sh
snakemake -npr --configfile config.yml -q
```

### Collect data
```sh
snakemake -npr --configfile config.yml -q all_fvices
```

### Run allelic imbalance f-VICE calculations
```sh
snakemake -npr --configfile config.yml -q all_snp2kmer
```