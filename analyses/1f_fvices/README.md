# f-VICE calculations using BMO predictions

## Usage
### Run BMO and generate V-plots on predicted bound motifs
```sh
snakemake -npr --configfile config.yml -q
```

### Aggregate results
```sh
snakemake -npr --configfile config.yml all_fvices -q
```