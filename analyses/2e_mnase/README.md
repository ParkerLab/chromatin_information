# MNase read mapping and motif analyses

## Usage
### Before running
Adjust the necessary paths in the `config.yml` file. The most critical are the BWA index (see `0_mapping` for more details) and the output directory for the sra tools `prefetch`, which defaults to your home directory. This may not be ideal if your home directory is on a separate partition.

### Excute pipeline
```sh
snakemake -npr --configfile config.yml -q
```