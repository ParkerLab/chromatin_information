# ATAC-seq datasets mapping to hg19
## Before running

1. Download the necessary fastq files using the `download_fastq.sh` script in the `scripts` folder in the root directory of this repository. It will populate the `../../data` directory with the input files.
2. Make a [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml) and place/link it at `../../data/bwa_index/hg19` relative to this directory.

## Usage

```sh
snakemake --configfile config.yml --cluster-config cluster.yml
```