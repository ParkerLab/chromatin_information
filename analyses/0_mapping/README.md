# ATAC-seq datasets mapping to hg19
## Before running

1. Make sure that `../../data/fastq` is populated with the necessary fastq files for outlined in the `config.yml` file and that they are formatted **exactly** as `{accession}.{1,2}.fastq`. Note that the file **must not be gzipped**.
2. Make a [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml) and place/link it at `../../data/bwa_index/hg19`.

## Usage

```{sh}
snakemake [options] --configfile config.yml --cluster-config cluster.yml
```