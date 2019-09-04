# TF-binding algorithms comparisons

## Usage
### Downsample BAM files, call peaks, and run TF-binding prediction algorithms (minus PIQ)
```sh
snakemake -npr --configfile config.yml -q
```
### Get AUC-PR and F1 results
This will generate a `sh.roc` and `sh.f1` file containing all the jobs that need to be submitted (one per line). You can submit the using a job scheduler (*e.g.* `qsub sh.f1`) or run them in serial (*e.g.* `sh sh.roc`).

**Important**: For unknown reasons, this snakemake script does not generate the two job files reliably and they end up empty. If that happens, you need to run the other two commands in the code block below until a non-empty output is generated for each (it may take *several* tries).
```sh
snakemake -npr --configfile config.yml all_evaluate -q

# Regenerate ROC and PR job file
snakemake --configfile config.yml all_evaluate -R get_roc
# Regenerate F1 job file
snakemake --configfile config.yml all_evaluate -R get_f1
```

### Aggregate results
```sh
snakemake -npr --configfile config.yml all_aggregate -q
```

### Compare BMO and PIQ
```sh
snakemake -npr --configfile config.yml all_piq -q
```