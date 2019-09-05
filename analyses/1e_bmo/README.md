# TF-binding algorithms comparisons

## Usage
### Downsample BAM files, call peaks, and run TF-binding prediction algorithms (minus PIQ)
```sh
snakemake -npr --configfile config.yml -q
```
### Get AUC-PR and F1 results
This will generate a `sh.roc` and `sh.f1` file containing all the jobs that need to be submitted (one per line). You can submit the using a job scheduler (*e.g.* `qsub sh.f1`) or run them in serial (*e.g.* `sh sh.roc`).

**Important**: For unknown reasons, this snakemake script does not generate the two job files reliably and one (or both) can end up empty. If that happens, you need to execute the respective commands in the code block below until a non-empty output is generated for each (it may take *several* tries).
```sh
snakemake -npr --configfile config.yml all_evaluate -q

# Check if the sh.roc and sh.f1 files were created
ls -hlrt
# Regenerate sh.roc
snakemake --configfile config.yml all_evaluate -R get_roc
# Regenerate sh.f1
snakemake --configfile config.yml all_evaluate -R get_f1
```

### Aggregate results
```sh
snakemake -npr --configfile config.yml all_aggregate -q
```

### Compare BMO and PIQ
```sh
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q all_metrics
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q aggregate_f1 aggregate_prauc
```