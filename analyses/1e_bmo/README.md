# TF-binding algorithms comparisons

## Usage (GM12878)
### Downsample BAM files, call peaks, and run TF-binding prediction algorithms
```sh
snakemake -npr --configfile config.yml -q
```
### Get AUC-PR and F1 results
This will generate a `sh.roc` and `sh.f1` file containing all the jobs that need to be submitted (one per line). You can submit the using a job scheduler (*e.g.* `qsub sh.f1`) or run them in serial (*e.g.* `sh sh.roc`).

**Important**: For unknown reasons, this snakemake script does not generate the two job files reliably and one (or both) can end up empty. If that happens, you need to execute the respective commands in the code block below until a non-empty output is generated for each (it may take *several* tries).
```sh
snakemake -npr --configfile config.yml all_evaluate -q

# Check if the sh.roc and sh.f1 files are not empty
wc -l sh.*
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
# Run PIQ
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q
# Compare methods
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q all_metrics
# Aggregate results
snakemake -npr -s rules/piq_bmo.smk --configfile config.yml -q aggregate_f1 aggregate_prauc
```

## Usage (HepG2)
### Downsample BAM files, call peaks, and run TF-binding prediction algorithms
```sh
snakemake -npr -s rules_hepg2/Snakefile --configfile rules_hepg2/config.yml -q
```
### Get AUC-PR and F1 results
This will generate a `sh.roc` and `sh.f1` file containing all the jobs that need to be submitted (one per line). You can submit the using a job scheduler (*e.g.* `qsub sh.f1`) or run them in serial (*e.g.* `sh sh.roc`).

**Important**: For unknown reasons, this snakemake script does not generate the two job files reliably and one (or both) can end up empty. If that happens, you need to execute the respective commands in the code block below until a non-empty output is generated for each (it may take *several* tries).
```sh
snakemake -npr -s rules_hepg2/Snakefile --configfile rules_hepg2/config.yml all_evaluate -q

# Check if the sh.roc and sh.f1 files are not empty
wc -l sh.*
# Regenerate sh.roc
snakemake -s rules_hepg2/Snakefile --configfile rules_hepg2/config.yml all_evaluate -R get_roc
# Regenerate sh.f1
snakemake -s rules_hepg2/Snakefile --configfile rules_hepg2/config.yml all_evaluate -R get_f1
```

### Aggregate results
```sh
snakemake -npr -s rules_hepg2/Snakefile --configfile rules_hepg2/config.yml all_aggregate -q
```

### Compare BMO and PIQ
```sh
# Run PIQ
snakemake -npr -s rules_hepg2/piq_bmo.smk --configfile rules_hepg/config.yml -q
# Compare methods
snakemake -npr -s rules_hepg/piq_bmo.smk --configfile rules_hepg/config.yml -q all_metrics
# Aggregate results
snakemake -npr -s rules_hepg/piq_bmo.smk --configfile rules_hepg/config.yml -q aggregate_f1 aggregate_prauc
```