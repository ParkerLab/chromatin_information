# The information landscape of human transcription factor and DNA interactions
_A companion repository_

Ricardo D’Oliveira Albanus, Yasuhiro Kyono, John Hensley, Arushi Varshney, Peter Orchard, Jacob O. Kitzman, Stephen C. J. Parker

## To-do
- [ ] Auxiliary scripts
    - [ ] conda environment
    - [ ] install_dependencies.R
    - [ ] download_fastq.sh
    - [ ] download_processed_files.sh

## Description
This repository will allow you to generate all the main figures from our manuscript.

## Usage
### Browsing processed results
The `notebooks_processed` folder contains a Rmd (R Markdown) file for each set of analyses. Just open these files in RStudio and knit the html output.
### Generating results from scratch
1. Go to each of the `analyses` subfolders **in order** and run the commands described in the README.
2. Go to the `notebooks_local` folder and execute the Rmd files. They should be able to knit the html output once step 1 is complete.

## Requirements
### Using our conda environment
Assuming that you have a Linux 64-bit system, download and install Anaconda 3:
```
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh
```
Create the environment:
```
conda env create --name chromatin_information --file environment.yaml
source activate chromatin_information
```

### Manually installing dependencies
All requirements (except for Slurm) are included in the conda virtual environment provided here (`environment.yml`). We highly recommend you to use it. If that is not an option, below are the core requirements and the versions we used:
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (5.5.0)
* [Slurm](https://slurm.schedmd.com) (16.05.09) or some other job scheduler. The cluster config files here are set up for slurm.
* [Rstudio](https://www.rstudio.com) (1.1.456). Required for the RMarkdown notebooks
* [R](https://www.r-project.org) (3.5.1). Note that each notebook will have different package requirements.
* [Python](https://www.python.org) (3.5.3)
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) (2.63)

## Contact
* For questions directly about the paper, the corresponding author is Dr. Stephen C.J. Parker (scjp@umich.edu).
* If anything here is (or feels) broken or you just have questions regarding how to set up the analyses, contact Ricardo Albanus directly at albanus@umich.edu.
