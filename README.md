# The information landscape of human transcription factor and DNA interactions
_A companion repository_

Ricardo D’Oliveira Albanus, Yasuhiro Kyono, John Hensley, Arushi Varshney, Peter Orchard, Jacob O. Kitzman, Stephen C. J. Parker

## Description
This repository will allow you to generate all the main figures from our manuscript.

## Usage
### Browsing processed results
The `notebooks_processed` folder contains a Rmd (R Markdown) file for each set of analyses. Just open these files in RStudio and knit the html output.
### Generating results from scratch
1. Go to each of the `analyses` subfolders **in order** and run the commands described in the README.
2. Go to the `notebooks_local` folder and execute the Rmd files. They should be able to knit the html output once step 1 is complete.

## Software requirements
### Using our conda environment
Assuming that you have a Linux 64-bit system, download and install Anaconda 3:
```
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh
```
Create the base environment:
```sh
conda env create --file env/base.yml
source activate chromatin_information
```
This will provide most of the command line software (python, R, samtools, *etc*.)

### Install other software
Fetch the necessary R packages
```sh
Rscript scripts/install_dependencies.R
```
Other required software can be found following the links below. The versions we used are in parentheses.
* [HINT](http://www.regulatory-genomics.org/hint) (RGT 0.12.1)
* [DNase2TF](https://sourceforge.net/projects/dnase2tfr) (1.0)
* [PIQ](http://piq.csail.mit.edu) (1.3)
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) (2.63)

### Manually installing core dependencies
* [R](https://www.r-project.org) (3.5.1)
* [Python](https://www.python.org) (3.5.3)
* [Rstudio](https://www.rstudio.com) (1.1.456). Required for the RMarkdown notebooks. 
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (5.5.0)
* [Slurm](https://slurm.schedmd.com) (16.05.09) or some other job scheduler. The cluster config files here are set up for slurm.


## Contact
* For questions directly about the paper, the corresponding author is Dr. Stephen C.J. Parker (scjp@umich.edu).
* If anything here is (or feels) broken or you just have questions regarding how to set up the analyses, contact Ricardo Albanus directly at albanus@umich.edu.
