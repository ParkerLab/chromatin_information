# The information landscape of human transcription factor and DNA interactions
_A companion repository_

Ricardo D’Oliveira Albanus, Yasuhiro Kyono, John Hensley, Arushi Varshney, Peter Orchard, Jacob O. Kitzman, Stephen C. J. Parker

## Description
This repository will allow you to generate all the main figures from our [manuscript](https://www.biorxiv.org/content/10.1101/777532v2).

## Usage
### Browsing processed results
1. Download the full processed data from zenodo and place it at the root of the cloned repository:
```sh
wget https://zenodo.org/record/3478583/files/chromatin_information_manuscript_data.tar.gz
tar -xzvf chromatin_information_manuscript_data.tar.gz
```
2. The `notebooks` folder contains a Rmd (R Markdown) file for each set of analyses. Just open these files in RStudio and knit the html output.

### Generating results from scratch
Go to each of the `analyses` subfolders **in order** and run the commands described in the README.

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
* [CENTIPEDE](http://centipede.uchicago.edu) (1.2)
* [DNase2TF](https://sourceforge.net/projects/dnase2tfr) (1.0)
* [HINT](http://www.regulatory-genomics.org/hint) (0.12.1, RGT 1.1.1)
* [PIQ](http://piq.csail.mit.edu) (1.3)
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) (2.63)

We bundled the core BMO scripts here, but if you wish to use it in your own projects, please go to the [official BMO repository](https://github.com/ParkerLab/BMO).

### Manually installing the core dependencies included in conda environment
In case you cannot use the environment, below are the core dependencies:
* [R](https://www.r-project.org) (3.5.1)
* [Python](https://www.python.org) (3.5.3)
* [Rstudio](https://www.rstudio.com) (1.1.456). Required for the RMarkdown notebooks. 
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) (2.26.0)
* [samtools](http://samtools.sourceforge.net) (1.9)
* [MACS2](https://github.com/taoliu/MACS) (2.1.1.20160309)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (5.5.0)
* [Slurm](https://slurm.schedmd.com) (16.05.09) or some other job scheduler. The cluster config files here are set up for slurm.


## Contact
* For questions directly about the paper, the corresponding author is Dr. Stephen C.J. Parker (scjp@umich.edu).
* If anything here is (or feels) broken or you just have questions regarding how to set up the analyses, contact Ricardo Albanus directly at albanus@umich.edu.
