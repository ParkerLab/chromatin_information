# The information landscape of human transcription factor and DNA interactions
_A companion repository_

## To-do
- [ ] Processed data
    - [ ] Motif scans
    - [ ] Motif f-VICEs
    - [ ] 6-mer f-VICEs
    - [ ] Asymmetry
    - [ ] Protein domains
    - [ ] FRAP times
- [ ] Auxiliary scripts
    - [ ] download_fastq.sh
- [ ] Analyses and notebooks
    - [x] 0_mapping
    - [x] 1bcd_vplots
    - [x] 1e_bmo
        - [x] PIQ
    - [ ] 1f_fvices
    - [ ] 2b_ctcf-cohesin
    - [ ] 2f_asymmetry
    - [ ] 3c_enrichments
    - [ ] 3efg_6mers

## Description
This repository will allow you to generate all the main figures from our manuscript.

## Usage
Each folder inside `figures` points to a main figure and its related supplementary figures. Follow the instructions inside each one to perform the analyses.

## Requirements
### Using our conda environment
Assuming that you have a 64-bit system, on Linux, download and install Anaconda 3:
```sh
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh
```
Then, create the environment:
```sh
conda env create --name chromatin_information --file environment.yaml
source activate chromatin_information
```

### Manually installing dependencies
All requirements (except for Slurm) are included in the conda virtual environment provided here (`environment.yml`). We highly recommend you to use it. If that is not an option, below are the core requirements:
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (5.5.0+)
* [Slurm](https://slurm.schedmd.com) (16.05.09+) or some other job scheduler. The cluster config files here are set up for slurm.
* [Rstudio](https://www.rstudio.com) (1.1.456+). Required for the RMarkdown notebooks
* [R](https://www.r-project.org) (3.5.1+). Note that each notebook will have different package requirements.
* [Python](https://www.python.org) (3.5.3+)

## Contact
* For questions directly about the paper, the corresponding author is Dr. Stephen C.J. Parker (scjp@umich.edu).
* If anything here is (or feels) broken or you just have questions regarding how to set up the analyses, contact me directly at albanus@umich.edu.
