# The information landscape of human transcription factor and DNA interactions
_A companion repository_

## Description
This repository will allow you to generate all the main figures from our manuscript.

## Usage
Each folder inside `figures` points to a main figure and its related supplementary figures. Follow the instructions inside each one to perform the analyses.

## Requirements
All requirements (except for Slurm) are included in the conda virtual environment provided here (`environment.yml`). We highly recommend you to use it. If that is not an option, below are the core requirements:
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (5.5.0+)
* [Slurm](https://slurm.schedmd.com) (16.05.09+) or some other job scheduler. The cluster config files here are set up for slurm.
* [Rstudio](https://www.rstudio.com) (1.1.456+)
* [R](https://www.r-project.org) (3.5.1+). Note that each notebook will have different package requirements.
* [Python 3](https://www.python.org)

## Contact
* For questions directly about the paper, the corresponding author is Dr. Stephen C.J. Parker (scjp@umich.edu).
* If anything here is (or feels) broken or you just have questions regarding how to set up the analyses, contact me directly at albanus@umich.edu.
