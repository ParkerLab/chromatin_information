cranURL <- "https://cran.r-project.org/src/contrib/Archive"
packages <- c(
    "data.table/data.table_1.12",
    "ggplot2_3.1.0",
    "scales/scales_0.5.0",
    "dplyr/dplyr_0.7.8",
    "tidyr/tidyr_0.8.2",
    "RColorBrewer/RColorBrewer_1.1-2",
    "PRROC/PRROC_1.3",
    "ROCR/ROCR_1.0-7",
    "MASS/MASS_7.3-50"
)

for(package in packages) {
    packageURL <- sprintf("%s/%s.tar.gz", cranURL, package)
    install.packages(packageURL, repos=NULL, type='source')   
}

# NGSplot
install.packages("doMC", dep=T)
install.packages("caTools", dep=T)
install.packages("utils", dep=T)

source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("Rsamtools")
biocLite("ShortRead")

# Knitr and Rmd notebooks
install.packages(DT)
install.packages(knitr)