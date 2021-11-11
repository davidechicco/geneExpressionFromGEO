setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

## Script that installs the R packages requested by the getGeneExpressionFromGEO() function

cat(":: Script that installs the R packages ::\n:: required by the GeneExpressionFromGEO R package ::\n\n")

# Here we install the CRAN missing packages
list.of.packages <- c("easypackages", "xml2", "markdown", "knitr", "rmarkdown", "pacman", "R.utils") # other packages
new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://utstat.toronto.edu/cran/")

# Here we install the Bioconductor missing packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "Biobase", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
if(length(bioCpackagesNotInstalled)) cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

cat("\n:: Additional information about the geneExpressionFromGEO R package ::\n:: can be found on the https://github.com/davidechicco/geneExpressionFromGEO website ::\n\n")
