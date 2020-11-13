# getGeneExpressionFromGEO #

`getGeneExpressionFromGEO`: an easy method to retrieve gene expression and its symbol from GEO accession and platform ID

## Summary ##

The `getGeneExpressionFromGEO()`is an easy function that reads in the Gene Expression Omnibus (GEO) code and the platform ID of a gene expression dataset, retrieves its data from GEO, (optional) retrieves the gene symbols of the dataset, and returns a simple dataframe table containing all the data. Platforms available: GPL11532, GPL23126, GPL6244, GPL80, GPL8300, GPL80, GPL96, GPL570, GPL571, GPL20115, GPL1293,  GPL6102, GPL6104, GPL6883, GPL6884, GPL13497, GPL14550, GPL17077, GPL6480. 

The GEO datasets are downloaded from the URL <https://ftp.ncbi.nlm.nih.gov/geo/series/>.
This function has been designed for beginners and users having limited experience with Bioconductor.

## Installation ##

To run `getGeneExpressionFromGEO`, you need to have the following programs and packages installed in your machine:

* R (version > 3.1)
* R Bioconductor packages **Biobase, annotate, GEOquery** with the following commands:

`if (!requireNamespace("BiocManager", quietly = TRUE))`
`install.packages("BiocManager")`
`listOfBiocPackages <- c("Biobase", "annotate", "GEOquery")`

`bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )`
`cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")`

`# check there's still something left to install`
`if( length(bioCpackagesNotInstalled) ) {`
   `BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])`
`}`

`library("easypackages")`
`libraries(listOfBiocPackages)`

After installing R and these functions, you can install the `getGeneExpressionFromGEO`, package by cloning this repository.

## Execution instructions ##

To run `getGeneExpressionFromGEO`, you just need to have the GEO accession code of the dataset you want to download, and the ID of its platform. 
The other two parameters are boolean. The first one allow you to decide if you want `getGeneExpressionFromGEO` to retrieve all the gene symbols of the probesets of the dataset, and assign them to the probesets.
The last parameter allows you to decide if you want the function to print messages during its operations or not.

# An example

We want to retrieve the gene expression dataset with GEO accession GSE3268  of the platform GPL96:

    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO("GSE3268", "GPL96", associateSymbolsToGenes, verbose)
    
## Contacts ##

The `getGeneExpressionFromGEO` package was developed by [Davide Chicco](https://www.DavideChicco.it). Questions should be
addressed to davidechicco(AT)davidechicco.it
