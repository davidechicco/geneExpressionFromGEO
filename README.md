# geneExpressionFromGEO #

`geneExpressionFromGEO`: an easy method to retrieve gene expression and its symbol from GEO accession code

## Summary ##

The `getGeneExpressionFromGEO()` function is an easy method that reads in the Gene Expression Omnibus (GEO) code of a gene expression dataset, retrieves its data from GEO, (optional) retrieves the gene symbols of the dataset, and returns a simple dataframe table containing all the data. Platforms available: GPL11532, GPL23126, GPL6244, GPL80, GPL8300, GPL80, GPL96, GPL570, GPL571, GPL20115, GPL1293,  GPL6102, GPL6104, GPL6883, GPL6884, GPL13497, GPL14550, GPL17077, GPL6480. 

The GEO datasets are downloaded from the URL <https://ftp.ncbi.nlm.nih.gov/geo/series/>.
This function has been designed for beginners and users having limited experience with Bioconductor.

## Installation ##

To run `getGeneExpressionFromGEO`, you need to have the following programs and packages installed in your computer:

* R (version > 3.1)
* R Bioconductor packages `Biobase, annotate, GEOquery`

You can install the `geneExpressionFromGEO` package and its dependencies from CRAN, and load it, with the following commands typed in the `R` terminal console:

    R
    install.packages("geneExpressionFromGEO", repos='http://cran.us.r-project.org')
    library("geneExpressionFromGEO")
    
If it is impossible to download the package and its dependencies from CRAN, you can download the `geneExpressionFromGEO` package from this GitHub repository, and then can execute the `install_packages.r` script that will install all the dependencies automatically:

    cd geneExpressionFromGEO
    R
    source("install_packages.r")
    
Afterwards,  you execute the `geneExpressionFromGEO.r` file from an R terminal:

    source("geneExpressionFromGEO.r")

## Execution instructions ##

To run `getGeneExpressionFromGEO`, the only parameter need to have is the GEO accession code of the dataset you want to download.
The other two parameters are boolean. The first one allow you to decide if you want `getGeneExpressionFromGEO` to retrieve all the gene symbols of the probesets of the dataset, and assign them to the probesets.
The last parameter allows you to decide if you want the function to print messages during its operations or not.

## An example ##

We want to retrieve the gene expression dataset with GEO accession [GSE3268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3268)  of the platform [GPL96](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96). Here are the commands we can use in a R shell environment:
    
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO("GSE3268",  associateSymbolsToGenes, verbose)
    
## Article
Additional information about this project will be available in the following peer-reviewed published article:

> Davide Chicco, "geneExpressionFromGEO: an R package to facilitate data reading from Gene Expression Omnibus (GEO)". [Microarray Data Analysis, Methods in Molecular Biology, volume 2401](https://doi.org/10.1007/978-1-0716-1839-4), Springer Protocols, New York City, New York, USA, 2021. In press.
    
## Contacts ##

The `geneExpressionFromGEO` package was developed by [Davide Chicco](https://www.DavideChicco.it). Questions should be
addressed to davidechicco(AT)davidechicco.it
