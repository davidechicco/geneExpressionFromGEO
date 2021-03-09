

#' Function that reads in the GEO code of a dataset, and returns the gene expression dataframe.
#' 
#' @param datasetGeoCode the GEO code of a dataset.
#' @param retrieveGeneSymbols a boolean flag stating if the function should retrieve the gene symbols or not.
#' @param verbose a boolean flag stating if helping messages should be printed or not
#' @return a gene expression dataset.
#' @examples
#' geneExpressionDF1 <- getGeneExpressionFromGEO("GSE3268", FALSE, FALSE)
getGeneExpressionFromGEO <- function(datasetGeoCode, retrieveGeneSymbols, verbose = FALSE) 
{

            GSE_code <- datasetGeoCode
            
            # check   URL
            checked_html_text <- "EMPTY_STRING"
            checked_html_text <- xml2::read_html("https://ftp.ncbi.nlm.nih.gov/geo/series/")
            
            if(all(checked_html_text == "EMPTY_STRING")) {
         
                    cat("The web url https://ftp.ncbi.nlm.nih.gov/geo/series/ is unavailable right now. Please try again later. The function will stop here\n")
                    return(NULL)
         
            } else {

            gset <- GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
            
            thisGEOplatform <- toString((gset)[[1]]@annotation)

            if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
            gset <- gset[[idx]]
            gene_expression <- as.data.frame(Biobase::exprs(gset))


            if(retrieveGeneSymbols == TRUE) {
                gene_expression$GeneSymbol <- ""


                # we retrieve the platform details
                platform_ann <- annotate::readGEOAnn(GEOAccNum = thisGEOplatform)
                platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
                
                if (verbose == TRUE) {
                    print("sort(names(platform_ann_df))")
                    print(sort(names(platform_ann_df)))
                }
                                
                # "gene_assignment
                platformsWithGene_assignmentField <- c("GPL11532", "GPL23126", "GPL6244")
                
                # "Gene Symbol"
                platformsWithGeneSpaceSymbolField <- c("GPL80", "GPL8300", "GPL80", "GPL96", "GPL570", "GPL571")
                
                # "gene_symbol"
                platformsWithGene_SymbolField <- c("GPL20115")
                
                # "symbol"
                platformsWithSymbolField <- c("GPL1293", "GPL6102", "GPL6104", "GPL6883", "GPL6884")
                 
                # "GENE_SYMBOL
                platformsWith_GENE_SYMBOL_Field <- c("GPL13497", "GPL14550", "GPL17077", "GPL6480")

                # if symbol
                if(thisGEOplatform %in% c(platformsWithGeneSpaceSymbolField, platformsWithGene_SymbolField, platformsWithSymbolField, platformsWith_GENE_SYMBOL_Field)   ) {
                
                        emptyGeneSymbol <- ""
                        FIRST_GENE_EXPRESSION_INDEX <- 2

                        if (verbose == TRUE)    cat("\n[start] loop for the association of the gene symbols to the probeset ID's: completed \n", sep="")
                        # we start from 2 because 1 is the label
                        for(k in FIRST_GENE_EXPRESSION_INDEX:nrow(gene_expression)) {

                        currentCompletionPerc <- k*100 / nrow(gene_expression)
                        kForPrint <- 1000
                        if (verbose == TRUE) if ((k %% kForPrint)==0) { cat(dec_two(currentCompletionPerc), "% ", sep="") }

                        if(thisGEOplatform %in% platformsWithGeneSpaceSymbolField) thisSymbol <- platform_ann_df[platform_ann_df$ID==rownames(gene_expression[k,]),]$"Gene Symbol"
                        if(thisGEOplatform %in% platformsWithGene_SymbolField) thisSymbol <- platform_ann_df[platform_ann_df$ID==rownames(gene_expression[k,]),]$"gene_symbol"
                        if(thisGEOplatform %in% platformsWithSymbolField) thisSymbol <- platform_ann_df[platform_ann_df$ID==rownames(gene_expression[k,]),]$"symbol"
                        if(thisGEOplatform %in% platformsWith_GENE_SYMBOL_Field) thisSymbol <- platform_ann_df[platform_ann_df$ID==rownames(gene_expression[k,]),]$"GENE_SYMBOL"
                        
                        gene_expression[k,]$GeneSymbol <- thisSymbol
                        
                        }
                        if (verbose == TRUE) cat("\n [end] loop for the association of the gene symbols to the probeset ID's \n ", sep="")

                }  else if(thisGEOplatform %in% platformsWithGene_assignmentField) {            # if assignment  
                
                
                FIRST_GENE_EXPRESSION_INDEX <- 2
                SIZE_SPLIT_STRING <- 2
                GENE_SYMBOL_INDEX <- 2

               if (verbose == TRUE)  cat("\n[start] loop for the association of the gene symbols to the probeset ID's: completed ", sep="")

                # we start from 2 because 1 is the label
                for(k in FIRST_GENE_EXPRESSION_INDEX:nrow(gene_expression)) {
                    
                    currentCompletionPerc <- k*100 / nrow(gene_expression)
                    kForPrint <- 1000
                    if (verbose == TRUE) if ((k %% kForPrint)==0) { cat(dec_two(currentCompletionPerc), "%  ", sep="") }

                    thisAssignment <- NULL
                    thisProbesetID <- NULL
                    thisProbesetID <- rownames(gene_expression[k,])
                    # cat("this probeset ID: ", thisProbesetID, "\t", sep="" )
                    thisAssignment <- platform_ann_df[platform_ann_df$ID==thisProbesetID, ]$"gene_assignment"

                    split_string <- strsplit(thisAssignment, "//")

                    if (length(split_string[[1]]) >= SIZE_SPLIT_STRING) {
                        thisGeneSymbol_temp <- NULL
                        thisGeneSymbol  <- NULL
                        thisGeneSymbol_temp <- split_string[[1]][GENE_SYMBOL_INDEX]
                        thisGeneSymbol <- gsub(" ", "", thisGeneSymbol_temp, fixed = TRUE)        
                    #  cat("\t this gene symbol: ", thisGeneSymbol, "\n", sep="")

                    gene_expression[k,]$GeneSymbol <- thisGeneSymbol

                    } else {

                        gene_expression[k,]$GeneSymbol <- thisAssignment

                    }
                }
                    if (verbose == TRUE) cat("\n [end] loop for the association of the gene symbols to the probeset ID's ", sep="")                
                } else {
                
                    if (verbose == TRUE) cat("\n\n[Impossible to perform gene symbol retrieval]\n")
                    if (verbose == TRUE) cat("We're sorry but the indicated platform (", thisGEOplatform, ") is not among the platforms included in this function.\nThe gene symbol retrieval cannot be performed.\nPlease visit the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", thisGEOplatform, " website for more information about this platform.\n\n", sep="")
                    gene_expression$GeneSymbol <- NULL
                
                }
            
            }
            
            return(gene_expression)
        }
}   
