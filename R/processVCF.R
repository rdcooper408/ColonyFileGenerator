
#' Load and Process a VCF
#'
#' This function loads in a VCf using the vcfR::read.vcfR. The vcf file is then
#' processed to match the Colony input format. This function returns an updated
#' VCF genotype dataframe with modified genotype notation. This dataframe is the
#' required input for the main function "writeColonyFile()".
#'
#' @param vcf.file Path to the input VCF file
#' @return A dataframe of the modified VCF genotype
#' @export
processVCF <- function(vcf.file){
    vcf <- vcfR::read.vcfR(vcf.file) # user vcfR::read.vcfR() to load in a filtered VCF file
    cat("Extracting Genotype DataFrame\n")
    vcf.df <- as.data.frame(vcf@gt) # extract genotype dataframe from vcf object
    vcf.df[,-1] <- as.character(gsub(pattern = "(^./.).*",replacement = "\\1" ,x = as.matrix(vcf.df[,-1]))) # removes the info fields after the GT data
    vcf.df[,-1] <- as.character(gsub(pattern = "(^.)\\|(.).*",replacement = "\\1/\\2" ,x = as.matrix(vcf.df[,-1]))) # removes info field after GT data for phased loci
    vcf.df<- vcf.df[,-1] # remove info column
    vcf.fix <- as.data.frame(vcfR::getFIX(vcf)) # extract locus information
    vcf.fix$chrom.pos <- paste(vcf.fix$CHROM, vcf.fix$POS, sep = ".") # create a chrom.position column (unique for each locus)
    row.names(vcf.df) <- vcf.fix$chrom.pos

    # Convert gt notation: Colony uses allelic counts (0,1,2,3, etc..) instead of the typical vcf format
    # Such that: "." => 0, 0 => 1, 1 => 2
    cat("Modifying Genotype Notation")
    vcf.df[,] <- as.character(gsub(pattern = "0/1",replacement = "1/2" ,x = as.matrix(vcf.df)))
    vcf.df[,] <- as.character(gsub(pattern = "1/1",replacement = "2/2" ,x = as.matrix(vcf.df)))
    vcf.df[,] <- as.character(gsub(pattern = "0/0",replacement = "1/1" ,x = as.matrix(vcf.df)))
    vcf.df[,] <- as.character(gsub(pattern = "\\./\\.",replacement = "0/0" ,x = as.matrix(vcf.df)))
    return(vcf.df)
}


#' @export
sumNonMiss <- function(x){sum(x != "0/0")}

#' Filter a Colony VCF df for sample missingness
#'
#' This function takes the input from the processVCF() function and filters for
#' sample missingness. So individuals that are missing more than perc.missing
#' are dropped from the output.
#'
#' @param vcf.df modified VCF dataframe from processVCF() function
#' @param perc.missing the allowable % (0-100) of loci for which a sample can be missing data
#' @return A VCF dataframe with some samples removed
#' @export
filterMissVCF <- function(vcf.df, perc.missing = 50){
    non0.df <- apply(X = vcf.df, 2, FUN = sumNonMiss)
    tot <- nrow(vcf.df)
    cutoff <- round(tot*((100-perc.missing)/100),0)
    n.drop <- sum(non0.df < cutoff)
    drop.id <- colnames(vcf.df[,non0.df < cutoff])
    cat(paste("Dropping ", n.drop, " samples with more than ", perc.missing,"% missing data\n",sep = "" ))
    keep <- which(non0.df >= cutoff)
    vcf.cut <- vcf.df[,keep]
    cat("Samples dropped:\n")
    print(drop.id)
    return(vcf.cut)
}


#' Reduce the length of column and row names
#'
#' This function cuts the column and row names of a dataframe.
#'
#' @param vcf.df any dataframe
#' @param num.col.char number of desired characters in column names
#' @param num.row.char number of desired characters in row names
#' @return A VCF dataframe with shortened column and row names
#' @export
shortNames <- function(vcf.df, num.row.char = 20, num.col.char = 20){
    row.names(vcf.df) <- gsub(pattern = paste("^.*(.{",num.row.char = 20,"}$)", sep = ""), replacement = "\\1", rownames(vcf.df)) # Take only the LAST 20 characters of loci names (max allowed in colony)
    colnames(vcf.df) <- gsub(pattern = paste("(^.{",num.col.char = 20,"}).*$", sep = ""), replacement = "\\1", colnames(vcf.df)) # Take only the FIRST 20 characters of sample names (max allowed in colony)
}


