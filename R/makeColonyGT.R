# This function converts a VCD data.frame into a colony-style genotype table
## This colony gt splits each sample/locus into two columns, one for each allele

#' @export
makeColonyGT <- function(vcf.input){
    print(paste("There are", ncol(vcf.input), "samples in your VCF"))
    print(paste("There are", nrow(vcf.input), "Loci in your VCF"))

    a1 <- a2 <- vcf.input
    rownames(a1) <- a1.names <- final.colnames <- (rownames(vcf.input))
    rownames(a2) <- a2.names <- paste(rownames(vcf.input), "-2", sep = "")

    a1[,] <- as.character(gsub(pattern = "(^\\d)/\\d",replacement = "\\1" ,x = as.matrix(a1)))
    a2[,] <- as.character(gsub(pattern = "^\\d/(\\d)",replacement = "\\1" ,x = as.matrix(a2)))

    aa <- rbind(a1, a2)
    aa <- aa[sort(rownames(aa)), ]
    gt <- t(aa)
}
