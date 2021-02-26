
#' Write a Colony Input Data File
#'
#' This function writes a data file for use in the software Colony (https://www.zsl.org/science/software/colony).
#' This function requires a modified vcf genotype dataframe that can be generated using the processVCF function.
#' Note: Colony does not accept row or column names longer than 20 characters. This function includes "num.char"
#' argument to cut these names.
#'
#'
#' @param colony.run.name Name for this colony run (included in the filename)
#' @param offspring.id List of all potential offspring sample names
#' @param father.id List of all potential father sample names
#' @param mother.id List of all potential mother sample names
#' @param vcf.input.df Modified genotype dataframe generated from the processVCF function
#' @param out.path The full path to the desired output folder, if it does not exist it will be created
#' @param num.char The desired length to cut rownames and colnames (must be <= 20 for Colony)
#' @param out.file.suffix The file ending that will be appended to the colony.run.name
#' @param n.runs  Number of replicate Colony runs
#' @param l.run  Length of Colony run. Give a value of 1, 2, 3, 4 to indicate short, medium, long, very long run.
#' @param seed <- 1234 The random seed for Colony
#' @param updateAF  0/1=Not updating/updating allele frequency
#' @param MD  2/1=Dioecious/Monoecious species
#' @param InBD  0/1=No inbreeding/inbreeding
#' @param ploid  0/1=Diploid species/HaploDiploid species
#' @param PM  0/1=Polygamy/Monogamy for males & females
#' @param clone  0/1=Clone inference =No/Yes
#' @param SibSizeScale  0/1=Full sibship size scaling =No/Yes
#' @param SibSizePrior  0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size
#' @param popAF  0/1=Unknown/Known population allele frequency
#' @param monitor  0/1=Monitor method by Iterate. Give a value of either 0 or 1 to indicate monitoring the intermediate results by iterate number or running time. Always choose value 0 for run without Windows GUI.
#' @param monInterval  Monitor interval in Iterate
#' @param AnalysisType  Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)
#' @param Precision  1/2/3=low/medium/high Precision for Full likelihood
#' @param MarkerTypes  Marker Type: 0 co-dominant / 1 dominant, use "0@" or "1@" to repeat the option for all loci
#' @param AlleleDropRate  Allelic Drop-out rate, use "x@" to repeat the option for all loci
#' @param OtherErrorRate  Other Error rate, use "x@" to repeat the option for all loci
#' @param ProbParentIncl  Probability parent included
#' @return Writes a Colony input data file into the specified directory
#' @export
writeColonyFile <- function(colony.run.name="colony_run1",
                            offspring.id=offspring.id,
                            father.id=father.id,
                            mother.id=mother.id,
                            vcf.input.df=vcf.input.df,
                            out.path="./",
                            out.file.suffix="colonyInputFile.txt",
                            num.char = 20,
                            seed=1234,
                            updateAF=0,
                            MD=2,
                            InBD=0,
                            ploid=0,
                            PM=0,
                            clone=0,
                            SibSizeScale=0,
                            SibSizePrior=0,
                            popAF=0,
                            n.runs=1,
                            l.run=1,
                            monitor=0,
                            monInterval=10000,
                            AnalysisType=1,
                            Precision=3,
                            MarkerTypes="0@",
                            AlleleDropRate="0.02@",
                            OtherErrorRate="0.02@",
                            ProbParentIncl=0.01
){


    # Cut the sample names to match the Colony Format (less than 20 characters)
    if(num.char > 20){stop("Colony requires that locus names and samples names be less than 20 characters. Please set num.char <= 20")}
    row.names(vcf.input.df) <- gsub(pattern = paste("^.*(.{",num.char,"}$)", sep = ""), replacement = "\\1", rownames(vcf.input.df)) # Take only the LAST 20 characters of loci names (max allowed in colony)
    colnames(vcf.input.df) <- gsub(pattern = paste("(^.{",num.char,"}).*$", sep = ""), replacement = "\\1", colnames(vcf.input.df)) # Take only the FIRST 20 characters of sample names (max allowed in colony)
    father.id.cut <- gsub(pattern = paste("(^.{",num.char,"}).*$", sep = ""), replacement = "\\1", father.id)
    mother.id.cut <- gsub(pattern = paste("(^.{",num.char,"}).*$", sep = ""), replacement = "\\1", mother.id)
    offspring.id.cut <- gsub(pattern = paste("(^.{",num.char,"}).*$", sep = ""), replacement = "\\1", offspring.id)

    # Set up genotype tables
    father.index <- which(colnames(vcf.input.df) %in% father.id.cut)
    mother.index <- which(colnames(vcf.input.df) %in% mother.id.cut)
    offspring.index <- which(colnames(vcf.input.df) %in% offspring.id.cut)

    n.father=length(father.index)
    n.mother=length(mother.index)
    n.offspring=length(offspring.index)
    n.loci=nrow(vcf.input.df)

    # Build the Colony Specific genotype tables (requires makeColonyGT() function)
    cat("Generating Father Genotype Table\n")
    father.gt <- makeColonyGT(vcf.input.df[,father.index])
    cat("\n\nGenerating Mother Genotype Table\n")
    mother.gt <- makeColonyGT(vcf.input.df[,mother.index])
    cat("\n\nGenerating Offspring Genotype Table\n")
    offspring.gt <- makeColonyGT(vcf.input.df[,offspring.index])

    # Write the output file
    out.file.name <- paste(colony.run.name, out.file.suffix, sep = ".")
    dir.create(out.path, showWarnings = FALSE) # create the output directory if it does not exist
    out.file <- paste(out.path,out.file.name, sep = "")
    cat(paste("\n\nWriting file to: ", out.file, "\n", sep = ""))

    cat(c(colony.run.name, "! Dataset name","\n"), file = out.file, sep = "\t", append = F)
    cat(c(out.file.name, "! Output file name","\n"), file = out.file, sep = "\t", append = T)
    cat(c(n.offspring, "! Number of offspring in the sample","\n"), file = out.file, sep = "\t", append = T)
    cat(c(n.loci, "! Number of Loci","\n"), file = out.file, sep = "\t", append = T)
    cat(c(seed, "! Seed for random number generator","\n"), file = out.file, sep = "\t", append = T)
    cat(c(updateAF, "! 0/1=Not updating/updating allele frequency","\n"), file = out.file, sep = "\t", append = T)
    cat(c(MD, "! 2/1=Dioecious/Monoecious species","\n"), file = out.file, sep = "\t", append = T)
    cat(c(InBD, "! 0/1=No inbreeding/inbreeding","\n"), file = out.file, sep = "\t", append = T)
    cat(c(ploid, "! 0/1=Diploid species/HaploDiploid species","\n"), file = out.file, sep = "\t", append = T)
    cat(c(PM, PM, "! 0/1=Polygamy/Monogamy for males & females","\n"), file = out.file, sep = "\t", append = T)
    cat(c(clone, "! 0/1=Clone inference =No/Yes","\n"), file = out.file, sep = "\t", append = T)
    cat(c(SibSizeScale, "! 0/1=Full sibship size scaling =No/Yes","\n"), file = out.file, sep = "\t", append = T)
    cat(c(SibSizePrior, "! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size","\n"), file = out.file, sep = "\t", append = T)
    cat(c(popAF, "! 0/1=Unknown/Known population allele frequency","\n"), file = out.file, sep = "\t", append = T)
    cat(c(n.runs, "! Number of runs","\n"), file = out.file, sep = "\t", append = T)
    cat(c(l.run, "! Length of run","\n"), file = out.file, sep = "\t", append = T)
    cat(c(monitor, "! 0/1=Monitor method by Iterate","\n"), file = out.file, sep = "\t", append = T)
    cat(c(monInterval, "! Monitor interval in Iterate","\n"), file = out.file, sep = "\t", append = T)
    cat(c("0", "! non-Windows version","\n"), file = out.file, sep = "\t", append = T)
    cat(c(AnalysisType, "! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)","\n"), file = out.file, sep = "\t", append = T)
    cat(c(Precision, "! 1/2/3=low/medium/high Precision for Full likelihood","\n"), file = out.file, sep = "\t", append = T)
    cat(c(row.names(vcf.input.df),"! Marker IDs (Loci)","\n"), file = out.file, sep = "\t", append = T)
    cat(c(MarkerTypes,"!  Marker Type: 0 co-dominant / 1 dominant, use '0@' or '1@' to repeat the option for all loci","\n"), file = out.file, sep = "\t", append = T)
    cat(c(AlleleDropRate,"! Allelic Drop-out rate, use 'x@' to repeat the option for all loci","\n"), file = out.file, sep = "\t", append = T)
    cat(c(OtherErrorRate,"! Other Error rate, use 'x@' to repeat the option for all loci","\n"), file = out.file, sep = "\t", append = T)

    write.table(x = offspring.gt, file = out.file, append = T, quote = F, sep = "\t", row.names = T, col.names = F)

    cat(c(ProbParentIncl, ProbParentIncl,"! Probability dad/mom is included","\n"), file = out.file, sep = "\t", append = T)
    cat(c(n.father, n.mother,"! Numbers of candidate males & females","\n","\n"), file = out.file, sep = "\t", append = T)

    write.table(x = father.gt, file = out.file, append = T, quote = F, sep = "\t", row.names = T, col.names = F)
    cat("\n", file = out.file, sep = "\t", append = T)

    write.table(x = mother.gt, file = out.file, append = T, quote = F, sep = "\t", row.names = T, col.names = F)
    cat("\n", file = out.file, sep = "\t", append = T)

    cat("
        0 0                                   ! Number of offspring with known father
        0 0                                   ! Number of offspring with known mother
        0                                    ! Number of known paternal sibships
        0                                    ! Number of known maternal sibships
        0                                    ! Number of offspring with known excluded fathers
        0                                    ! Number of offspring with known excluded mothers
        0                                    ! Number of offspring with known excluded paternal sibships
        0                                    ! Number of offspring with known excluded maternal sibships", file = out.file, sep = "\t", append = T)
}

