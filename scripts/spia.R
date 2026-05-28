############################################################################
#
# Compute SPIA test
#  Input: 1. x: a matrix with a column for each cell line and a row for each SNP
#            in the SPIA format (use toSPIAData before)
#         2. row.names: specify if the fisrt column contains the name of the SNPs
#         3. test.prob: specify if probabilistic test has to be performed
#         4. test.param: specify the parameters of the test
#            - test.param$Pmm: SNP probability of mismatch in a matching population (dafault 0.1)
#            - test.param$nsigma: area limit for Pmm
#            - test.param$Pmm_nonM: SNP probability of mismach in a non matching population (dafault 0.6)
#            - test.param$nsigma_nonM: area limit for Pmm_nonM
#            - test.param$PercValidCall: percentage of valid call to consider the test valid
#         5. verbose: print verbose information
#  Output: a matric with a line for each cell line and with columns with the
#          informationss about distances 
#
#############################################################################
# rm(list=ls());verbose<-TRUE; options(max.print=1000);getOption("max.print")


install.packages("https://cran.r-project.org/src/contrib/Archive/SPIAssay/SPIAssay_1.1.0.tar.gz")

## set your working dir
# setwd("HandsOn_SPIA(1)"); getwd();

library(SPIAssay)
library(dplyr)

SPIAtable_from_vcf <- function(vcf_file) {
    #
    # TEST 1
    # 
    # read two times the vcf file, first for the columns names, second for the data
    tmp_vcf<-readLines(vcf_file)
    tmp_vcf_data<-read.table(vcf_file, stringsAsFactors = FALSE)

    # filter for the columns names
    tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
    vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
    names(tmp_vcf_data)<-vcf_names
    #
    # TEST 2
    #
    # Read the VCF file
    vcf_data <- read.table(vcf_file, header = TRUE, sep = "\t", as.is = TRUE)
    
    # Extract relevant columns (assuming the first column is SNP ID and the rest are genotypes)
    snp_ids <- vcf_data[, 1]
    genotypes <- vcf_data[, -1]
    
    # Create a matrix for SPIA input
    spia_matrix <- as.matrix(genotypes)
    rownames(spia_matrix) <- snp_ids
    
    return(spia_matrix)
}

args <- commandArgs(trailingOnly = TRUE)
input_table <- args[1]

## set param for high MAF SNPs
spia_parameters <- list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)

## load genotype data
genotype_data <- read.table(input_table, header = TRUE, sep = "\t", as.is = TRUE ); 
dim(genotype_data); head(genotype_data);

## load high MAF SNPs unique IDs
selected_snps <- read.delim("SPIA_selected_SNPs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
length(selected_snps); head(selected_snps);

##TASK 1
## select genotype data based on SNP IDs
genotype_data_selected <- genotype_data[selected_snps,]
dim(genotype_data_selected);

## run test
SPIA_selected <-  SPIATest(x = genotype_data_selected, row.names = FALSE, test.param = spia_parameters) 

test_result <- SPIA_selected$SPIAresult
return test_result