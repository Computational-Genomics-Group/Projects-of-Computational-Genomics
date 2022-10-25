#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem
#TASK2: Calculates the minor allele frequency q for each SNP
#TASK3: Return the minor allele frequency of each SNP as a vector of numeric values with names corresponding 
#to the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the 
#input matrix
#Suggestion: 
#• It might be useful to use the function table() and to convert it in a data.frame(). 
#• Sometime you might have 0 subject with genotype aa… c
SNPdata <- read.table("C:\\Users\\Piermarco\\Documents\\GitHub\\BigDataBuona\\Projects-of-Computational-Genomics\\Project1\\SNPdatasmall.txt",
                      header = TRUE, sep = "")
qcalculation <- function(SNPdata){
  for(col in 1:ncol(SNPdata)) {
    AA <- 0
    Aa <- 0
    aa <- 0
    for(row in 1:nrow(SNPdata)) {
      i <- SNPdata[row,col]
      if(i == 0){AA = AA + 1}
      if(i == 1){Aa = Aa + 1}
      if(i == 2){aa = aa + 1}
    }
    a_sum <- aa * 2 + Aa
    sum_all <- (aa + Aa + AA) * 2 
    q <- a_sum / sum_all
    print(row)
    print(q)
    row.names(SNPdata) <- row
    vec[row] <- c(q)
  }
 
  #both same result 
  "a_sum = aa + (Aa / 2)
  sum_all = aa + Aa + AA
  q = a_sum / sum_all"
}
#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem 
#TASK2: Compute a HWE test for each SNP given as input 
#• By calculating the chi^2 obs from the data
#• By computing the p value using the function pchisq (DO NOT use directly the chisq.test() function)
#TASK3: Return the HWE test p-values for each SNP as a vector of numeric values with names corresponding to 
#the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the input 
#matrix
#Suggestion: be careful when you use pchisq(). The probability it gives as output by default is P[X ≤ chi^2


HWEtest <- (SNPdata){
  
}

VARIANTanalysis <- function(filepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01){
  
}