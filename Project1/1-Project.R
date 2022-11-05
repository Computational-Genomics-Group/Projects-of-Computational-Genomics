#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem
#TASK2: Calculates the minor allele frequency q for each SNP
#TASK3: Return the minor allele frequency of each SNP as a vector of numeric values with names corresponding 
#to the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the 
#input matrix
#Suggestion: 
#• It might be useful to use the function table() and to convert it in a data.frame(). 
#• Sometime you might have 0 subject with genotype aa… c
# SNPdata <- read.table("SNPdata.txt",header = TRUE, sep = "") * Just to visualize the data
# print(SNPdata[1,1])
qcalculation <- function(SNPdata){
  N <- ncol(SNPdata)
  res <- c(1:nrow(SNPdata))
  for(row in 1:nrow(SNPdata)) {
    AA <- 0
    Aa <- 0
    aa <- 0
    for(col in 1:ncol(SNPdata)) {
      ind_val <- SNPdata[row,col]
      if(ind_val == 0){AA = AA + 1}
      if(ind_val == 1){Aa = Aa + 1}
      if(ind_val == 2){aa = aa + 1}
    }
    
    q <- (aa*2 + Aa) / (2*N)
    res[row] <- q
  }
  
  return(res)
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
  N <- ncol(SNPdata)
  
 
  HWE_pvalue <- c(1:nrow(SNPdata))

  names(HWE_pvalue) <- row.names(SNPdata)
  
  for (row in 1:nrow(SNPdata)) {
    
    AA <- 0
    Aa <- 0
    aa <- 0
    
    for (col in 1:ncol(SNPdata)) {
      
      cell <- SNPdata[row,col]
      if(cell == 0){AA = AA + 1}
      if(cell == 1){Aa = Aa + 1}
      if(cell == 2){aa = aa + 1}
    }
    
    q = (aa*2+Aa)/(2*N)
    p = 1 - q
    Chi = ((AA-N*p^2)^2 )/(N*p^2) + ((Aa-2*N*p*q)^2 )/(2*N*p*q) + ((aa-N*q^2)^2 )/(N*q^2)
    
    pvalue <- pchisq(Chi, 1, ncp = 0, lower.tail = FALSE)
    
    HWE_pvalue[row] = pvalue
    
  }
  return(HWE_pvalue)
  
}

VARIANTanalysis <- function(filepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01){
  
}
