#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem
#TASK2: Calculates the minor allele frequency q for each SNP
#TASK3: Return the minor allele frequency of each SNP as a vector of numeric values with names corresponding 
#to the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the 
#input matrix
#Suggestion: 
#• It might be useful to use the function table() and to convert it in a data.frame(). 
#• Sometime you might have 0 subject with genotype aa… consider this possible

qcalculation <- function(SNPdata) {
  
  # num of chromo
  N = dim(SNPdata)[2]
  
  # initialization of output vector
  q_vec = c(1:nrow(SNPdata))*0
  # labeling of the vector
  names(q_vec)<-row.names(SNPdata)
  
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
    
    q_est = (aa*2+Aa)/(2*N)
    q_vec[row] = q_est
    
  }
  return(q_vec)
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

HWEtest <- function(SNPdata){
  
  N = dim(SNPdata)[1]
 
  # initialization of output vector
  HWE_pvalue_vec = c(1:nrow(SNPdata))*0
  
  # labeling of the vector
  names(HWE_pvalue_vec)<-row.names(SNPdata)
  
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
    
    q_est = (aa*2+Aa)/(2*N)
    p_est = 1 - q_est
    Chi_est = ((AA-N*p_est^2)^2 )/(N*p_est^2) + ((Aa-2*N*p_est*q_est)^2 )/(2*N*p_est*q_est) + ((aa-N*q_est^2)^2 )/(N*q_est^2)
    
    pvalue_est <- pchisq(Chi_est, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    
    HWE_pvalue_vec[row] = pvalue_est
    
  }
  return(HWE_pvalue_vec)
  
}

VARIANTanalysis <- function(filepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01){
  
  # read data from file path
  SNPdata <- read.table(filepath, header = TRUE, sep = "")
  ## FILTERING OF THE RAW DATA
  
  # Selection of control group
  SNPdata_ctrl = SNPdata[ , c(indCTRL:2000)]
  
  # Check HWE assumption of th control group
  HWE_pvalue_ctrl <- HWEtest(SNPdata_ctrl)
  
  # Filter data with  HWE > HWEalpha 
  # !!!WARNING LOGIC CHECK I HAVE CHANGED FROM --->SNPdata_filt = SNPdata[HWE_pvalue_ctrl > HWEalpha, ] TO -->!!!
  #with 0 condition take 0 element
  vec <- which(HWE_pvalue_ctrl > HWEalpha)
  SNPdata_filt = SNPdata_ctrl[which(HWE_pvalue_ctrl > HWEalpha), ]
  
  # Chech MAF assumption
  q_vec <- qcalculation(SNPdata_filt)
  
  # Filter data  with MAF > MAFth 
  SNPdata_filt = SNPdata_filt[q_vec > MAFth, ]
  
  # calculation of the Chi_2 of the association test from the data
  
  N = dim(SNPdata_filt)[2]
  
  # initialization of output vector
  HWE_pvalue_vec = c(1:nrow(SNPdata_filt))*0
  
  # labeling of the vector
  names(HWE_pvalue_vec)<-row.names(SNPdata_filt)
  
  for (row in 1:nrow(SNPdata_filt)) {
    
    AA <- 0
    Aa <- 0
    aa <- 0
    
    for (col in 1:ncol(SNPdata_filt)) {
      
      cell <- SNPdata_filt[row,col]
      if(cell == 0){AA = AA + 1}
      if(cell == 1){Aa = Aa + 1}
      if(cell == 2){aa = aa + 1}
    }
    
    q_est_filt = (aa*2+Aa)/(2*N)
    p_est_filt = 1 - q_est
    Chi_est_filt = ((AA-N*p_est^2)^2 )/(N*p_est^2) + ((Aa-2*N*p_est*q_est)^2 )/(2*N*p_est*q_est) + ((aa-N*q_est^2)^2 )/(N*q_est^2)
    
  # calutation of the p_value
    
    P_value_filt <- pchisq(Chi_est_filt, 1, lower.tail = FALSE)
  
  
  }
  }
VARIANTanalysis(path, 1201)


## ---------------------------------------------------------------------------------------------------------------------
## Test of the functions
setwd("C:/Users/Piermarco/Documents/GitHub/BigDataBuona/Projects-of-Computational-Genomics/Project1")
#setwd("C:/Users/Alessio/OneDrive - Università degli Studi di Padova/Bio_Eng/ComputationalGenomics/LAB/Groupwork1")
SNPdata <- read.table("SNPdata.txt",header = TRUE, sep = "")
path <- "C:/Users/Piermarco/Documents/GitHub/BigDataBuona/Projects-of-Computational-Genomics/Project1/SNPdata.txt"
q_vec <- qcalculation(SNPdata)
HWE_pvalue_vec <- HWEtest(SNPdata)

