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

<<<<<<< Updated upstream
VARIANTanalysis <- function(filepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01){
=======
##################################
# FILTERING IDEA

# vector containing the height of 4 person
v <- c(150,160,170,180)
v
# name of the person
names(v) <- c("John","Sara","Paul","Carl")
v
# keep in the vector 'v' only person with height > 165
v = v [v>165]
v

# I will apply the same idea in the following function

################################
## Here we checked that by filtering the matrix the label in raw will remain correct
# toy example

SNP_sample = SNPdata[1:2,]
SNP_sample_new_1 = SNP_sample[c(TRUE, FALSE),]
##############################


VARIANTanalysis <- function(filepath=path, indCTRL=1201:2000, MAFth = 0.01, HWEalpha = 0.01){
>>>>>>> Stashed changes
  
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
  
<<<<<<< Updated upstream
  
  }
  }
VARIANTanalysis(path, 1201)


## ---------------------------------------------------------------------------------------------------------------------
## Test of the functions
setwd("C:/Users/Piermarco/Documents/GitHub/BigDataBuona/Projects-of-Computational-Genomics/Project1")
#setwd("C:/Users/Alessio/OneDrive - Università degli Studi di Padova/Bio_Eng/ComputationalGenomics/LAB/Groupwork1")
=======
  # Make the Contingency table and expected occurencies
  N = dim(SNPdata_patient_filt)[1]
  N = dim(SNPdata_ctrl_filt)[1]
  chromosome <- matrix(ncol = 8)
  for (row in 1:nrow(SNPdata)) {
    
    AA_patient <- 0
    Aa_patient <- 0
    aa_patient <- 0
    AA_ctrl <- 0
    Aa_ctrl <- 0
    aa_ctrl <- 0
      
      for (col in 1:ncol(SNPdata_ctrl_filt)) {
        
        cell <- SNPdata_ctrl_filt[row,col]
        if(cell == 0){AA_ctrl = AA_ctrl + 1}
        if(cell == 1){Aa_ctrl = Aa_ctrl + 1}
        if(cell == 2){aa_ctrl = aa_ctrl + 1}
      }
    
      for (col in 1:ncol(SNPdata_patient_filt)) {
      
      cell <- SNPdata_patient_filt[row,col]
      if(cell == 0){AA_patient = AA_patient + 1}
      if(cell == 1){Aa_patient = Aa_patient + 1}
      if(cell == 2){aa_patient = aa_patient + 1}
      }
    contingency_table <- matrix(c(AA_patient, Aa_patient, aa_patient, AA_ctrl, Aa_ctrl, aa_ctrl), nrow = 2, ncol = 3, byrow = TRUE)
    #colnames(contingency_table) <- c("AA", "Aa", "aa")
    #row.names(contingency_table) <- c("patients", "controls")
    N <- AA_patient + Aa_patient + aa_patient + AA_ctrl + Aa_ctrl + aa_ctrl
    
    expected_occurrences <- matrix(c(((AA_patient + AA_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                     ((Aa_patient + Aa_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                     ((aa_patient + aa_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                     ((AA_patient + AA_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl)),
                                     ((Aa_patient + Aa_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl)),
                                     ((aa_patient + aa_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl))),
                                     nrow = 2, ncol = 3, byrow = TRUE)
    
    Chi_obs = ((contingency_table[1,1] - expected_occurrences[1,1])^2/ expected_occurrences[1,1]) + 
      ((contingency_table[1,2] - expected_occurrences[1,2])^2/ expected_occurrences[1,2]) +
      ((contingency_table[1,3] - expected_occurrences[1,3])^2/ expected_occurrences[1,3]) +
      ((contingency_table[2,1] - expected_occurrences[2,1])^2/ expected_occurrences[2,1]) +
      ((contingency_table[2,2] - expected_occurrences[2,2])^2/ expected_occurrences[2,2]) +
      ((contingency_table[2,3] - expected_occurrences[2,3])^2/ expected_occurrences[2,3])
    
    pvalue_est <- pchisq(Chi_obs, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    qvalue_est <- qchisq(Chi_obs, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    
    chromosome <- rbind(c(contingency_table[2,1],contingency_table[2,2],contingency_table[2,3],
                      contingency_table[1,1],contingency_table[1,2],contingency_table[1,3],
                      pvalue_est, qvalue_est))
  }
  
 
contingency_table <- matrix(c(AA_patient, Aa_patient, aa_patient, AA_ctrl, Aa_ctrl, aa_ctrl), nrow = 2, ncol = 3, byrow = TRUE)


N <- AA_patient + Aa_patient + aa_patient + AA_ctrl + Aa_ctrl + aa_ctrl
expected_occurrences <- matrix(c(((AA_patient + AA_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                 ((Aa_patient + Aa_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                 ((aa_patient + aa_ctrl)/ N * (AA_patient + Aa_patient + aa_patient)),
                                 ((AA_patient + AA_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl)),
                                 ((Aa_patient + Aa_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl)),
                                 ((aa_patient + aa_ctrl)/ N * (AA_ctrl + Aa_ctrl + aa_ctrl))),
                                 nrow = 2, ncol = 3, byrow = TRUE)

Chi_obs = ((contingency_table[1,1] - expected_occurrences[1,1])^2/ expected_occurrences[1,1]) + 
        ((contingency_table[1,2] - expected_occurrences[1,2])^2/ expected_occurrences[1,2]) +
        ((contingency_table[1,3] - expected_occurrences[1,3])^2/ expected_occurrences[1,3]) +
        ((contingency_table[2,1] - expected_occurrences[2,1])^2/ expected_occurrences[2,1]) +
        ((contingency_table[2,2] - expected_occurrences[2,2])^2/ expected_occurrences[2,2]) +
        ((contingency_table[2,3] - expected_occurrences[2,3])^2/ expected_occurrences[2,3])

pvalue_est <- pchisq(Chi_obs, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

colnames(chrosome) <- c("AA_ctrl", "Aa_ctrl", "aa_ctrl", "AA_case", "Aa_case", "aa_case", "pval", "qval")

}
## ---------------------------------------------------------------------------------------------------------------------
## Test of the functions
setwd("C:/Users/Piermarco/Documents/GitHub/BigDataBuona/Projects-of-Computational-Genomics/Project1")
>>>>>>> Stashed changes
SNPdata <- read.table("SNPdata.txt",header = TRUE, sep = "")
path <- "C:/Users/Piermarco/Documents/GitHub/BigDataBuona/Projects-of-Computational-Genomics/Project1/SNPdata.txt"
q_vec <- qcalculation(SNPdata)
HWE_pvalue_vec <- HWEtest(SNPdata)

