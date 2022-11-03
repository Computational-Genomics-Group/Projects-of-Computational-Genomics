#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem
#TASK2: Calculates the minor allele frequency q for each SNP
#TASK3: Return the minor allele frequency of each SNP as a vector of numeric values with names corresponding 
#to the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the 
#input matrix
#Suggestion: 
#• It might be useful to use the function table() and to convert it in a data.frame(). 
#• Sometime you might have 0 subject with genotype aa… c

qcalculation <- function(SNPdata) {
  calcq <- function(d) {
    N = length(d)
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    q = (aa * 2 + Aa) / (2 * N)
    ##cat(AA,Aa,aa,N,q,"\n") #debug
    return(q)
  }
  
  a = apply(SNPdata, 1, calcq) # apply calcq to rows
  return(as.data.frame(a))
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
                                                                                                 

HWTest <- function(SNPdata) {
  calcchi <- function(d) {
    # 1. Calculate observed frequencies
    # 2. Calculate expected frequencies
    # 3. Calculate chi_value
    # 4. calculate p-value
    
    # d=t(SNPdata[1,])
    N = length(d)
    
    ## 1. Calculate observed values
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    O = c(AA, Aa, aa)
    q = (aa * 2 + Aa) / (2 * N)
    p = 1 - q
    cat(O,"\n")
    
    ## 2. Expected frequencies --------------------------------------------
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    #    m = rep(N, 3)
    #    E = c(AA, Aa, aa)
    #    E = round(E * m / N, digits = 0) # vector containing con#catenation of expected frequencies for control and patients
    #    cat(E,"\n")
    #    
    #    ## 3. Calculate chi-squared --------------------------------------------
    #    # method 1
    #    chi_squared = sum((O - E) ^ 2 / E)
    
    # method2
    prob_expec=c(p^2,2*p*q,q^2)
    chi_squared=sum((O-N*prob_expec)^2/(N*prob_expec))
    
    ## 4. Calculate pvalues -------------------------------------------------
    pvalue <- pchisq(chi_squared, 1, lower.tail = FALSE)
    
    
    return(c(chi_squared, pvalue))
  }
  
  ## Apply calcchi to rows of SNPdata
  # Transpose result in order to have chr as rows and [chisquare,pvalue] as columns
  a = as.data.frame(t(apply(SNPdata, 1, calcchi)))
  return(a)
}

VARIANTanalysis <-
  function(filepath,
           indCTRL,
           MAFth = 0.01,
           HWEalpha = 0.01) {
  ##### Input and parameters setup
  SNPdata <- read.table("SNPdata.txt", header = TRUE, sep = "\t")
  SNPdata = SNPdata[qcalculation(SNPdata) > MAFth |
                      HWTest(SNPdata)[2] > HWEalpha, ]
  `%notin%` <- Negate(`%in%`)
  N = dim(SNPdata)[2]
  indPATI = 1:N
  indPATI = indPATI[indPATI %notin% indCTRL]
  
  ##### calcchi and calcallele function definition
  t=HWTest(SNPdata[indCTRL])
  chisq = t[1]
  pvalues = t[2]
  maf = qcalculation(SNPdata)
  
  
  # "AA_ctrl","Aa_ctrl","aa_ctrl","AA_case","Aa_case","aa_case"
  calcallele <- function(d, indCTRL, indPATI) {
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    O = c(AA, Aa, aa)
    return(O)
  }
  
  O = as.data.frame(cbind(
    t(apply(SNPdata[, indCTRL], 1, calcallele)),
    t(apply(SNPdata[, indPATI], 1, calcallele))
  ))
  
  ##cal q values following Benjamini-Hochberg Procedure
  r=order(pvalues,decrease=FALSE)
  qvalues=pvalues*N/r
  
  ## Final Output
  out = as.data.frame(cbind(O, pvalues, qvalues))
  colnames(out) <-
    c("AA_ctrl",
      "Aa_ctrl",
      "aa_ctrl",
      "AA_case",
      "Aa_case",
      "aa_case",
      "pval",
      "qval")
  return(out)
}


SNPdata <- read.table("SNPdata.txt", header = TRUE, sep = "\t")
vec_q<-as.data.frame(qcalculation(SNPdata))
HWE_pvalue_vec <- HWEtest(SNPdata)
varanal=as.data.frame(VARIANTanalysis("SNPdata.txt",1201:2000,MAFth = 0.05))
