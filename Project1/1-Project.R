# TASK1 -------------------------------------------------------------------
# Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
# provided in stem 


qcalculation <- function(SNPdata) {
  calcq <- function(d) {
    N = length(d)
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    q = min( ((aa * 2 + Aa) / (2 * N)),
             ((AA * 2 + Aa) / (2 * N)),
             na.rm = TRUE
    )
    return(q)
  }
  
  out = as.data.frame(apply(SNPdata, 1, calcq)) 
  colnames(out)="MAF"
  return(out)
}


# TASK2 -------------------------------------------------------------------
# Compute a HWE test for each SNP given as input 
# • By calculating the chi^2 obs from the data
# • By computing the p value using the function pchisq (DO NOT use directly the chisq.test() function)
                                                                                                 

HWEtest <- function(SNPdata) {
  calcp<- function(d) {
    N = length(d)
    
    ## Calculate observed values -------------------------------------------
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    O = c(AA, Aa, aa)
    
    ## Expected frequencies ------------------------------------------------
    q = (aa * 2 + Aa) / (2 * N)
    p = 1 - q
    prob_expec=c(p^2,2*p*q,q^2)
    
    ## ChiSquared ----------------------------------------------------------
    chi_squared=sum((O-N*prob_expec)^2/(N*prob_expec))
    
    ## pvalue --------------------------------------------------------------
    pvalue <- pchisq(chi_squared, 1, lower.tail = FALSE)
    
    return(pvalue)
  }
  
  out = as.data.frame(apply(SNPdata, 1, calcp))
  colnames(out)=c("pvalue")
  return(out)
}

# TASK3 -------------------------------------------------------------------
# Return the HWE test p-values for each SNP as a vector of numeric values with names corresponding to 
# the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the input 
# matrix
# Suggestion: be careful when you use pchisq(). The probability it gives as output by default is P[X ≤ chi^2

VARIANTanalysis <-
  function(filepath,
           indCTRL,
           MAFth = 0.01,
           HWEalpha = 0.01) {
  ## Functions definitions 
  calcallele <- function(d, indCTRL, indPATI) {
    AA = length(d[d == 0])
    Aa = length(d[d == 1])
    aa = length(d[d == 2])
    O = c(AA, Aa, aa)
    return(O)
  }
  
  calcchip<- function(d, indCTRL, indPATI) {
    # return chisquared and pvalue 
    N = length(d)
    O <- c(calcallele(d[indCTRL]),    
           calcallele(d[indPATI]))   
    rows=c(rep(length(indCTRL),3),
           rep(length(indPATI),3)) 
    columns=rep(c(O[1]+O[4], #AA
                  O[2]+O[5], #Aa
                  O[3]+O[6]) #aa
                ,2)
    E= round(rows*columns/N,digits = 0)           
    chi_squared=sum((O-E)^2/E,na.rm = TRUE)
    if (O[3]!=0 | O[6]!=0)
      df=2
    else
      df=1
    pvalue <- pchisq(chi_squared, df, lower.tail = FALSE)
    
    return(c(chi_squared, pvalue))
  }
  
  `%notin%` <- Negate(`%in%`)
  
  ## Input and parameters setup ----------------------------------------------
  SNPdata <- read.table("SNPdata.txt", header = TRUE, sep = "\t")
  SNPdata = SNPdata[qcalculation(SNPdata) > MAFth & 
                    HWEtest(SNPdata[,indCTRL]) > HWEalpha , ]
  N = dim(SNPdata)[2]
  indPATI = (1:N)[(1:N) %notin% indCTRL]
  
  ## Calculations ------------------------------------------------------------
  O = as.data.frame(cbind(                            #observed data
        t(apply(SNPdata[, indCTRL], 1, calcallele)),
        t(apply(SNPdata[, indPATI], 1, calcallele))
      ))
  
  t=as.data.frame(t((apply(SNPdata, 1, calcchip, indCTRL=indCTRL, indPAT=indPATI))))
  chisquared=t[,1]
  pvalues = t[,2]
  
  r = order(unlist(pvalues), decreasing = FALSE)      
  qvalues = pvalues * N / r
  
  ## Final Output
  out = as.data.frame(cbind(O, pvalues, qvalues))
  colnames(out) <-c("AA_ctrl","Aa_ctrl","aa_ctrl","AA_case","Aa_case","aa_case","pval","qval")
  return(as.data.frame(out))
}


# tests -------------------------------------------------------------------
SNPdata <- read.table("SNPdata.txt", header = TRUE, sep = "\t")
vec_maf <- qcalculation(SNPdata = SNPdata)
vec_HWE_chi_pvalue <- HWEtest(SNPdata)
vec_VarAnal = VARIANTanalysis("SNPdata.txt", 1201:2000, MAFth = 0.05)

