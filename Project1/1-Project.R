#TASK1: Take as input a numeric data matrix that is supposed to have the same format of the genetic data 
#provided in stem
#TASK2: Calculates the minor allele frequency q for each SNP
#TASK3: Return the minor allele frequency of each SNP as a vector of numeric values with names corresponding 
#to the SNP IDs (chromosome_position, e.g. “chr1_11169676”) with the same order they had in the 
#input matrix
#Suggestion: 
#• It might be useful to use the function table() and to convert it in a data.frame(). 
#• Sometime you might have 0 subject with genotype aa… c

qcalculation <- function(SNPdata){
  calcq<- function(d){
    t=table(d)  # if #1,#2,#3>0 then return integer[3], if #3=0-> (2<->aa does not exist) it return  returns a integer[2]
    AA=t[1]
    Aa=t[2]
    aa=if(is.na(t[3])) 0 else t[3] # if integer[2] apply correction adding 0
    N=length(d)
    q=(aa*2+Aa)/(2*N)
    #cat(AA,Aa,aa,N,q,"\n") #debug
    return(q)
  }
  
  out=apply(SNPdata,2, calcq) # apply calcq to rows
  return(out)
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
  calcchi<- function(d){
    # 1. Calculate observed frequencies
    # 2. Calculate expected frequencies
    # 3. Calculate chi_value
    # 4. calculate p-value
    
    N=length(d)
    
    ## 1. Calculate observed values for control patients -------------------
    # control subjects 
    t=table(d[1201:2000])
    AA=t[1]
    Aa=t[2]
    aa=if(is.na(t[3])) 0 else t[3] 
    O=c(AA,Aa,aa)
    q_obs_ctrl=(aa*2+Aa)/(2*N)
    p_obs_ctrl=1-q_obs_ctrl
    cat("O1",O,"\n")
    
    # patients
    t=table(d[1:1200])
    AA=t[1]
    Aa=t[2]
    aa=if(is.na(t[3])) 0 else t[3] 
    O=c(O,c(AA,Aa,aa))  # vector containing concatenation of observed frequencies for control and patients
    q_obs_pat=(aa*2+Aa)/(2*N)
    p_obs_pat=1-q_obs_pat
    cat("O2",O,"\n")
    
    
    ## 2. Expected frequencies --------------------------------------------
    t=table(d)
    AA=t[1]
    Aa=t[2]
    aa=if(is.na(t[3])) 0 else t[3] 
    t=c(AA,Aa,aa)
    cat("t",t,"\n")
    E=c(t,t)
    m=c(rep(800,3),rep(1200,3))
    E=round(E*c(rep(2000-1200,3),rep(1200,3))/N,digits = 0) # vector containing concatenation of expected frequencies for control and patients
    E=round(E*m/N,digits = 0) # vector containing concatenation of expected frequencies for control and patients
    cat("E",E,"\n")
    
    ## 3. Calculate chi-squared --------------------------------------------
    # method 1
    chi_squared=sum((O-E)^2/E)
    cat("chi",chi_squared,"\n")
    
    ## method2
    #prob_expec=c(p_obs_ctrl^2,2*p_obs_ctrl*q_obs_ctrl,q_obs_ctrl^2,
    #             p_obs_pat^2,2*p_obs_pat*q_obs_pat,q_obs_pat^2)
    #chi_squared=sum((O-N*prob_expec)^2/(N*prob_expec))
    
    ## 4. Calculate pvalues -------------------------------------------------
    pvalue <- pchisq(chi_squared,1,lower.tail = FALSE)
    cat("p",pvalue,"\n")
    
    
    return(c(chi_squared, pvalue))
  }
  
  ## Apply calcchi to rows of SNPdata
  # Transpose result in order to have chr as rows and [chisquare,pvalue] as columns
  out=as.data.frame(t(apply(SNPdata,1, calcchi))) 
  return(out)
}

VARIANTanalysis <- function(filepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01){
  SNPdata <- read.table("SNPdata.txt", header = TRUE, sep = "\t")
  SNPdata=SNPdata[qcalculation(SNPdata)>MAFth | HWTest(SNPdata)[2]>HWEalpha,]
  nam<-c("AA_ctrl","Aa_ctrl","aa_ctrl","AA_case","Aa_case","aa_case","pval","qval")
  
  res=as.data.frame(HWTest(SNPdata))
  return(res)
}