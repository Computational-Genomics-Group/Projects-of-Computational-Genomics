MvAplot <- function(exprData, 
                    pdffilename){
  namesamples=colnames(exprData)
  nsamples=dim(exprData)[2]
  
  pdf(pdffilename)
  par(mfrow = c(2,2))
  
  for(i in 2:nsamples){
    M=log2(exprData[,1]/exprData[,i])
    A=0.5*log2(exprData[,1]*exprData[,i])
    plot(A, M,xlab = "A",ylab = "M"
         ,main = sprintf("%s vs %s",namesamples[1],namesamples[i])
         ,pch=20
         )
    abline(h=0)
  }
  
  dev.off()
}

TMMnorm <- function(exprData = RawTranscript,
                    annot = RawTranscript_annot,
                    Atrim = c(0,8),
                    Mtrim = 0.02){
  samples.names=as.data.frame.vector(colnames(exprData))
  samples.number=dim(exprData)[2]
  genes.names=as.data.frame.vector(row.names(exprData))
  genes.number=dim(exprData)
  genes.level=as.data.frame.vector(rowSums(exprData))
  genes.length=annot$Length
  samples.depth=as.data.frame.vector(colSums(exprData))
  
  # Scale data by their sequencing depth and multiply by 10^6
  exprData.scaled=exprData/t(samples.depth)*(10^6)
  
  # Calculate SF with respect to sample1
  # by trimming the most extreme values of A
  # by taking trimmed means of M values
  samples.sizefactor=as.data.frame.vector(
    apply(
      exprData.scaled,
      MARGIN = 2,
      (function(x) return(sum(x*genes.length)))
    )
  )
  
  MvAplot(exprData.scaled,"mva.pdf")
  
  # Calculate M and A (using scale factor in order to make possible comparisons between different samples)
  # use instead of n1 -> n1/S1
  exprData.norm=exprData.scaled
  SF.vector=rep(NaN,samples.number-1)
  for(i in 2:samples.number){
   # M=log2(exprData.norm[,1]/samples.sizefactor[1,]*(samples.sizefactor[i,]/exprData.norm[,i]))
   # A=0.5*(log2(exprData.norm[,1]/samples.sizefactor[1,])+log2(exprData.norm[,i]/samples.sizefactor[i,]))
    M=log2(exprData.scaled[,1]/exprData.scaled[,i])
    A=0.5*(log2(exprData.scaled[,1]*exprData.scaled[,i]))
    temp=sort(M)
    M.boundaries=c(temp[ceiling(1+Mtrim*length(M))],
                   temp[floor(length(M)*(1-Mtrim))])
    M.trimmed=M[(A>Atrim[1] & A<Atrim[2]) & (M>M.boundaries[1] & M<M.boundaries[2])]
    SF.vector[i-1]=mean(M.trimmed)
    exprData.norm[,i]=exprData.scaled[,i]+2^SF
    
  }
  
  MvAplot(exprData.norm,"mvanorm.pdf")
  
  
  ### part 2
  # scale the genes (in original scale) by their lenght and multiply by 10^3
  exprData.scaled.lenghtgene=exprData/genes.length*(10^3)
  
  
  return(matrix(c(SF.vector, exprData.scaled.lenghtgene))
}








# TESTS

exprData=read.table("raw_trascr_count.txt")
annot=read.delim("raw_trascr_count_annot.txt", header = TRUE ,sep = "\t")
MvAplot(exprData, "MvAPlots.pdf")
output=TMMnorm(exprData,annot)
osf=output[1,]
odt=output[2,]

