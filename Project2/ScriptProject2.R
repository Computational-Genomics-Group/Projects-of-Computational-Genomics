MvAplot <- function(exprData, 
                    pdffilename){
  namesamples=colnames(exprData)
  nsamples=dim(exprData)[2]
  
  pdf(pdffilename)
  par(mfrow = c(2,2))
  
  for(i in 2:nsamples){
    M=log2(exprData[,1]/exprData[,i])
    A=0.5*log2(exprData[,1]*exprData[,i])
    plot(A, M,xlab = "A",ylab = "M",
         main = sprintf("%s vs %s",namesamples[1],namesamples[i]))
    abline(h=0)
  }
  
  dev.off()
}









# TESTS

exprData=read.table("raw_trascr_count.txt")
MvAplot(exprData, "MvAPlots.pdf")

