MvAplot <- function(exprData,
                    pdffilename) {
  namesamples = colnames(exprData)
  nsamples = dim(exprData)[2]
  
  pdf(pdffilename)
  par(mfrow = c(2, 2))
  
  for (i in 2:nsamples) {
    M = log2(exprData[, 1] / exprData[, i])
    A = 0.5 * log2(exprData[, 1] * exprData[, i])
    plot(
      A,
      M,
      xlab = "A",
      ylab = "M"
      ,
      main = sprintf("%s vs %s", namesamples[1], namesamples[i])
      ,
      pch = 20
    )
    abline(h = 0)
  }
  
  dev.off()
}

# TMMnorm ######################################################################
TMMnorm <- function(exprData = RawTranscript,
                    annot = RawTranscript_annot,
                    Atrim = c(0, 8),
                    Mtrim = 0.02) {
  #samples.names=as.data.frame.vector(colnames(exprData))
  #genes.number=dim(exprData)
  #genes.names=as.data.frame.vector(row.names(exprData))
  #genes.level=as.data.frame.vector(rowSums(exprData))
  
  # Scale data by their sequencing depth and multiply by 10^6 ##################
  samples.number = dim(exprData)[2]
  samples.depth = as.data.frame.vector(
    colSums(exprData)
    )
  exprData.scaled = exprData / t(samples.depth) * (10 ^ 6)
  genes.length = annot$Length
  
  # Calculate SizeFactor #######################################################
  samples.sizefactor = as.data.frame.vector(
    apply(exprData.scaled
          ,MARGIN = 2
          ,(function(x) return(sum(x * genes.length)))
          )
    )
  
  MvAplot(exprData.scaled, "mva.pdf")
  
  # Calculate SF ScaleFactor ###################################################
  # as average of M
  # * by trimming the most extreme values of A
  # * by taking trimmed means of M values
  SF.vector = rep(NaN, samples.number - 1) 
  
  for (i in 2:samples.number) {
    # M=log2(exprData.norm[,1]/samples.sizefactor[1,]*(samples.sizefactor[i,]/exprData.norm[,i]))
    # A=0.5*(log2(exprData.norm[,1]/samples.sizefactor[1,])+log2(exprData.norm[,i]/samples.sizefactor[i,]))
    ## remember -> exprData.scaled = exprData / t(samples.depth) * (10 ^ 6)
    ## remember -> sample.depth = colSums(exprData)
    M = log2(exprData.scaled[, 1] / exprData.scaled[, i])
    A = 0.5 * (log2(exprData.scaled[, 1] * exprData.scaled[, i]))
    temp = sort(M)
    M.boundaries = c(temp[ceiling(1 + Mtrim * length(M))],
                     temp[floor(length(M) * (1 - Mtrim))])
    M.trimmed = M[(A > Atrim[1] & A < Atrim[2])
                  & 
                  (M > M.boundaries[1] & M < M.boundaries[2])
                ]
    SF.vector[i - 1] = mean(M.trimmed)
    
    ## trying to normalize data
    exprData.norm[, i] = exprData.scaled[, i] + 2 ^ SF
    
  }
  
  MvAplot(exprData.norm, "mvanorm.pdf")
  
  
  ### part 2
  # scale the genes (in original scale) by their length and multiply by 10^3
  exprData.scaled.lenghtgene = exprData / genes.length * (10 ^ 3)
  
  
  # return(
  #  matrix(c(SF.vector, exprData.scaled.lenghtgene))
  # )
}








# TESTS

exprData = read.table("raw_trascr_count.txt")
annot = read.delim("raw_trascr_count_annot.txt",
                   header = TRUE ,
                   sep = "\t")
MvAplot(exprData, "MvAPlots.pdf")
output = TMMnorm(exprData, annot)
osf = output[1, ]
odt = output[2, ]
