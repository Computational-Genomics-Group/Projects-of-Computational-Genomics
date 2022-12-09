MvAplot <- function(exprData,
                    pdffilename) {
  samples.name = colnames(exprData)
  samples.number = dim(exprData)[2]
  
  pdf(pdffilename)
  par(mfrow = c(2, 2))
  
  for (i in 2:samples.number) {
    M = log2(exprData[, 1]) - log2(exprData[, i])
    A = 0.5 * (log2(exprData[, 1]) + log2(exprData[, i]))
    plot(
      A,
      M,
      xlab = "log-average | A",
      ylab = "log-ratio | M",
      main = sprintf("%s vs %s", samples.name[1], samples.name[i]),
      pch = 20,
      type = "p"
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
  # Scale data by their sequencing depth and multiply by 10^6 ##################
  samples.name = colnames(exprData)
  samples.number = dim(exprData)[2]
  samples.seqdepth = as.data.frame.vector(colSums(exprData, na.rm = TRUE))
  exprData.scaled = exprData / t(samples.seqdepth) * (10 ^ 6)
  
  # Calculate SizeFactor #######################################################
  genes.length = annot$Length
  samples.sizefactor = as.data.frame.vector(
    apply(exprData.scaled
          , MARGIN = 2
          , (function(x) return(sum(x * genes.length)))
    )
  )
  
  MvAplot(exprData.scaled, "output/mva.pdf")
  
  # Calculate SF ScaleFactor ###################################################
  # as average of M
  # * by trimming the most extreme values of A
  # * by taking trimmed means of M values
  SF.vector = as.data.frame.vector(rep(NA,samples.number-1), length.out=samples.number - 1,row.names = colnames(exprData)[2:samples.number], col.names= "SF") # init vector
  exprData.norm = exprData.scaled*0-1 # init exprData.norm 
  
  for (i in 2:samples.number) {
    if(i!=1){
      ## remember -> exprData.scaled = exprData / t(samples.depth) * (10 ^ 6)
      ## remember -> sample.depth = colSums(exprData)
      M = log2(exprData.scaled[, 1]) - log2(exprData.scaled[, i])
      A = 0.5 * (log2(exprData.scaled[, 1]) + log2(exprData.scaled[, i]))
      
      M.Atrimmed = M[(A > Atrim[1] & A < Atrim[2])] ## most extreme values of A
      SF = mean(M.Atrimmed,trim=Mtrim) ## trim extreme values of M 
      SF.vector[i - 1,] = SF
      
      rm(M.Atrimmed,M,A)
      
      ## trying to normalize data
      exprData.norm[i] = exprData.scaled[, i] + 2 ^ SF
    }else{
      exprData.norm[1]=exprData.scaled[1]
    }
  }
  
  MvAplot(exprData.norm, "output/mvanorm.pdf")
  
  
  ### part 2
  # scale the genes (in original scale) by their length and multiply by 10^3
  exprData.scaled = exprData / genes.length * (10 ^ 3)
  
  
  return(
   c(SF.vector, exprData.scaled)
  )
}


# TESTS

exprData = read.table("raw_trascr_count.txt")
annot = read.delim("raw_trascr_count_annot.txt",
                   header = TRUE ,
                   sep = "\t")
MvAplot(exprData, "output/MvAPlots.pdf")
output = TMMnorm(exprData, annot)
osf = output[1,]
odt = output[2,]


i=2
Atrim = c(0, 8)
Mtrim = 0.02 
