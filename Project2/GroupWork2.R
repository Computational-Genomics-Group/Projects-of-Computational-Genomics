# MvAPlot ######################################################################
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
  samples.seqdepth = as.vector(colSums(exprData, na.rm = TRUE))
  exprData.scaled = exprData / t(samples.seqdepth) * (10 ^ 6)
  
  # Calculate SizeFactor #######################################################
  genes.length = annot$Length
  samples.sizefactor = apply(exprData.scaled
          , MARGIN = 2
          , (function(x) return(sum(x * genes.length)))
    )
  
  MvAplot(exprData.scaled, "output/mvabeforenorm.pdf")
  
  # Calculate SF ScaleFactor ###################################################
  # as average of M
  # * by trimming the most extreme values of A
  # * by taking trimmed means of M values
  SF.vector = as.vector(rep(NA,samples.number-1))
  exprData.norm = exprData.scaled*0-1 # init exprData.norm 
  
  for (i in 1:samples.number) {
    if(i!=1){
      M = log2(exprData.scaled[, 1]) - log2(exprData.scaled[, i])
      A = 0.5 * (log2(exprData.scaled[, 1]) + log2(exprData.scaled[, i]))
      
      M.Atrimmed = M[(A > Atrim[1] & A < Atrim[2])] ## most extreme values of A
      SF = mean(M.Atrimmed,trim=Mtrim) ## trim extreme values of M 
      SF.vector[i - 1] = SF
      
      rm(M.Atrimmed,M,A)
      
      ## trying to normalize data
      exprData.norm[i] = exprData.scaled[, i] + 2 ^ SF
    }else{
      exprData.norm[1]=exprData.scaled[1]
    }
  }
  
  MvAplot(exprData.norm, "output/mvaafternorm.pdf")
  
  
  # scale the genes (in original scale) by their length and multiply by 10^3
  exprData.scaled = exprData / genes.length * (10 ^ 3)
  
  return(list(
    SF.vector
    , exprData.scaled
  ))
}



#################################################################################
DEbyEdgeR <- function(rawdat, groups, alpha = 0.05) {
  indGroup <- c(
    grep("Group1", colnames(DATA)),
    grep("Group2", colnames(DATA))
  )
  lg <- length(indGroup)
  
  indCTRL <- grep("CTRL", colnames(DATA))
  lc <- length(indCTRL)
  
  target = matrix( rep("Patients", dim(DATA)[2]) )
  target[indCTRL, ] = "CTRL"
  rownames(target)=colnames(DATA)
  colnames(target) = c("Groups")
  
  Group = factor(target)
  design = model.matrix( ~ 0 + Group)
  
  y <- DGEList(counts = DATA)
  y <- calcNormFactors(y)
  SF <- y$samples
  
  y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design) 
  #summary(fit)
  
  
  Confronti <- makeContrasts(Treatment = "GroupCTRL-GroupPatients", levels = design)
  RES <- glmLRT(fit, contrast = Confronti[, "Treatment"])
  # The first coloumn of RES reports the log_Fold_Change, i.e.:
  # log2(Normalized_data_average_GroupSARSCoV2 / Normalized_data_average_GroupMock)
  
  G=dim(DATA)[1]
  G0=0.8*G
  
  alpha = 0.05
  selected = sum(RES$table$PValue < alpha)
  
  EFP=min(selected, G0*alpha)
  ETP=max(0,selected-EFP)
  ETN=min(G-selected, G0-EFP)
  EFN=G-EFP-ETP-ETN
  FDR=EFP/selected
  
  pvalues.vector=RES$table$PValue
  pvalues.rank=order(pvalues.vector)
  qvalues.benjhock=pvalues.vector*G/pvalues.rank
  
  # Output
  out=list(
    as.vector(
      c(selected,
        EFP,
        EFN,
        FDR),
      )
    , matrix(c(pvalues.vector
               ,qvalues.benjhock),
             ncol = 2,
             byrow=FALSE,  
             dimnames=list(row.names(RES$table),
                           c("pvalues","qvalues.BenjHock"))
             )
    )
  return(out)
}