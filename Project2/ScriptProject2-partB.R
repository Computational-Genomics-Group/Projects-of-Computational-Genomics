DEbyEdgeR <- function(rawdat, groups, alpha = 0.05) {
  indGroup <- c(
    grep("Group1", colnames(DATA)),
    grep("Group2", colnames(DATA))
  )
  lg <- length(indGroup)

  indCTRL <- grep("CTRL", colnames(DATA))
  lc <- length(indCTRL)

  target = as.data.frame(
    rep("Patients", dim(DATA)[2]), row.names = colnames(DATA)
  )
  target[indCTRL, ] = "CTRL"
  colnames(target) = c("Groups")

  Group = factor(target$Groups)
  design = model.matrix( ~ 0 + Group)

  y <- DGEList(counts = DATA)    # y is an object of type DGE
  y <- calcNormFactors(y)   # This calculates the SF using the TMM normalization !!!
  SF <- y$samples

  y <- estimateGLMCommonDisp(y, design, verbose = TRUE) #phi common to the entire dataset
  y <- estimateGLMTrendedDisp(y, design) #phi depends on mu
  y <- estimateGLMTagwiseDisp(y, design) #phi is gene specific
  fit <- glmFit(y, design) #finally the model fit (that accounts for raw NB data and scaling factors and seq. depth)
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

  #tablehptest=as.data.frame(matrix(c(ETN,EFP,G0,EFN,ETP,G-G0,G-selected,selected,G),nrow=3,byrow = TRUE),row.names = c("H0","H1"," "))
  #colnames(tablehptest)=c("H0","H1"," ")
  #View(tablehptest)

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
      mode = "any"
    ),
    matrix(c(pvalues.vector,qvalues.benjhock),
          ncol = 2,
          byrow=FALSE,
          dimnames=list(row.names(RES$table),
                        c("pvalues","qvalues.BenjHock"))
    )
  )
  return(out)
}

# Tests
library(edgeR)
DATA <-
  read.table(
    "raw_trascr_count.txt",
    sep = "\t",
    row.names = 1,
    header = TRUE
  )

rawdat = read.table("raw_trascr_count.txt")
out=DEbyEdgeR(rawdat, groups)
out
