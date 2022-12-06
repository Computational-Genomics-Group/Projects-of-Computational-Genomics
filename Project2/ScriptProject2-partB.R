library(edgeR)

DATA <-
  read.table(
    "raw_trascr_count.txt",
    sep = "\t",
    row.names = 1,
    header = TRUE
  )
indGroup <-
  c(grep("Group1", colnames(DATA)), grep("Group2", colnames(DATA)))
indCTRL <- grep("CTRL", colnames(DATA))
lg <- length(indGroup)
lc <- length(indCTRL)

target = as.data.frame(rep("Patients", dim(DATA)[2]), row.names = colnames(DATA))
colnames(target) = c("Groups")
target[indCTRL, ] = "CTRL"
group = factor(target$Groups)
design = model.matrix( ~ 0 + group)


y <- DGEList(counts = DATA)    # y is an object of type DGE
y <-
  calcNormFactors(y)   # This calculates the SF using the TMM normalization !!!
SF <- y$samples

y <- 
  estimateGLMCommonDisp(y, design, verbose = TRUE) #phi common to the entire dataset
y <- estimateGLMTrendedDisp(y, design) #phi depends on mu
y <- estimateGLMTagwiseDisp(y, design) #phi is gene specific
fit <-
  glmFit(y, design) #finally the model fit (that accounts for raw NB data and scaling factors and seq. depth)
summary(fit)

Confronti <-
  makeContrasts(Treatment = "groupCTRL-groupPatients", levels = design)
RES <- glmLRT(fit, contrast = Confronti[, "Treatment"])
# The first coloumn of RES reports the log_Fold_Change, i.e.:
# log2(Normalized_data_average_GroupSARSCoV2 / Normalized_data_average_GroupMock)

alpha = 0.05

selected = sum(RES$table$PValue < alpha)
G=dim(DATA)[1]
G0=0.8*G

EFP=min(selected, G0*alpha)
ETP=max(0,selected-EFP)
ETN=min(G-selected, G0-EFP)
EFN=G-EFP-ETP-ETN
FDR=EFP/selected

tablehptest=as.data.frame(matrix(c(ETN,EFP,G0,EFN,ETP,G-G0,G-selected,selected,G),nrow=3,byrow = TRUE),row.names = c("H0","H1"," "))
colnames(tablehptest)=c("H0","H1"," ")


salpha = selected.Group + selected.CTRL
EFP = min(sum(DATA) * 0.8 * alpha, salpha)
EFN = sum(DATA) - EFP - max(0, salpha - EFP) - min(sum(DATA) - salpha, sum(DATA) *
                                                     0.8 - EFP)
FDR = EFP / salpha


pvalues.vector=RES$table$PValue
pvalues.rank=order(pvalues.vector)
qvalues.benjhock=pvalues.vector*G/pvalues.rank

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
##################################################################################

DEbyEdgeR <- function(rawdat, groups, alpha = 0.05) {
  return(c(selected.genes.alpha, EFP))
}





alpha = 0.05
rawdat = read.table("raw_trascr_count.txt")
groups = colnames(rawdat)

DEbyEdgeR(rawdat, groups)