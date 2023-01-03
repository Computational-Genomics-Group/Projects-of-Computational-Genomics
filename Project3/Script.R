
# PACKAGES ----------------------------------------------------------------

if(!require(compositions)) 
  install.packages(compositions)
library(compositions)

if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
library(metagenomeSeq)


if (!requireNamespace("GUniFrac", quietly = TRUE))  install.packages("GUniFrac")
library("GUniFrac")

if (!requireNamespace("glmnet", quietly = TRUE))  install.packages("glmnet")
if (!requireNamespace("devtools", quietly = TRUE))  install.packages("devtools")
if (!requireNamespace("Matrix", quietly = TRUE))  install.packages("Matrix")
library(devtools)
if (!requireNamespace("mbImpute", quietly = TRUE))  install_github("ruochenj/mbImpute/mbImpute R package")
library(mbImpute)
library(glmnet)
library(Matrix)

# ReadData ----------------------------------------------------------------
filename="Genus_otu_table.txt"
DATA=read.table(filename,header = TRUE, sep = "\t")
#DATA = DATA[1:4, 1:3]

# DATA TRANSFORMATIONS ----------------------------------------------------

data.alr=as.matrix(alr(DATA));
data.clr=as.matrix(clr(DATA));


# NORMALIZATION -----------------------------------------------------------
## CSS
metaSeqObject = newMRexperiment(DATA)
metaSeqObject_CSS = cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
CSS.readcounts=as.matrix(MRcounts(metaSeqObject_CSS,norm = TRUE,log = TRUE))


## GMPR
sizefactor.GMPR=GMPR(t(DATA))


# IMPUTATION --------------------------------------------------------------
imputation.CSS = mbImpute(otu_tab = CSS.readcounts)


# DA analysis -------------------------------------------------------------
# ALDEX




# ANCOM

