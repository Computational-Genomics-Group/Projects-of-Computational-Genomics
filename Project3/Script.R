# PACKAGES ----------------------------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE))  install.packages("devtools")
library(devtools)

if(!require(compositions)) install.packages(compositions)
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
if (!requireNamespace("GUniFrac", quietly = TRUE))  install.packages("GUniFrac")
if (!requireNamespace("glmnet", quietly = TRUE))  install.packages("glmnet")
if (!requireNamespace("Matrix", quietly = TRUE))  install.packages("Matrix")
if (!requireNamespace("mbImpute", quietly = TRUE))  install_github("ruochenj/mbImpute/mbImpute R package")

library(GUniFrac)
library(metagenomeSeq)
library(glmnet)
library(Matrix)
library(compositions)
library(mbImpute)
library(devtools)

# ReadData ----------------------------------------------------------------
filename="Genus_otu_table.txt"
DATA=as.matrix(read.table(filename,header = TRUE, sep = "\t"))
#DATA = DATA[1:4, 1:3]

# DATA TRANSFORMATIONS ----------------------------------------------------
data.alr=as.matrix(alr(DATA));
data.clr=as.matrix(clr(DATA));
data.ilr=as.matrix(ilr(DATA));


# OUR DATA TRANSFORMATIONS
# Clr transformation with pseudocounts
data.pseudocounts=DATA+1
np<-dim(data.pseudocounts)[2] #number of samples
nt<-dim(data.pseudocounts)[1] #number of taxa
clrtransform_pseudo<-data.pseudocounts #just to initialize clrtransform_pseudo
for (i in (1:np)) {
  den<-(prod(data.pseudocounts[,i]))^(1/nt) #geometric mean of column i
  clrtransform_pseudo[,i]<-log2(data.pseudocounts[,i]/den) #clr transformation of column i
}
head(clrtransform_pseudo[,1:10])


## CLR transformation without pseudocounts
np<-dim(DATA)[2] #number of samples
nt<-dim(DATA)[1] #number of taxa
clrtransform<-DATA #just to initialize clrtransform
for (i in (1:np)) {
  x<-DATA[,i]
  x[which(x==0)]<-NA
  den<-(prod(x,na.rm=TRUE)^(1/length(which(!is.na(x))))) #geometric mean of column i (excluding 0)
  clrtransform[,i]<-log2(x/den) #clr transformation of column i
  clrtransform[which(is.na(x)),i]<-0
}
head(DATA)[1:5,1:10]
head(clrtransform)[1:5,1:10]


# NORMALIZATION -----------------------------------------------------------
## CSS
data.alr=as.matrix(alr(DATA));
data.clr=as.matrix(clr(DATA));
data.ilr=as.matrix(ilr(DATA));
metaSeqObject_CSS = newMRexperiment(DATA)
metaSeqObject_CSS = cumNorm(metaSeqObject_CSS, p=cumNormStat(metaSeqObject_CSS))
data_CSS=as.matrix(MRcounts(metaSeqObject_CSS,norm = FALSE,log = TRUE))

metaSeqObject_CSS.clr = newMRexperiment(data.clr)
metaSeqObject_CSS.clr = cumNorm(metaSeqObject_CSS.clr, p=cumNormStat(metaSeqObject_CSS.clr))
data_CSS.clr=as.matrix(MRcounts(metaSeqObject_CSS.clr,norm = FALSE,log = TRUE))

metaSeqObject_CSS.alr = newMRexperiment(data.alr)
metaSeqObject_CSS.alr = cumNorm(metaSeqObject_CSS.alr, p=cumNormStat(metaSeqObject_CSS.alr))
data_CSS.alr=as.matrix(MRcounts(metaSeqObject_CSS.alr,norm = TRUE,log = TRUE))

metaSeqObject_CSS.ilr = newMRexperiment(data.ilr)
metaSeqObject_CSS.ilr = cumNorm(metaSeqObject_CSS.ilr, p=cumNormStat(metaSeqObject_CSS.ilr))
data_CSS.ilr=as.matrix(MRcounts(metaSeqObject_CSS.ilr,norm = TRUE,log = TRUE))

head(DATA[,1:10])
head(data_CSS[,1:10])
head(data_CSS.alr[,1:10])
head(data_CSS.clr[,1:10])
head(data_CSS.ilr[,1:10])

## GMPR
GMPR_factors<- GMPR(OTUmatrix = t(DATA), min_ct = 5, intersect_no = 4) #see help for parameters meaning
data.GMPR<- t(t(DATA)/GMPR_factors)

GMPR_factors[1:10]
head(feature_table_gen)[1:5,1:10]
head(feature_table_gen_gmpr)[1:5,1:10]

# IMPUTATION --------------------------------------------------------------
label_samples<-read.table("metadata_table.txt",sep="\t",header=T)
imputation.CSS = mbImpute(condition = label_samples$DiseaseState, otu_tab = data.CSS, parallel = TRUE,metadata = label_samples, ncores=10, unnormalized = TRUE)
imputation.GMPR= mbImpute(condition = label_samples$DiseaseState, otu_tab = DATA, parallel = TRUE, metadata = label_samples, ncores=10, unnormalized = TRUE)

plotly
# DA analysis -------------------------------------------------------------
# ALDEX
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")
library(ALDEx2)
aldex_results<- aldex(reads = feature_table_gen, conditions =label_samples$DiseaseState, mc.samples = 128,test = "t")


# ANCOM

