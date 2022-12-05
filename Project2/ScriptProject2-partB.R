library(edgeR)

DATA<-read.table("raw_trascr_count.txt",sep="\t",row.names=1,header=TRUE)
indGroup<-c(grep("Group1",colnames(DATA)),grep("Group2",colnames(DATA)))
indCTRL<-grep("CTRL",colnames(DATA))
lg<-length(indGroup)
lc<-length(indCTRL)

target=as.data.frame(rep("Patients", dim(DATA)[2]), row.names=colnames(DATA))
colnames(target)=c("Groups")
target[indCTRL,]="CTRL"
group=factor(target$Groups)
design=model.matrix(~0+group)

#design = matrix(0, nrow=dim(DATA)[2], ncol = 2)
#colnames(design)=c("Group","CTRL")
#rownames(design)=colnames(DATA)
#design[indGroup1,1]=1
#design[indGroup2,1]=1
#design[indCTRL,2]=1

#target=matrix("",nrow=nrow(design), ncol=2)
#for(i in 1:nrow(design)){
#  if(design[i,1]==1){
#    o="Group"
#  }else{
#    o="CTRL"
#  }
#  target[i,1]=rownames(design)[i]
#  target[i,2]=o
#}


#TODO write design file
#fwrite(target,"design.txt",sep="\t", col.names=c("Sample","Group"), row.names=FALSE)


#Group <- factor(as.data.frame.matrix(target)$V2)
#target=model.matrix(~0+Group)
#print(Group)


y <- DGEList(counts=DATA)    # y is an object of type DGE
y <- calcNormFactors(y)   # This calculates the SF using the TMM normalization !!!
SF<-y$samples

y <- estimateGLMCommonDisp(y,design, verbose=TRUE) #phi common to the entire dataset
y <- estimateGLMTrendedDisp(y,design) #phi depends on mu
y <- estimateGLMTagwiseDisp(y,design) #phi is gene specific
fit <- glmFit(y,design) #finally the model fit (that accounts for raw NB data and scaling factors and seq. depth) 
summary(fit)

Confronti<-makeContrasts(Treatment ="groupCTRL-groupPatients",levels=design)
RES<-glmLRT(fit,contrast=Confronti[,"Treatment"])
# The first coloumn of RES reports the log_Fold_Change, i.e.: 
# log2(Normalized_data_average_GroupSARSCoV2 / Normalized_data_average_GroupMock)

alpha=0.05

selected=DATA[RES$table$PValue<alpha,]
selected.Group=sum(selected[,indGroup])
selected.CTRL=sum(selected[,indCTRL])
notselected=DATA[!RES$table$PValue<alpha,]
notselected.Group=sum(notselected[,indGroup])
notselected.CTRL=sum(notselected[,indCTRL])


salpha=selected.Group+selected.CTRL
EFP=min(sum(DATA)*0.8*alpha, salpha)
EFN=sum(DATA)-EFP-max(0,salpha-EFP)-min(sum(DATA)-salpha, sum(DATA)*0.8-EFP)
FDR=EFP/salpha
##################################################################################

DEbyEdgeR <- function(rawdat, groups, alpha=0.05){
  return(c(selected.genes.alpha, EFP))
}





alpha=0.05
rawdat = read.table("raw_trascr_count.txt")
groups= colnames(rawdat)
  
DEbyEdgeR(rawdat,groups)