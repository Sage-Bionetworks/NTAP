##test out PC analysis before copying to crossDataComp.R file

library(ggbiplot)
##let's analyze with NCATS
source("../../bin/ncatsSingleAgentScreens.R")
source("../../bin/RNASeqData.R")

#' quick function to do PC analysis
do.pc<-function(mat){
  zr=which(apply(mat,1,var,na.rm=T)==0)
  zc=which(apply(mat,2,var,na.rm=T)==0)
  if(length(zr)>0)
    mat=mat[-zr,]
  if(length(zc)>0)
    mat=mat[,-zc]
  pcp=prcomp(t(mat),scale=T,center=T)

  pcp
}


computePCAcor<-function(pca1,pca2,names1,names2){
  min.pc=min(10,min(ncol(pca1$x),ncol(pca2$x)))
  
  sd1=pca1$sdev[1:min.pc]
  sd2=pca2$sdev[1:min.pc]
  names(sd1)<-names(sd2)<-colnames(pca1$x)[1:min.pc]
  cmat=cor(pca1$x[,1:min.pc],pca2$x[,1:min.pc])
  pheatmap(cmat,cluster_rows = F, cluster_cols = F,
           annotation_row=data.frame(Variance=sd1/sum(sd1)), labels_row=paste(names1,names(sd1)),
           annotation_col=data.frame(Variance=sd2/sum(sd2)),labels_col=paste(names2,names(sd2)),
           filename=paste(names1,'vs',names2,'PCcorrelations.png',sep='_'))
  return(cmat)
}


genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)

ncatsMat<-getValueForAllCells("FAUC")
ncatsMat[which(is.na(ncatsMat),arr.ind=TRUE)]<-0.0
ncatsReMat=getRecalculatedAUCMatrix()

ncats.cells=intersect(colnames(genCodeMat),colnames(ncatsMat))
ncats.genotype=dfiles$entity.sampleGenotype[match(ncats.cells,dfiles$entity.sampleName)]
names(ncats.genotype)<-ncats.cells

rna.pc=do.pc(genCodeMat[,ncats.cells])
ncats.pc=do.pc(ncatsMat[,ncats.cells])
computePCAcor(rna.pc,ncats.pc,'RNA','NCATS')

png('PCA_plots_of_RNASeq.png')
p1=ggbiplot(rna.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=1:2,groups=ncats.genotype)+ggtitle('RNA PC1,2')
print(p1)
dev.off()
png("PCA_plots_of_NCATS.png")
p2=ggbiplot(ncats.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=2:1,groups=ncats.genotype)+ggtitle('NCATS PC2,1')
print(p2)
dev.off()

##now do ncats rescored
ncatsReMat[which(is.na(ncatsReMat),arr.ind=T)]<-0.0
ncatsr.pc=do.pc(ncatsReMat[,ncats.cells])

computePCAcor(rna.pc,ncatsr.pc,'RNA','NCATS_rescored')

png("PCA_plots_of_NCATS_rescored.png")
p2=ggbiplot(ncatsr.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=2:1,groups=ncats.genotype)+ggtitle('NCATS Rescored PC2,1')
print(p2)
dev.off()

#but also compare to ctrp
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ccleData.R")

ctrpMat=as.data.frame(getCtrpScreensAsMatrix())
ctrpMat[which(is.na(ctrpMat),arr.ind=T)]<-0.0

ccle.tpm<-getCCLEDataTPM(removeDupes=TRUE)

com.cells=intersect(colnames(ctrpMat),colnames(ccle.tpm))


ccle.pc=do.pc(ccle.tpm[,com.cells])
ctrp.pc=do.pc(ctrpMat[,com.cells])

computePCAcor(ccle.pc,ctrp.pc,'CCLE','CTRP')

#now do biplots

#also get re-normalized CTRP
ctrpReMat<-ctrpDoseResponseCurve(FALSE,TRUE)
ctrpReMat[which(is.na(ctrpReMat),arr.ind=T)]<-0.0
ctrpr.pc=do.pc(ctrpReMat[,com.cells])
computePCAcor(ccle.pc,ctrpr.pc,'CCLE','CTRP_rescored')


png('PCA_plots_of_CCLE_RNASeq.png')
p1=ggbiplot(ccle.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=1:2)+ggtitle('CCLE PC1,2')
print(p1)
dev.off()

png('PCA_plots_of_CTRP.png')
p1=ggbiplot(ctrp.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=c(4,1))+ggtitle('CCLE PC4,1')
print(p1)
dev.off()

png('PCA_plots_of_CTRP_rescored.png')
p1=ggbiplot(ctrpr.pc,var.axes=F,scale = 0,labels=ncats.cells,choices=c(4,1))+ggtitle('CCLE PC4,1')
print(p1)
dev.off()



