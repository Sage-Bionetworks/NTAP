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

#' Compute correlation between two PCA options derived from prcomp
#' @param pca1 first PC object
#' @param pca2
#' @param names1 descriptor of first object
#' @param names2 descriptor of secon object
#' @return matrix of all correlations
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

#'Performs correlation testing between PCA and an arbitrary matrix
#'assuming both have the same number of columns 
computePCAMatCor<-function(pca1,mat2,names1,names2){
  #first just get the variation
  sd1=pca1$sdev
  names(sd1)<-colnames(pca1$x)
  
  over=intersect(rownames(pca1$x),colnames(mat2))
  #now compute the correlations
  cmat=cor(pca1$x[over,],t(mat2[,over]))
  cmat[which(is.na(cmat))]<-0.0
  
  #plot correlations of first two components
  high.cor=union(which(abs(cmat[1,])>0.8),which(abs(cmat[2,])>0.8))
  plot(cmat[1,],cmat[2,],
       main=paste(names1,'PCs 1 vs 2 correlated with',names2),
       xlab='PC1',ylab='PC2')
  
  return(cmat)
}

#we want to correlate drug response with gene expression
genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)
genePathMat <- read.table(synGet('syn5689231')@filePath)
ncatsMat<-getValueForAllCells("FAUC")
ncatsMat[which(is.na(ncatsMat),arr.ind=TRUE)]<-0.0
ncatsReMat=getRecalculatedAUCMatrix()

ncats.cells=intersect(colnames(genCodeMat),colnames(ncatsMat))
ncats.genotype=dfiles$entity.sampleGenotype[match(ncats.cells,dfiles$entity.sampleName)]
names(ncats.genotype)<-ncats.cells

rna.pc=do.pc(genCodeMat[,ncats.cells])
ncats.pc=do.pc(ncatsMat[,ncats.cells])

targs<-ncatsDrugTargets()
tvals=as.character(targs$Target)
names(tvals)<-targs$Drug

c1=computePCAMatCor(rna.pc,ncatsMat,'RNA','NCATS')

sset=colnames(c1)[grep('inib',colnames(c1))]

pheatmap(c1[,sset],annotation_col = data.frame(Target=tvals[sset]),cluster_rows = F,filename='RNAPcsVsDrugSubset.png')

c2=computePCAMatCor(ncats.pc,genCodeMat,'NCATS','RNA')
dt=which(colnames(c2)%in%targs$Target)
tcounts=as.numeric(table(targs$Target))
names(tcounts)<-names(table(targs$Target))
dt=intersect(names(sort(tcounts)),colnames(c2))

pheatmap(t(c2[,dt]),cluster_cols = F,cluster_rows = F,cellwidth=10,cellheight=10,
         annotation_row = data.frame(NumDrugs=tcounts),
         filename='drugTargetsVsDrugPCs.png')

ddt=intersect(names(sort(tcounts[which(tcounts>5)])),colnames(c2))
pheatmap(t(c2[,ddt]),cluster_cols = F,cluster_rows = F,cellwidth=10,cellheight=10,
         annotation_row = data.frame(NumDrugs=tcounts),
         filename='fiveOrMoredrugTargetsVsDrugPCs.png')

##now do ncats rescored
ncatsReMat[which(is.na(ncatsReMat),arr.ind=T)]<-0.0
ncatsr.pc=do.pc(ncatsReMat[,ncats.cells])

#computePCAGenecor(rna.pc,ncatsr.pc,'RNA','NCATS_rescored')

#now do biplots


this.script="https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-25/drug_trans_pca_analysis.R"

for(file in list.files('.'))
  if(length(grep('png',file))>0){
#    f=File(file,parentId='')
#    synStore(f,activityName='PCA correlation',
#             used=list(list(url=this.script,wasExecuted=TRUE)))
  }
    


