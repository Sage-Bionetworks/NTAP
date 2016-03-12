##test out PC analysis before copying to crossDataComp.R file

library(ggbiplot)
library(mHG)
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


getPCAEnrichment<-function(pca,names,enrichList){
  require(mHG)
  #for each element in the loadings, rank and compute enrichment of all functions
  res<-lapply(colnames(pca$rotation),function(x){
    pcname=x
    x=pca$rotation[,x]
    ordered.list=names(sort(abs(x),decreasing=F))
    pvals=lapply(enrichList,function(y){
          vec=rep(0,length(ordered.list))
          vec[which(ordered.list%in%unlist(geneIds(y)))]<-1
          t <- mHG.test(vec)
          t$p.value
    })
    pvals=unlist(pvals)
    return(data.frame(Component=rep(pcname,length(pvals)),
                      Target=unlist(lapply(enrichList,setName)),
                      Drugs=unlist(lapply(enrichList,function(x) paste(geneIds(x),collapse=','))),
                      Pvalue=pvals,AdjustedP=p.adjust(pvals,method='BH')))
  })
  df=do.call('rbind',res)
  write.table(df,paste(names,'pcaLoadingRankedminHGEnrichment.csv',sep='_'),sep=',')
  return(df)
  
}

genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)

ncatsMat<-getValueForAllCells("FAUC")
ncatsMat[which(is.na(ncatsMat),arr.ind=TRUE)]<-0.0
ncatsReMat=getRecalculatedAUCMatrix()

ncats.cells=intersect(colnames(genCodeMat),colnames(ncatsMat))
ncats.genotype=dfiles$entity.sampleGenotype[match(ncats.cells,dfiles$entity.sampleName)]
names(ncats.genotype)<-ncats.cells

rna.pc=do.pc(genCodeMat[,ncats.cells])


library(GSVA)
library(GSVAdata)
data(c2BroadSets)

eids<-as.data.frame(fread('../../data/HugoGIDsToEntrez_DAVID.txt',sep='\t',header=T))

eidx=match(rownames(genCodeMat),eids$From)
entrez.mat=genCodeMat[,ncats.cells]
entrez.mat<-entrez.mat[which(!is.na(eidx)),]
rownames(entrez.mat)<-eids[eidx[which(!is.na(eidx))],'To']

##now do pc of entrez mat
entrez.rna.pc=do.pc(entrez.mat)
pc.gene.enrich=getPCAEnrichment(entrez.rna.pc,'geneEpressionBroadEnrichment',c2BroadSets)



ncats.pc=do.pc(ncatsMat[,ncats.cells])
test.targ.list<-ncatsDrugTargets()

##create a similar list...
ncats.targ.list=GeneSetCollection(lapply(as.character(unique(test.targ.list$Target)), 
                                         function(x) {
                                           drugs=as.character(unique(test.targ.list[which(test.targ.list$Target==x),'Drug']))
                                           gs=GeneSet()
                                           geneIds(gs)<-drugs
                                           setName(gs)<-x
                                           return(gs)}))



ncats.targ.enrich<-getPCAEnrichment(ncats.pc,'ncatsDrugTargetEnrich',ncats.targ.list)

##now do ncats rescored
ncatsReMat[which(is.na(ncatsReMat),arr.ind=T)]<-0.0
ncatsr.pc=do.pc(ncatsReMat[,ncats.cells])
ncats.targ.enrich<-getPCAEnrichment(ncats.pc,'rescoredNcatsDrugTargetEnrich',ncats.targ.list)


synid='syn5685756'

#but also compare to ctrp
#source("../../bin/ctrpSingleAgentScreens.R")
#source("../../bin/ccleData.R")

#ctrpMat=as.data.frame(getCtrpScreensAsMatrix())
#ctrpMat[which(is.na(ctrpMat),arr.ind=T)]<-0.0

#ccle.tpm<-getCCLEDataTPM(removeDupes=TRUE)

#com.cells=intersect(colnames(ctrpMat),colnames(ccle.tpm))


#ccle.pc=do.pc(ccle.tpm[,com.cells])
#ctrp.pc=do.pc(ctrpMat[,com.cells])
#
#sapply(colnames(ccle.pc$rotation)[1:10],function(x){
#  write.table(names(sort(ccle.pc$rotation[,x]^2)),file=paste('ccleGeneExpression',x,'.txt',sep=''),row.names=F,col.names=F,quote=F)
#})

#ctrpr.pc=do.pc(ctrpReMat[,com.cells])

this.script="https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-24/drug_trans_pca_analysis.R"


##what if we combine the tpm values? 
#com.genes<-intersect(rownames(genCodeMat),rownames(ccle.tpm))
#com.mat<-cbind(genCodeMat[com.genes,],ccle.tpm[com.genes,])
#require(limma)
#origin<-c(rep("PNF",ncol(genCodeMat)),rep("CCLE",ncol(ccle.tpm)))

#design<-model.matrix(~factor(origin))
#colnames(design)<-c("CCLE","PNF")
#fit=lmFit(com.mat,design)
#fit <- eBayes(fit)


#
#names(res)<-unique(test.targ.list$Target)
#eid.match=eids[match(rownames(com.mat),eids[,1]),2]
#eid.rows=which(!is.na(eid.match))
#vnames=eid.match[eid.rows]

#ent.rna<-com.mat[eid.rows,]
#ent.rna<-ent.rna[-which(duplicated(vnames)),]
#rownames(ent.rna)=vnames[-which(duplicated(vnames))]

#es<-gsva(as.matrix(ent.rna),c2BroadSets,rnaseq=T,no.bootstraps=10)



