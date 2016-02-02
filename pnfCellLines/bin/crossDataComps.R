##begin to try to correlate drug sensitivity data with existing data


source("../../bin/drugSensData.R")
source('../../bin/CNVData.R')
source('../../bin/RNASeqData.R')
library(ggplot2)

##first compare drug sensitivity data with RNA/CNV
drugRna<-function(valname='MAXR',gene='NF1',useGencode=F,doLog=F,
                  collapseAllCounts=FALSE,proteinCoding=FALSE,qthresh=0.1,
                  pthresh=0.001,doPlot=T){
  if(useGencode)
    rnaMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
  else
    rnaMat<-rnaKallistoMatrix(useCellNames=TRUE)
  
  ##get genotype
  gt.cell<-synTableQuery('SELECT "Sample Name","Sample Genotype" FROM syn5014742')@values
  
  drugMat<-getValueForAllCells(valname)
  
  ##now look for drug/rna correlations!
  cells=intersect(colnames(drugMat),colnames(rnaMat))
  
  rnaMat<-rnaMat[,cells]
  drugMat<-drugMat[,cells]
  
  rnavars<-which(apply(rnaMat,1,var)==0)
  nzvMat=rnaMat
  if(length(rnavars)>0)
    nzvMat=rnaMat[-rnavars,]
  nf1Mat=nzvMat[grep(paste("^",gene,".EN",sep=''),rownames(nzvMat)),]
  
  if(proteinCoding){
    pc=grep("protein_coding",rownames(nf1Mat))
    if(length(pc)>0){
      print(paste("Restricting",length(pc),"transcripts to only those that encode proteins"))
      
      nf1Mat=nf1Mat[pc,]
    }
  }
  if(collapseAllCounts && nrow(nf1Mat)>1){
    #allgenes<-unique(sapply(rownames(nf1Mat),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))
    nf1Mat<-t(colSums(nf1Mat))
    
  }
  
  if(nrow(nf1Mat)==0)
    return(NULL)
  if(doLog)
    nf1Mat=log2(nf1Mat+0.01)
  all.cors=apply(nf1Mat,1,function(x){
    apply(drugMat,1,function(y)
      cor(x,y,use='pairwise.complete.obs'))})
  
  require(psych)
  all.cor.ps=1-pnorm(fisherz(all.cors))
  ##all p-values, since sample size is so darn low
  #all.cor.ps=apply(nf1Mat,1,function(x){
  #  apply(drugMat,1,function(y)
  #    cor.test(x,y,use='pairwise.complete.obs')$p.value)})
  
  p.corrected<-apply(all.cor.ps,2,p.adjust,method='fdr')
 # names(p.corrected)<-names(all.cor.ps)
  sigs=which(p.corrected<qthresh,arr.ind=T)
  uncorrected=FALSE
  if(nrow(sigs)==0 && !is.na(pthresh)){
    uncorrected=TRUE
    sigs=which(all.cor.ps<pthresh,arr.ind=T)
  }
      if(nrow(sigs)==0)
    return(NULL)
  
  ##now can we plot those values? 
  dvals=unique(sigs[,1])
  print(paste('Found',length(dvals),'drugs that correlate with',gene,'expression'))
  drug<-c()
  exp<-c()
  drugname<-c()
  trans<-c()
  gt=c()
  for(i in 1:nrow(sigs)){
    r=sigs[i,]
    drug<-c(drug,drugMat[r[1],cells])
    gt=c(gt,gt.cell[match(cells,gt.cell[,1]),2])
    dn=rownames(drugMat)[r[1]]    
    if(is.na(dn))
      dn='NA_TARG'
  
    drugname=c(drugname,rep(dn,length(cells)))
    exp=c(exp,unlist(nf1Mat[r[2],cells]))
    tname=rownames(nf1Mat)[r[2]]
    if(is.null(tname))
      tname=gene
    trans=c(trans,rep(tname,length(cells)))
  }
  df=data.frame(DrugVals=drug,RNAExpr=exp,DrugName=drugname,Transcript=trans,Genotype=gt)
  #names(df)[1]=paste('Drug',valname,'Value',sep='')
  if(doPlot){
    png(paste(valname,'CorrelatedWith',ifelse(proteinCoding,'protCoding',''),gene,ifelse(collapseAllCounts,'byGene','byTrans'),ifelse(doLog,'log2',''),'expression.png',sep=''))
    p<-ggplot(df,aes(x=RNAExpr,y=DrugVals))+geom_point(aes(colour=DrugName,shape=Genotype),size=6,alpha=0.8)+geom_line(aes(colour=DrugName,linetype=Transcript))
    p<-p+ggtitle(paste(gene,'Expression for',valname,'\nvalues,',ifelse(uncorrected,paste('p<',pthresh),paste('q',qthresh))))
    print(p)
    dev.off()
  }
  return(df)
}

##now collapse by gene
##first compare drug sensitivity data with RNA/CNV
drugGene<-function(valname='MAXR',gene='NF1',useGencode=F){
  if(useGencode)
    rnaMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
  else
    rnaMat<-rnaKallistoMatrix(useCellNames=TRUE)
  drugMat<-getValueForAllCells(valname)
  
  ##now look for drug/rna correlations!
  cells=intersect(colnames(drugMat),colnames(rnaMat))
  rnavars<-which(apply(rnaMat[,cells],1,var)==0)
  nzvMat=rnaMat[-rnavars,]
  nf1Mat=nzvMat[grep(paste("^",gene,".EN",sep=''),rownames(nzvMat)),]
  geneSums<-colSums(nf1Mat[,cells])
  #all.cors=apply(nf1Mat[,cells],1,function(x){
  geneCor<-apply(drugMat[,cells],1,function(y)
      cor(geneSums,y,use='pairwise.complete.obs'))
  
  ##all p-values, since sample size is so darn low
  #all.cor.ps=apply(nf1Mat[,cells],1,function(x){
  genePs<-apply(drugMat[,cells],1,function(y)
      cor.test(geneSums,y,use='pairwise.complete.obs')$p.value)
  
  p.corrected<-p.adjust(genePs)#apply(all.cor.ps,2,p.adjust)
  
  sigs=which(p.corrected<0.1)  
  if(length(sigs)==0)
    return(NULL)
  ##now can we plot those values? 
  #dvals=unique(sigs[,1])
  print(paste('Found',length(sigs),'drugs that correlate with',gene,'expression'))
  drug<-c()
  exp<-c()
  drugname<-c()
  trans<-c()
  for(r in sigs){
    #
    drug<-c(drug,drugMat[r,cells])
    drugname=c(drugname,rep(rownames(drugMat)[r],length(cells)))
    exp=c(exp,geneSums)#unlist(nf1Mat[r[2],cells]))
    trans=c(trans,rep(gene,length(cells)))
  }
  df=data.frame(DrugVals=drug,RNAExpr=exp,DrugName=drugname,Transcript=trans)
  #names(df)[1]=paste('Drug',valname,'Value',sep='')
  png(paste(valname,'CorrelatedWith',gene,'Geneexpression.png',sep=''))
  p<-ggplot(df,aes(x=RNAExpr,y=DrugVals))+geom_point(aes(colour=DrugName,shape=Transcript))+geom_line(aes(colour=DrugName,shape=Transcript))
  p<-p+ggtitle(paste(gene,'Expression for',valname,'values, corrected p<0.1'))
  print(p)
  dev.off()
}