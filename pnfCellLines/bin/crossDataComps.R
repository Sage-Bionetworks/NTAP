##begin to try to correlate drug sensitivity data with existing data


source("../../bin/drugSensData.R")
#source('../../bin/CNVData.R')
source('../../bin/RNASeqData.R')
source("../../bin/singleDrugAnalysis.R")

library(ggplot2)

drugPlsr<-function(valname='MAXR',useGencode=T,collapseAllCounts=TRUE,proteinCoding=TRUE){
  require(pls)
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
  nzvMat=rnaMat[-rnavars,]
  if(collapseAllCounts)
    nzvMat<-nzvMat[grep('protein_coding',rownames(nzvMat)),]
  df=data.frame(t(rbind(drug=x,nzvMat)))  
  res<-plsr(drug~.,df)
}

#'Basic correlation analysis comparing individual gene expression with
#'drug value of interest (e.g. MAXR, FAUC)
#'@param
#'@param
#'@param
#'@param
#'@param
#'@param
#'@param
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



#'Generic function to compute RNA/Drug normalization
#'@param drugMat - matrix of drug AUC values by cell line
#'@param rnaMat - matrix of rna expressionv values by cell line
#'@param prefix - to name files
#'@param sampleCombos - if sampling from all possible combos, provide number here
#'@alpha alpha parameter to use as double sigma input
#'@return list of files summarizing results
computeDrugRnaNormalizedCor<-function(drugMat,rnaMat,prefix='',sampleCombos=NA,alpha=2.5,doPlot=TRUE){
# if(require(parallel))
  #  lapply<-function(x,...) mclapply(x,...)

  #collect overlap, reshape matrices
  cells=intersect(colnames(drugMat),colnames(rnaMat))
  print(paste('found',length(cells),'cells common between both drug and expression experiments'))
  rnaMat<-rnaMat[,cells]
  drugMat<-drugMat[,cells]

  ##remove zero variance genes/drugs
  zv=which(apply(rnaMat,1,var)==0)
  if(length(zv)>0)
    rnaMat<-rnaMat[-zv,]
  zv=which(apply(drugMat,1,var)==0)
  if(length(zv)>0)
    drugMat<-drugMat[-zv,]

  ##enumerate all combos (or sample for later)
  gene.drug.combos=do.call('rbind',lapply(rownames(rnaMat),function(x) cbind(Gene=rep(x,nrow(drugMat)),Drug=rownames(drugMat))))
  #add in sampling option to estimate distribution to avoid computing all

  if(!is.na(sampleCombos)){
    if(sampleCombos>nrow(gene.drug.combos)){
	print('Sample size is greater than all permutations, just doing all')
        idx=1:nrow(gene.drug.combos)
    }else{
    idx=sample(1:nrow(gene.drug.combos),sampleCombos)
  }}
  else
    idx=1:nrow(gene.drug.combos)
  #make into list to multi-core
  gdc<-lapply(idx,function(x) gene.drug.combos[x,])

  ##get all combos
  res=lapply(gdc,function(x,rnaMat,drugMat){
    #print(paste(x,collapse=','))
    ret=list(Pearson=NA,fisherZ=NA,Overlap=NA)

    g=as.numeric(rnaMat[x[[1]],])
    d=as.numeric(drugMat[x[[2]],])

    na1=intersect(which(!is.na(g)),which(!is.nan(g)))
    na2=intersect(which(!is.na(d)),which(!is.nan(d)))
    try(ret$Overlap<-length(intersect(na1,na2)))

    try(ret$Pearson<-stats::cor(g,d,use='p'))
    try(ret$fisherZ<-fzCor(g,d))
    return(ret)
  },rnaMat,drugMat)

  full.res=data.frame(do.call('rbind',gdc),do.call('rbind',res))

  full.res$fisherZ<-as.numeric(full.res$fisherZ)
  full.res$Overlap<-as.numeric(full.res$Overlap)
  full.res$Pearson<-as.numeric(full.res$Pearson)

  full.res$doubleSigma=dSigTransform(as.numeric(full.res$fisherZ),alpha=alpha)

  if(doPlot){
      require(ggplot2)

  ##plot

      p<-ggplot(full.res)+geom_point(mapping=aes(x=Pearson,y=doubleSigma,colour=Overlap),stat='identity')
      p<-p+scale_colour_gradientn(colours=rainbow(5))+ggtitle(paste('Pearson R vs. transformed (alpha=',alpha,') \nfor',ifelse(is.na(sampleCombos),'all',sampleCombos),prefix,'\nDrug-Transcript Combinations'))

      basename=paste(prefix,'alpha',alpha,'pearsonVsDoubleSig',ifelse(is.na(sampleCombos),'all',sampleCombos),'Pairs',sep='_')

      pngname=paste(basename,'pdf',sep='.')
      pdf(pngname)
      print(p)
      dev.off()
  }
  #write to file
  tabname=paste(basename,'tab',sep='.')
  write.table(full.res,file=tabname,row.names=F,col.names=T,sep='\t')
  return(paste(basename,c('pdf','tab'),sep='.'))

}


##now collapse by gene
##first compare drug sensitivity data with RNA/CNV
drugGene<-function(valname='MAXR',gene='NF1',useGencode=F,by.pval=TRUE){
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
  print(geneSums)
  #all.cors=apply(nf1Mat[,cells],1,function(x){
  geneCor<-apply(drugMat[,cells],1,function(y)
      cor(geneSums,y,use='pairwise.complete.obs'))

  ##all p-values, since sample size is so darn low
  #all.cor.ps=apply(nf1Mat[,cells],1,function(x){
  genePs<-apply(drugMat[,cells],1,function(y){
      p<-1.0
      try(p<-cor.test(geneSums,y,use='pairwise.complete.obs')$p.value)
      return(p)})

  p.corrected<-p.adjust(genePs)#apply(all.cor.ps,2,p.adjust)

  if(by.pval){
    sigs=which(p.corrected<0.1)
    if(length(sigs)==0)
      sigs=which(genePs<0.01)
  }
  else{
    sigs=which((rank(abs(geneCor))/length(geneCor))>0.99)
  }
  ##now can we plot those values?
  #dvals=unique(sigs[,1])
  print(paste('Found',length(sigs),'drugs that correlate with',gene,'expression'))
  drug<-c()
  exp<-c()
  drugname<-c()
  trans<-c()
  for(r in sigs){
    drug<-c(drug,drugMat[r,cells])
    drugname=c(drugname,rep(rownames(drugMat)[r],length(cells)))
    exp=c(exp,geneSums)#unlist(nf1Mat[r[2],cells]))
    trans=c(trans,rep(gene,length(cells)))
  }
  df=data.frame(DrugVals=drug,RNAExpr=exp,DrugName=drugname,Transcript=trans)
  #names(df)[1]=paste('Drug',valname,'Value',sep='')
  png(paste(valname,'CorrelatedWith',gene,'Geneexpression',ifelse(by.pval,'','top1Percent'),'.png',sep=''))
  p<-ggplot(df,aes(x=RNAExpr,y=DrugVals))+geom_point(aes(colour=DrugName,shape=Transcript))+geom_line(aes(colour=DrugName,shape=Transcript))
  p<-p+ggtitle(paste(gene,'Expression for',valname,'values, ',ifelse(by.pval,'corrected p<0.1','top1%')))
  print(p)
  dev.off()
  return(df)
}


ds.dist<-read.table(synGet('syn5705377')@filePath,header=T)
qvals <-quantile(ds.dist$doubleSigma,c(0.005,0.995),na.rm=T)
z.vals<-quantile(as.numeric(as.character(ds.dist$fisherZ)),c(0.005,0.995),na.rm=T)
##now collapse by gene
##first compare drug sensitivity data with RNA/CNV
drugGeneNorm<-function(valname='MAXR',gene='NF1',useGencode=F,by.pval=TRUE){
  if(useGencode)
    rnaMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
  else
    rnaMat<-rnaKallistoMatrix(useCellNames=TRUE)
  drugMat<-getValueForAllCells(valname)
  targs<-ncatsDrugTargets()
    

  ##now look for drug/rna correlations!
  cells=intersect(colnames(drugMat),colnames(rnaMat))
  rnavars<-which(apply(rnaMat[,cells],1,var)==0)
  nzvMat=rnaMat[-rnavars,]
  gvals<-grep(paste("^",gene,".EN",sep=''),rownames(nzvMat))
  if(length(gvals)==0)
    return(data.frame(DrugVals=NA,RNAExpr=NA,DrugName=NA,Transcript=NA))

  nf1Mat=nzvMat[gvals,]
  geneSums<-colSums(nf1Mat[,cells])
 # print(geneSums)
  #all.cors=apply(nf1Mat[,cells],1,function(x){
  geneCor<-data.frame(do.call('rbind',lapply(rownames(drugMat),function(y){
        c(Drug=y,Gene=gene,
          R=cor(geneSums,drugMat[y,cells],use='pairwise.complete.obs'),
                    FZ=fzCor(geneSums,drugMat[y,cells]))})))
        
  geneCor$dsTrans=dSigTransform(as.numeric(as.character(geneCor$FZ)),alpha=1)
  ##all p-values, since sample size is so darn low
  #all.cor.ps=apply(nf1Mat[,cells],1,function(x){
  navals=which(is.na(geneCor$dsTrans))
  if(length(navals)>0)
    geneCor=geneCor[-navals,]
  
  ##first filter by quantile
  sigs=c(which(geneCor$dsTrans<qvals[1]),which(geneCor$dsTrans>qvals[2]))
  sigs1=c(which(as.numeric(as.character(geneCor$FZ))<z.vals[1]),
         which(as.numeric(as.character(geneCor$FZ))>z.vals[2]))
  sigs2=which((rank(abs(geneCor$dsTrans))/nrow(geneCor))>0.995)
  
  if(length(sigs1)>10)
    sigs<-intersect(sigs1,sigs2)
  else
    sigs<-sigs1
  ##now can we plot those values?
  #dvals=unique(sigs[,1])
  #print(paste('Found',length(sigs),'drugs that correlate with',gene,'expression'))
  drug<-c()
  exp<-c()
  drugname<-c()
  trans<-c()
  for(r in sigs){
    drug<-c(drug,drugMat[r,cells])
    drugname=c(drugname,rep(paste(rownames(drugMat)[r],' (',targs$Target[match(rownames(drugMat)[r],targs$Drug)],')',sep=''),length(cells)))
    exp=c(exp,geneSums)#unlist(nf1Mat[r[2],cells]))
    trans=c(trans,rep(gene,length(cells)))
  }
  df=data.frame(DrugVals=drug,RNAExpr=exp,DrugName=drugname,Transcript=trans)
  #names(df)[1]=paste('Drug',valname,'Value',sep='')
  png(paste(valname,'CorrelatedWith',gene,'Geneexpression',ifelse(by.pval,'','top1Percent'),'.png',sep=''))
  p<-ggplot(df,aes(x=RNAExpr,y=DrugVals))+geom_point(aes(colour=DrugName,shape=Transcript))+geom_line(aes(colour=DrugName,shape=Transcript))
  p<-p+ggtitle(paste(gene,'Expression for',valname,'values, ',ifelse(by.pval,'corrected p<0.1','top1%')))
  print(p)
  dev.off()
  return(df)
}


#'
#'Lastly to do a global correlation of all transcripts/genes across all drugs
#'we first need to test normalization
testNormalization<-function(geneMat,drugMat,alpha=c(0.1,1,10,100),prefix=''){

}

#'
#'Now do large correlation matrix
globalCorrelationAnalysis<-function(geneMat,drugMat,alpha,prefix=''){}
