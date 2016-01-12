source('../../bin/RNASeqData.R')

#first get PNF data
tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,useCellNames=TRUE)
ecounts.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,metric='est_counts',useCellNames=TRUE)

##now get the ccle data
synq=synapseQuery("select id,name from entity where parentId=='syn2325154'")

bcfiles<-synq[grep('bias_corrected.sf',synq[,1]),]

##now download all data
alldat<-lapply(bcfiles[1:200,2],function(x){
  print(paste('Getting/loading',x))
  tab<-read.table(synGet(x)@filePath)
  colnames(tab)<-c('Gene','Length','tpm','est_counts')
  return(tab)
})

fnames<-sapply(bcfiles[1:200,1],function(x) gsub('_quant_bias_corrected.sf','',fixed=T,x))

names(alldat)<-fnames
#extract TPMs and ecounts from sailfish files
tpms<-sapply(alldat,function(x){
  tpm=x[,3]
  names(tpm)<-x[,1]
  tpm
})

ecounts<-sapply(alldat,function(x){
  ecounts=x[,4]
  names(ecounts)<-x[,1]
  return(ecounts)
})

#get enst values 
tnames=sapply(rownames(tpm.mat),function(x) {
  arr=unlist(strsplit(x,split='.',fixed=T))
  arr[grep('ENST',arr)]})
  
#now do basic matching
t.idx=match(rownames(tpms),tnames)

##create larger matrix
all.tpms<-cbind(tpms[which(!is.na(t.idx)),],tpm.mat[t.idx[which(!is.na(t.idx))],])
all.ecounts<-cbind(ecounts[which(!is.na(t.idx)),],ecounts.mat[t.idx[which(!is.na(t.idx))],])


doCorPlots<-function(df,prefix=''){
  require(pheatmap)
    #first compute correlation
  tcor=cor(df)
  pheatmap(tcor,main=paste(prefix,'Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'pearsonCors.png',sep='_'))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.png',sep='_'))
  
  #third: quanile normalization
  
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.png',sep='_'))
  
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
  pheatmap(tcor,main=paste(prefix,'Q-normed Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'qnormPearsonCor.png',sep='_'))
}

##not looking great
doCorPlots(all.tpms)
doCorPlots(all.ecounts)
