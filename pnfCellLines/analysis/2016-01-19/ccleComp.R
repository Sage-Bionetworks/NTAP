source('../../bin/RNASeqData.R')

#first get PNF data
tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,useCellNames=TRUE)
ecounts.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,metric='est_counts',useCellNames=TRUE)

###collapse by protein-coding
pcvals<-grep('protein_coding',rownames(tpm.mat))
prot_coding.tpm=tpm.mat[pcvals,]

##then also collapse by gene
all.genes<-unique(sapply(rownames(tpm.mat),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))
pc.genes<-unique(sapply(rownames(prot_coding.tpm),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))


##now get the ccle data
synq=synapseQuery("select id,name from entity where parentId=='syn2325154'")

bcfiles<-synq[grep('bias_corrected.sf',synq[,1]),]

##now download all data
alldat<-lapply(bcfiles[,2],function(x){
  print(paste('Getting/loading',x))
  tab<-read.table(synGet(x)@filePath)
  colnames(tab)<-c('Gene','Length','tpm','est_counts')
  return(tab)
})

fnames<-sapply(bcfiles[,1],function(x) gsub('_quant_bias_corrected.sf','',fixed=T,x))

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
  pheatmap(tcor,main=paste(prefix,'Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'pearsonCors.pdf',sep='_'))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.pdf',sep='_'))
  
  #third: quanile normalization
  
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.pdf',sep='_'))
  
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
  pheatmap(tcor,main=paste(prefix,'Q-normed Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'qnormPearsonCor.pdf',sep='_'))
}



##add method to do some statistics on the correlations
doCorStats<-function(df,prefix='',minCor=0.6){
  #first compute correlation
    tcor=cor(df)
  ##select those stats that are of interest:
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))
  sel.vals<-tcor#tcor[mvals,]
  write.table(sel.vals,paste(prefix,'pearsonCors.tsv',sep=''),sep='\t')
  
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minPearsonCor',minCor,sep=''))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  #mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))
  sel.vals<-tcor#[mvals,]
  
  write.table(sel.vals,paste(prefix,'spearmanCors.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minSpearmanCor',minCor,sep=''))
  
    #third: quanile normalization
  
 
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
#  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))

  sel.vals<-tcor#[mvals,]
  write.table(sel.vals,paste(prefix,'qnormedPearson.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minQnormedPearsonCor',minCor,sep=''))              
  
}

##not looking great
doCorStats(all.tpms,'tpm')
doCorStats(all.ecounts,'estCounts')

for(file in list.files('.')[grep('tsv',list.files('.'))])
  synStore(File(file,parentId='syn5594111'),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-01-19/ccleComp.R',wasExecuted=TRUE)))

