source('../../bin/RNASeqData.R')

#first get PNF data
tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,useCellNames=TRUE)
#ecounts.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,metric='est_counts',useCellNames=TRUE)

###collapse by protein-coding
pcvals<-grep('protein_coding',rownames(tpm.mat))
prot_coding.tpm=tpm.mat[pcvals,]

##then also collapse by gene
all.genes<-unique(sapply(rownames(tpm.mat),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))
pc.genes<-unique(sapply(rownames(prot_coding.tpm),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))


##now get the ccle data= this time from cindy

ccle.tpm<-read.table(synGet('syn5616092')@filePath,sep=',',row.names=1,header=T)
#ccle.ec<-read.table(synGet('syn5616077')@filePath,sep=',',row.names=1,header=T)
  

#now do basic matching
map<-read.table('../../data/gencodeGeneTranscriptMap.csv',sep=',',header=T)
gms=data.frame(Gene=pc.genes,Ens=sapply(as.character(map$target_id[match(pc.genes,map$gene)]),function(x) unlist(strsplit(unlist(strsplit(x,split='|',fixed=T))[2],split='.',fixed=T))[1]))

ccle.idx=match(gms[,2],rownames(ccle.tpm))

tpm.genes<-t(sapply(pc.genes,function(x) colSums(prot_coding.tpm[grep(paste(x,'.ENST',sep=''),rownames(prot_coding.tpm)),])))

##create larger matrix
all.tpms<-cbind(tpm.genes[which(!is.na(ccle.idx)),],ccle.tpm[ccle.idx[which(!is.na(ccle.idx))],])

#all.ecounts<-cbind(ecounts[which(!is.na(t.idx)),],ecounts.mat[t.idx[which(!is.na(t.idx))],])


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
doCorStats<-function(df,prefix='',minCor=0.25){
  #first compute correlation
    tcor=cor(df)
  ##select those stats that are of interest:
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))
  sel.vals<-tcor[mvals,]
  write.table(sel.vals,paste(prefix,'pearsonCors.tsv',sep=''),sep='\t')
  
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minPearsonCor',minCor,sep=''))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))
  sel.vals<-tcor[mvals,]
  
  write.table(sel.vals,paste(prefix,'spearmanCors.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minSpearmanCor',minCor,sep=''))
  
    #third: quanile normalization
  
 
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))

  sel.vals<-tcor[mvals,]
  write.table(sel.vals,paste(prefix,'qnormedPearson.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minQnormedPearsonCor',minCor,sep=''))              
  
}

##not looking great
doCorStats(all.tpms,'tpm')
#doCorStats(all.ecounts,'estCounts')
