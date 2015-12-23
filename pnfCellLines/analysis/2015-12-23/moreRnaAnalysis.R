source("../../bin/RNASeqData.R")

##first we have to store matrices

#run these once
#counts=rnaKallistoMatrix(TRUE,'est_counts')
#tpm=rnaKallistoMatrix(TRUE,'tpm')

#now run again
counts=rnaKallistoMatrix(FALSE,'est_counts')
tpm=rnaKallistoMatrix(FALSE,'tpm')
samp.mappings<-synTableQuery('SELECT "Sample Name","RNA-Seq Data","Sample Genotype" FROM syn5014742')@values


##now begin to do basic plotting
library(pheatmap)
library(ggbiplot)
plotMostVariableTranscripts<-function(count.mat,num=100,metric='tpm'){
  samp.names=samp.mappings$`Sample Name`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  samp.gen=samp.mappings$`Sample Genotype`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names
  var.vals<-apply(count.mat,1,var,na.rm=T)
  print(num)
  var.mat<-count.mat[order(var.vals,decreasing=T)[1:num],]
  pheatmap(var.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=samp.gen),filename=paste('top',num,'mostVariableGenesBy',metric,'.png',sep=''))
  
}

getGTAssociatedTranscripts<-function(count.mat,num=100,metric='tpm'){
  samp.names=samp.mappings$`Sample Name`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  samp.gen=samp.mappings$`Sample Genotype`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  names(samp.gen)<-samp.names
  orig.samp.gen<-samp.gen
  colnames(count.mat)<-samp.names
  ##make it abinary test
  samp.gen[which(samp.gen=='+-')]='++'
  ##now for each gene, evaluate significance of linear model predicting genotype
  lmvals<-apply(count.mat,1,function(x)
    summary(lm(Expr~Genotype,data=data.frame(Genotype=as.factor(samp.gen),Expr=x))))
  pvals<-sapply(lmvals,function(x){
    fs=x$fstatistic
    pf(fs[1],fs[2],fs[3],lower.tail=F)[1]
  })
  fdr.vals<-p.adjust(pvals,method='BH')
  #sig<-which(fdr.vals<0.15)
  #names(sig)<-sapply(names(sig),function(x) paste(unlist(strsplit(x,split='.',fixed=T))[1:2],collapse='.'))
  ##now plo those associated
  
  var.mat<-count.mat[order(pvals,decreasing=F)[1:num],]
  pheatmap(var.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=orig.samp.gen),filename=paste('top',num,'mostGTCorrelatedGenesBy',metric,'.png',sep=''))
  
}

plotNF1Transcripts<-function(count.mat,metric='tpm'){
  samp.names=samp.mappings$`Sample Name`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  samp.gen=samp.mappings$`Sample Genotype`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  names(samp.gen)<-samp.names
  orig.samp.gen<-samp.gen
  colnames(count.mat)<-samp.names
  nf1=grep("^NF1.ENST0",rownames(count.mat))
  print(paste("Found",length(nf1),'NF1 transcripts'))
  n.mat<-count.mat[nf1,]
  pheatmap(n.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=orig.samp.gen),filename=paste('NF1TranscriptsBy',metric,'.png',sep=''))
  
}

plotPCA<-function(count.mat,metric='tpm'){
  samp.names=samp.mappings$`Sample Name`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  samp.gen=samp.mappings$`Sample Genotype`[match(colnames(count.mat),samp.mappings$`RNA-Seq Data`)]
  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names
  
  zv<-which(apply(count.mat,1,var)==0)
  if(length(zv)>0)
    count.mat=count.mat[-zv,]
  pn<-prcomp(t(count.mat),center=T,scale=T)
  png(paste(metric,'valuesinRNASeqData.png',sep=''))
  p<-ggbiplot(pn,groups=samp.gen,var.axes=F)
  print(p)
  dev.off()
  
}


#plotNF1Transcripts(log2(tpm+0.01),'tpm')
#plotNF1Transcripts(log2(counts+0.01),'est_counts')

getGTAssociatedTranscripts(log2(tpm+0.01),100,'tpm')
getGTAssociatedTranscripts(log2(counts+0.01),100,'est_counts')

