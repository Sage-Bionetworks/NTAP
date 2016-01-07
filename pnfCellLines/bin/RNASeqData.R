###
###
##
## RNASeqData.R
## Basic set of library functions designed to collect and manage RNA-Seq data for NTAP PNF cell lines
##
###

rna.dir<-'syn5562003'
grch38.dir<-'syn5579783'
gencode.dir='syn5579785'
require(synapseClient)
require(data.table)
synapseLogin()

samp.mappings<-synTableQuery('SELECT "Sample Name","RNA-Seq Data","RNA-Seq Data (Gencode)", "Sample Genotype" FROM syn5014742')@values


rnaKallistoFiles<-function(){
  rfiles<-synapseQuery(paste('select * from entity where parentId=="',grch38.dir,'"',sep=''))
  rfiles<-rfiles[which(!is.na(rfiles$entity.sampleName)),]
  all.files<-lapply(rfiles$entity.id,function(x)
    as.data.frame(fread(synGet(x)@filePath))
  )
  names(all.files)<-rfiles$entity.id
  return(all.files)
}

rnaGencodeKallistoFiles<-function(){
  rfiles<-synapseQuery(paste('select * from entity where parentId=="',gencode.dir,'"',sep=''))
  rfiles<-rfiles[which(!is.na(rfiles$entity.sampleName)),]
  all.files<-lapply(rfiles$entity.id,function(x)
    as.data.frame(fread(synGet(x)@filePath))
  )
  names(all.files)<-rfiles$entity.id
  return(all.files)
}

rnaAnnotations<-function(){
  rfiles<-synapseQuery(paste('select id,sampleName from entity where parentId=="',grch38.dir,'"',sep=''))
  rfiles<-rfiles[which(!is.na(rfiles$entity.sampleName)),]
  samps<-rfiles$entity.sampleName
  names(samps)<-rfiles$entity.id
  return(samps)
}

rnaKallistoMatrix<-function(buildFromFiles=FALSE,metric='tpm',useCellNames=FALSE){
  ##metric is either 'tpm' or 'est_counts'
  if(!metric%in%c('tpm','est_counts')){
    print(paste(metric,'is not a valid kallisto output'))
    return(NULL)
  }
  if(buildFromFiles){
        allfiles<-rnaKallistoFiles()
        all.quants=sapply(allfiles,function(x){
          quants<-x[,metric]
          names(quants)=apply(x[,1:2],1,paste,collapse='.')
          return(quants)
        })
        fname=paste('kallistoDerived',metric,'RNASeq_values.tsv',sep='')
        write.table(all.quants,file=fname,sep='\t')
        sf=File(fname,parentId=grch38.dir)
        synStore(sf,used='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/RNASeqData.R')
        tab=all.quants

    }else{
      if(metric=='est_counts')
        tab=read.table(synGet('syn5562376')@filePath)
      else if(metric=='tpm')
        tab=read.table(synGet('syn5562378')@filePath)

    }
  if(useCellNames){
    annotes=rnaAnnotations()
    colnames(tab)<-annotes[colnames(tab)]
  }
  return(tab)
}

rnaGencodeKallistoMatrix<-function(buildFromFiles=FALSE,metric='tpm',useCellNames=FALSE){
  ##metric is either 'tpm' or 'est_counts'
  if(!metric%in%c('tpm','est_counts')){
    print(paste(metric,'is not a valid kallisto output'))
    return(NULL)
  }
  if(buildFromFiles){
        allfiles<-rnaGencodeKallistoFiles()
        all.quants=sapply(allfiles,function(x){
          quants<-x[,metric]
          names(quants)=apply(x[,c(6,1,8)],1,paste,collapse='.')
          return(quants)
        })
        fname=paste('kallistoDerived',metric,'RNASeq_values.tsv',sep='')
        write.table(all.quants,file=fname,sep='\t')
        sf=File(fname,parentId=gencode.dir)
        synStore(sf,used='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/RNASeqData.R')
        tab=all.quants

    }else{
      if(metric=='est_counts')
        tab=read.table(synGet('syn5580378')@filePath)
      else if(metric=='tpm')
        tab=read.table(synGet('syn5580347')@filePath)

    }
  if(useCellNames){
    annotes=rnaAnnotations()
    colnames(tab)<-annotes[colnames(tab)]
  }
  return(tab)
}


getSampleNames<-function(rseqfiles){
    samp.names=samp.mappings$`Sample Name`[match(rseqfiles,samp.mappings$`RNA-Seq Data`)]
    if(all(is.na(samp.names)))
        samp.names=samp.mappings$`Sample Name`[match(rseqfiles,samp.mappings$`RNA-Seq Data (Gencode)`)]
    return(samp.names)
}

getSampleGenotype<-function(rseqfiles){
    samp.gen=samp.mappings$`Sample Genotype`[match(rseqfiles,samp.mappings$`RNA-Seq Data`)]
    if(all(is.na(samp.gen)))
        samp.gen=samp.mappings$`Sample Genotype`[match(rseqfiles,samp.mappings$`RNA-Seq Data (Gencode)`)]
    return(samp.gen)
}
                                        #this code was originally run on 2015/12/23, but now incorporated into this file
plotTranscriptsOfGene<-function(gene='NF1',count.mat,metric='tpm',dolog=FALSE,ttype=c()){
    require(pheatmap)

    samp.names=getSampleNames(colnames(count.mat))
    samp.gen=getSampleGenotype(colnames(count.mat))

    names(samp.gen)<-samp.names
  orig.samp.gen<-samp.gen
  colnames(count.mat)<-samp.names
  nf1=grep(paste("^",gene,".ENST0",sep=''),rownames(count.mat))
  print(paste("Found",length(nf1),gene,'transcripts'))
    n.mat<-count.mat[nf1,]
    zvals<-which(apply(n.mat,1,function(x) all(x==0)))
    if(length(ttype)>0){
        righttype=sapply(rownames(n.mat),function(x){
            gn=unlist(strsplit(x,split='.',fixed=T))
            gn[length(gn)]%in%ttype})
        n.mat=n.mat[which(righttype),]
    }
    if(length(zvals)>0 && length(zvals)!=nrow(n.mat))
        n.mat=n.mat[-zvals,]
    if(dolog)
        if(metric=='tpm')
            n.mat=log2(n.mat+0.01)
        else
            n.mat=log2(n.mat+1)
    fname=paste(gene,'TranscriptsBy',ifelse(dolog,'Log2',''),metric,paste(ttype,collapse='_'),'.png',sep='')
    pheatmap(n.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=orig.samp.gen),filename=fname)

}

plotPCA<-function(count.mat,metric='tpm',ttype=c()){
    require(ggbiplot)
    samp.names=getSampleNames(colnames(count.mat))
    samp.gen=getSampleGenotype(colnames(count.mat))

  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names

  zv<-which(apply(count.mat,1,var)==0)
  if(length(zv)>0)
      count.mat=count.mat[-zv,]

    if(length(ttype)>0.0){
        righttype=sapply(rownames(count.mat),function(x){
            gn=unlist(strsplit(x,split='.',fixed=T))
            gn[length(gn)]%in%ttype})
        count.mat=count.mat[which(righttype),]
    }
  pn<-prcomp(t(count.mat),center=T,scale=T)
  png(paste(metric,'valuesin',paste(ttype,collapse='_'),'RNASeqData.png',sep=''))
  p<-ggbiplot(pn,groups=samp.gen,var.axes=F)
  print(p)
  dev.off()

}


calcDiffEx<-function(){

}
