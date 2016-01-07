###
###
##
## RNASeqData.R
## Basic set of library functions designed to collect and manage RNA-Seq data for NTAP PNF cell lines
##
###

rna.dir<-'syn5562003'

require(synapseClient)
require(data.table)

synapseLogin()

rnaKallistoFiles<-function(){
  rfiles<-synapseQuery(paste('select * from entity where parentId=="',rna.dir,'"',sep=''))
  rfiles<-rfiles[which(!is.na(rfiles$entity.sampleName)),]
  all.files<-lapply(rfiles$entity.id,function(x)
    as.data.frame(fread(synGet(x)@filePath))
  )
  names(all.files)<-rfiles$entity.id
  return(all.files)
}

rnaAnnotations<-function(){
  rfiles<-synapseQuery(paste('select id,sampleName from entity where parentId=="',rna.dir,'"',sep=''))
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
        sf=File(fname,parentId=rna.dir)
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
  

