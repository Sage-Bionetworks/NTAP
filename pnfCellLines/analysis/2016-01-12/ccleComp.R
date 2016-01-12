source('../../bin/RNASeqData.R')

tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE)
ecounts.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,metric='est_counts')

##now get the ccle data
synq=synapseQuery("select id,name from entity where parentId=='syn2325154'")

bcfiles<-synq[grep('bias_corrected.sf',synq[,1]),]

##now download all data
alldat<-sapply(bcfiles[,2],function(x){
  tab<-read.table(synGet(x)@filePath)
  colnames(tab)<-c('Gene','Length','tpm','est_counts')
  return(tab)
})

fnames<-sapply(bcfiles[,1],function(x) gsub('_quant_biased_corrected.sf','',x))



