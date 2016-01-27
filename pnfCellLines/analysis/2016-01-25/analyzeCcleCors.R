##get significanly correlated cell lines

require(synapseClient)
require(pheatmap)
require(preprocessCore)
require(parallel)
synapseLogin()

##first get all CCLE Data
synq=synapseQuery("select id,name from entity where parentId=='syn2325154'")

bcfiles<-synq[grep('bias_corrected.sf',synq[,1]),]

##now download all data
alldat<-mclapply(bcfiles[,2],function(x){
  print(paste('Getting/loading',x))
  tab<-read.table(synGet(x)@filePath)
  colnames(tab)<-c('Gene','Length','tpm','est_counts')
  return(tab)
})

fnames<-sapply(bcfiles[,1],function(x) gsub('_quant_bias_corrected.sf','',fixed=T,x))

names(alldat)<-fnames

tpms<-sapply(alldat,function(x){
  tpm=x[,3]
  names(tpm)<-x[,1]
  tpm
})

source("../../bin/RNASeqData.R")
##get our cell line data
tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=FALSE,useCellNames=TRUE)
#get enst values
tnames=sapply(rownames(tpm.mat),function(x) {
  arr=unlist(strsplit(x,split='.',fixed=T))
  arr[grep('ENST',arr)]})

#now do basic matching
t.idx=match(rownames(tpms),tnames)

##create larger matrix
all.tpms<-cbind(tpms[which(!is.na(t.idx)),],tpm.mat[t.idx[which(!is.na(t.idx))],])


#then get the correlation values we're interested in
res=synapseQuery("select * from entity where parentId=='syn5594111'")

for(qtile in c(0.95,0.99)){
    for(i in 1:nrow(res)){
        fname=res[[i,'entity.name']]
        synid=res[[i,'entity.id']]
        tab<-read.table(synGet(synid)@filePath)
                                        #get median 90th quantile
        corthresh=median(apply(tab,1,quantile,qtile))
        corcells=unique(unlist(apply(tab,1,function(x) names(which(x>corthresh)))))

        if(length(grep('qnormed',fname))>0)
            mod.tpm=normalize.quantiles(as.matrix(all.tpms))
        else
            mod.tpm=all.tpms
        colnames(mod.tpm)<-colnames(all.tpms)
        modcorcells<-sapply(corcells,function(x) paste('^',x,'$',sep=''))
                                        #now compute the correlation of the corcells only
        newcor=sapply(modcorcells,function(x){
            sapply(modcorcells,function(y) {cor(mod.tpm[,grep(x,colnames(mod.tpm))],
                                               mod.tpm[,grep(y,colnames(mod.tpm))],
                                                method=ifelse(length(grep('spearman',fname))>0,
                                                    'spearman','pearson'))})
        })

        pheatmap(newcor,
                 file=paste('top',1-qtile,'correlatedCellsFrom',gsub('.tsv','.pdf',fname),sep=''),
                 cellheight=10,cellwidth=10)
    }
}
