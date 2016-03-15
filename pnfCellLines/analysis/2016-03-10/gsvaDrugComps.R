
##compare drug response to gene expression pathway mappings

source("../../bin/ncatsSingleAgentScreens.R")
require(reshape2)
require(ggplot2)
##lastly let's look at GSVA clustrers
gsva=as.matrix(read.table(synGet('syn5689231')@filePath))
colnames(gsva)[1]='ipNF05.5 (mixed clone)'
colnames(gsva)[11]='ipNF05.5 (single clone)'

gsva.pvals=as.matrix(read.table(synGet('syn5689235')@filePath))
sig.enrich=which(apply(gsva.pvals,1,mean)<0.1)
gsva=gsva[sig.enrich,]

maxr.vals=getValueForAllCells("MAXR")
maxr.range<-apply(maxr.vals,1,function(x) max(x,na.rm=T)-min(x,na.rm=T))

auc.vals<-getValueForAllCells('FAUC')

this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-03-10/gsvaDrugComps.R'

#drug.by.cell=acast(getRecalculatedAUCTab(),Cell~Drug,value.var='AUC',fun.aggregate=function(x) mean(x,na.rm=T))
drug.by.cell <- t(auc.vals[names(sort(maxr.range,decreasing=T))[1:200],])

gsva.cells=intersect(colnames(gsva),rownames(drug.by.cell))
gsva.drugcor=cor(t(gsva[,gsva.cells]),drug.by.cell[gsva.cells,])

#highcor=which(abs(gsva.drugcor)>0.9,arr.ind=T)
allTops<-sapply(colnames(gsva.drugcor),function(x){
  ##first lets select by threshold
  tops<-names(which(abs(gsva.drugcor[,x])>0.95))
#  if(length(tops)>1)
#    tops=tops[which(apply(gsva[tops,],1,function(x) max(x,na.rm=T)-min(x,na.rm=T))>1.5)]
  if(length(tops)>12)
    tops<-names(sort(abs(gsva.drugcor[tops,x]),decreasing=T))[1:12]
  tops<-intersect(tops,rownames(gsva))
  if(length(tops)==0)
    return(NULL)
  restab<-do.call('rbind',lapply(tops,function(z) 
      data.frame(Pathway=rep(z,length(gsva.cells)),Cell=gsva.cells,Enrichment=gsva[z,gsva.cells],AUC=drug.by.cell[gsva.cells,x])))
  p<-ggplot(data.frame(restab))+geom_line(aes(x=AUC,y=Enrichment,col=Pathway))+ggtitle(paste('Pathways most correlated with',x))
  pdf(paste(gsub('/','',x),'gsvaCorrelation.pdf',sep='_'))
  print(p)
  hist(gsva.drugcor[,x],xlab='Pearson R',main=paste("GSVA Pathway Correlations with",x))
  dev.off()
})


##now do the same by pathway
allPaths<-sapply(rownames(gsva.drugcor),function(z){
  ##first lets select by threshold
  tops<-names(which(abs(gsva.drugcor[z,])>0.95))
  if(length(tops)>12)
    tops<-names(sort(abs(gsva.drugcor[z,tops]),decreasing=T))[1:12]
#  if(length(tops)>1)
  #  tops=tops[which(apply(gsva[tops],1,function(x) max(x,na.rm=T)-min(x,na.rm=T))>1.5)]
  
  tops<-intersect(tops,colnames(drug.by.cell))
  if(length(tops)<1)
    return(NULL)
  restab<-do.call('rbind',lapply(tops,function(x) 
    data.frame(Drug=rep(x,length(gsva.cells)),Cell=gsva.cells,Enrichment=gsva[z,gsva.cells],AUC=drug.by.cell[gsva.cells,x])))
  p<-ggplot(data.frame(restab))+geom_line(aes(x=AUC,y=Enrichment,col=Drug))+ggtitle(paste('Drugs most correlated with',z))
  pdf(paste(gsub('/','',z),'DrugCorrelation.pdf',sep='_'))
  print(p)
  hist(gsva.drugcor[z,],xlab='Pearson R',main=paste("Drug AUC Correlations with",z))
  dev.off()
})

path.cors=list.files(".")[grep('DrugCorrelation',list.files('.'))]
drug.cors=list.files(".")[grep('gsvaCorrelation',list.files('.'))]

used.list=Activity(executed=list(list(url=this.script)),used=list(list(entity='syn5689231'),list(entity='syn5689235')))

for(a in path.cors){
  synStore(File(path=a,parentId='syn5758476'),activity=used.list)
}

for(a in drug.cors){
  synStore(File(path=a,parentId='syn5758475'),activity=used.list)
}
