#drug-target correlation

source("../../bin/crossDataComps.R")


##get all drug-target pairs

drug.to.targ<-getValueForAllCells('target')[,1]

alltargs<-setdiff(unique(unlist(sapply(drug.to.targ,strsplit,', '))),NA)


computeDrugTargCorrelation<-function(valname='MAXR',doLog=T,qthresh=0.1,pthresh=0.001){
  ##now compute correlation for each drug
  all.res<-lapply(alltargs,function(x) drugRna(valname=valname,gene=x,useGencode=T,
                                               doLog=doLog,collapseAllCounts=T,
                                               proteinCoding=T,qthresh=qthresh,
                                               pthresh=pthresh,doPlot=FALSE))
  
  ##now figure out which of those drugs shows up as targets
  drugs<-sapply(all.res,function(x) unique(unlist(as.character(x$DrugName))))
  names(drugs)<-alltargs
  
  over.res=sapply(alltargs,function(x) {
    dts=names(drug.to.targ)[grep(x,drug.to.targ)]
    dcors=drugs[[x]]
    over=intersect(dcors,dts)
    list(DrugsThatTarget=length(dts),DrugsThatCorrelate=length(dcors),overlap=length(over))
  })
  ##now write out table
  fname=paste('target',ifelse(doLog,'log2',''),'expressionVs',valname,'.csv',sep=',')
  write.table(over.res,fname,sep=',',row.names=F)
  return(over.res)
  
}

for(val in c('TAUC','MAXR'))
  for(dol in c(TRUE,FALSE))
    computeDrugTargCorrelation(valname=val,doLog=dol)