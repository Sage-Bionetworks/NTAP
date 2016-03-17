##testing new drug comparisons
source("../../bin/crossDataComps.R")

nf1cor.t=drugGeneNorm(valname='FAUC',gene='NF1',useGencode=TRUE,by.pval=F)
nf1cor.t=drugGeneNorm(valname='FAUC',gene='MNX1',useGencode=TRUE,by.pval=F)
nf1cor.t=drugGeneNorm(valname='FAUC',gene='TP53',useGencode=TRUE,by.pval=F)

targs<-ncatsDrugTargets()
targs<-lapply(setdiff(unique(as.character(targs$Target)),NA),function(g){
  res=drugGeneNorm(valname='FAUC',gene=g,useGencode=TRUE,by.pval=F)
  num.over<-grep(g,res$DrugName)
  print(num.over)
  if(length(num.over)==0)
    return(0)
  else
    return(length(unique(res$DrugName[num.over])))
})
