##compare to CTP dataa
source('../../bin/drugSensData.R')


require(plyr)
require(nplr)


all.aucs=lapply(dfiles$entity.sampleName,function(x) doseResponseCurve(x,NA,FALSE))

names(all.aucs)<-dfiles$entity.sampleName

df=do.call("rbind",all.aucs)
#df<-as.data.frame(all.aucs)
#df<-apply(all.aucs,2,unlist)
#rownames(df)<-rownames(all.aucs)

write.table(df,file='aucRecalculatedFromNCATSscreens_nplr.txt',sep='\t',row.names=F)

synStore(File('aucRecalculatedFromNCATSscreens_nplr.txt',parentId='syn5522627'),executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-05/recalcAucFromNTAP.R')
