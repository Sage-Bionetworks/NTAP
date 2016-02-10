##compare to CTP dataa
source('../../bin/drugSensData.R')


require(plyr)
require(nplr)


all.aucs=sapply(dfiles$entity.sampleName,function(x) doseResponseCurve(x,NA,FALSE))
df<-as.data.frame(all.aucs)
df<-apply(all.aucs,2,unlist)
rownames(df)<-rownames(all.aucs)

write.table(df,file='aucRecalculatedFromNCATSscreens_nplr.txt',sep='\t')

synStore(File('aucRecalculatedFromNCATSscreens_nplr.txt',parentId='syn5522627'),executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-05/recalcAucFromNTAP.R')
