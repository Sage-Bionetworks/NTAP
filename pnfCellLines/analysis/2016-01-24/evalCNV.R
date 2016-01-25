##evalCnvData.


source("../../bin/clusterCNVBySample.R")
highthresh=idDiffVals(segdata2,thresh=-1)

midthresh=idDiffVals(segdata2,thresh=-0.5)
midthresh=idDiffVals(segdata2,thresh=-0.5,byChrom=T)

afiles=list.files('./')
sapply(afiles[grep('.png',afiles)], function(x){
  nf=File(x,parentId='syn5605333')
  synStore(nf,executed=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/clusterCNVBySample.R'),list(url="https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-01-24/evalCNV.R")))
})