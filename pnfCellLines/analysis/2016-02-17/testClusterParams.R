##now we want to test the clustering parameters

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")

##first get original auc vals
orig.ncats<-getValueForAllCells("FAUC")
orig.ctrp<-getCtrpScreensAsMatrix()

##now get the rescored values
rescored.ncats<-acast(getRecalculatedAUCTab(),Drug~Cell,value.var="AUC",fun.aggregate = mean)
rescored.ctrp<-acast(ctrpDoseResponseCurve(FALSE),Drug~Cell,value.var="AUC",fun.aggregate = mean)

##do drug-target clustering
ncats.targs<-ncatsDrugTargets()
ctrp.targs<-ctrpDrugTargets()

#get clusters, targets, and enrichment
for (cheight in c(0.1,0.5,1,5)){
  orig.n.clust=getClusters(orig.ncats,cheight,byCol=FALSE,doubleSigAlpha=1)
  orig.c.clust=getClusters(orig.ctrp,cheight,byCol=FALSE,doubleSigAlpha=10)
  rescored.n.clust=getClusters(orig.ncats,cheight,byCol=FALSE,doubleSigAlpha=1)
  rescored.c.clust=getClusters(orig.ctrp,cheight,byCol=FALSE,doubleSigAlpha=10)
  
  ##now do p-value enrichment
  n.enrich=clusterEnrichment(orig.n.clust,ncats.targs)
  c.enrich=clusterEnrichment(orig.c.clust,ctrp.targs)
  rn.enrich=clusterEnrichment(rescored.n.clust,ncats.targs)
  rc.enrich=clusterEnrichment(rescored.c.clust,ctrp.targs)
  
  ##now write files
  write.table(n.enrich,file=paste('ncatsOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  write.table(c.enrich,file=paste('ctrpOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
 
  write.table(rn.enrich,file=paste('rescored_ncatsOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  write.table(rc.enrich,file=paste('rescored_ctrpOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  
  
}
##are we seeing meaningful target enrichment?
#View(subset(c.enrich[order(c.enrich$Cluster),],FDR<0.05))

#yes!!  but still need to do fisher z transform and double sigmoid.  