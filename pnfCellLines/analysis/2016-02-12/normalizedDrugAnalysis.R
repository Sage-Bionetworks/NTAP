
##double checking these scripts to ensure they map back to the BROAD data

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")

#plotMostVariableAUCs()


##first get original auc vals
orig.ncats<-getValueForAllCells("FAUC")
orig.ctrp<-getCtrpScreensAsMatrix()

#rescored.ncats<-getRecalculatedAUCMatrix()

rescored.ncats<-acast(getRecalculatedAUCTab(),Drug~Cell,value.var="AUC",fun.aggregate = mean)
rescored.ctrp<-acast(ctrpDoseResponseCurve(FALSE),Drug~Cell,value.var="AUC",fun.aggregate = mean)

##do drug-target clustering
ncats.targs<-ncatsDrugTargets()
ctrp.targs<-ctrpDrugTargets()

#get clusters, targets, and enrichment
orig.n.clust=getClusters(orig.ncats,1,byCol=FALSE,doubleSigAlpha=10)
orig.c.clust=getClusters(orig.ctrp,0.8,byCol=FALSE,doubleSigAlpha=10)

##now do p-value enrichment
n.enrich=clusterEnrichment(orig.n.clust,ncats.targs)
c.enrich=clusterEnrichment(orig.c.clust,ctrp.targs)

##are we seeing meaningful target enrichment?
View(subset(c.enrich[order(c.enrich$Cluster),],FDR<0.05))

#yes!!  but still need to do fisher z transform and double sigmoid.  