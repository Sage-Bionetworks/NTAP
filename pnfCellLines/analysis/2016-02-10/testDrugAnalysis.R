

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")

plotMostVariableAUCs()


##first get original auc vals
orig.ncats<-getValueForAllCells("FAUC")
orig.ctrp<-getCtrpScreensAsMatrix()

rescored.ncats<-getRecalculatedAUCMatrix()


rescored.ctrp<-ctrpDoseResponseCurve(FALSE)

##do drug-target clustering
ncats.targs<-ncatsDrugTargets()
ctrp.targs<-ctrpDrugTargets()

#get clusters, targets, and enrichment
orig.n.clust=getDrugClusters(orig.ncats,4,T)
orig.c.clust=getDrugClusters(orig.ctrp,30,T)

##now do p-value enrichment
n.enrich=clusterEnrichment(orig.n.clust,ncats.targs)
c.enrich=clusterEnrichment(orig.c.clust,ctrp.targs)

##do cell-genotype clustering
