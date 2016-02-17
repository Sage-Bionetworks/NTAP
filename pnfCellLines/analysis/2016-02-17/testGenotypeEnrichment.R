##now let's test the genotype enrichment of the clusters of cells by AUC

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")

##first get original auc vals
orig.ncats<-getValueForAllCells("FAUC")
#orig.ctrp<-getCtrpScreensAsMatrix()

##now get the rescored values
rescored.ncats<-acast(getRecalculatedAUCTab(),Drug~Cell,value.var="AUC",fun.aggregate = mean)
#rescored.ctrp<-acast(ctrpDoseResponseCurve(FALSE),Drug~Cell,value.var="AUC",fun.aggregate = mean)

##do drug-genotype targeting
cell.genotypes<-dfiles[,c('entity.sampleName','entity.sampleGenotype')]
names(cell.genotypes)<-c('Cell','Genotype')

for (cheight in c(0.1,0.5,1,5)){
  orig.n.clust=getClusters(orig.ncats,h=NA,k=3,byCol=TRUE,doubleSigAlpha=20)
  rescored.n.clust=getClusters(orig.ncats,h=NA,k=3,byCol=TRUE,doubleSigAlpha=40)

  n.enrich=clusterEnrichment(orig.n.clust,cell.genotypes,byDrug=FALSE,feature='Genotype')
  rn.enrich=clusterEnrichment(rescored.n.clust,cell.genotypes,byDrug=FALSE,feature='Genotype')
}