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
stats<-list()
for (cheight in c(0.1,0.5,1)){
  orig.n.clust=getClusters(orig.ncats,cheight,byCol=FALSE,doubleSigAlpha=1)
  orig.c.clust=getClusters(orig.ctrp,cheight,byCol=FALSE,doubleSigAlpha=10)
  rescored.n.clust=getClusters(orig.ncats,cheight,byCol=FALSE,doubleSigAlpha=1)
  rescored.c.clust=getClusters(orig.ctrp,cheight,byCol=FALSE,doubleSigAlpha=10)
  
  ##now do p-value enrichment
  n.enrich=clusterEnrichment(orig.n.clust,ncats.targs)
  c.enrich=clusterEnrichment(orig.c.clust,ctrp.targs)
  rn.enrich=clusterEnrichment(rescored.n.clust,ncats.targs)
  rc.enrich=clusterEnrichment(rescored.c.clust,ctrp.targs)
  ##first NCATS
  stats$Dataset<-c(stats$Dataset,rep("NCATS",2))
  stats$Alpha<-c(stats$Alpha,rep(1,2))
  stats$ClusterHeight=c(stats$ClusterHeight,rep(cheight,2))
  #original
  stats$NumClusters=c(stats$NumClusters,length(unique(n.enrich$Cluster)))
  stats$CurveCalc=c(stats$CurveCalc,'Original')
  stats$NumSig=c(stats$NumSig,length(which(n.enrich$FDR<0.05)))
  #nplr
  stats$NumClusters=c(stats$NumClusters,length(unique(rn.enrich$Cluster)))
  stats$CurveCalc=c(stats$CurveCalc,'NPLR')
  stats$NumSig=c(stats$NumSig,length(which(rn.enrich$FDR<0.05))) 
  ##then CTRP
  stats$Dataset<-c(stats$Dataset,rep("CTRP",2))
  stats$Alpha<-c(stats$Alpha,rep(10,2))
  stats$ClusterHeight=c(stats$ClusterHeight,rep(cheight,2))
  #original
  stats$NumClusters=c(stats$NumClusters,length(unique(c.enrich$Cluster)))
  stats$CurveCalc=c(stats$CurveCalc,'Original')
  stats$NumSig=c(stats$NumSig,length(which(c.enrich$FDR<0.05)))
  #nplr
  stats$NumClusters=c(stats$NumClusters,length(unique(rc.enrich$Cluster)))
  stats$CurveCalc=c(stats$CurveCalc,'NPLR')
  stats$NumSig=c(stats$NumSig,length(which(rc.enrich$FDR<0.05))) 
  
  
  ##now write files
  write.table(n.enrich,file=paste('ncatsOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  write.table(c.enrich,file=paste('ctrpOriginal_ds10_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
 
  write.table(rn.enrich,file=paste('rescored_ncatsOriginal_ds1_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  write.table(rc.enrich,file=paste('rescored_ctrpOriginal_ds10_clusHeight',cheight,'clusterPvals.csv',sep=''),sep=',')
  
  
}
write.table(as.data.frame(stats),file='ClusterStatsByHeight.csv',sep=',',row.names=F)

afiles=list.files('.')
cluster.dir='syn5674273'
csvs=afiles[grep("csv",afiles)]
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-17/testClusterParams.R'
ncats.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ncatsSingleAgentScreens.R'
ctrp.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ctrpSingleAgentScreens.R'
analysis.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/singleDrugAnalysis.R'
for(csv in csvs){
  #first check to see if we're using CTRP or ncats, and original vs. rescored
  if(length(grep('ncats',csv))>0){
        
    if(length(grep('rescored',csv))>0){
      uf='syn5637634'
    }else{
      uf='syn5522627'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ncats.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='drug AUC Clustering')
    
  }else if(length(grep('ctrp',csv))>0){
    if(length(grep('rescored',csv))>0){
      uf='syn5622708'
    }else{
      uf='syn5632189'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ctrp.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='drug AUC Clustering')
    
  }else{
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ncats.script,wasExecuted=TRUE),
                          list(url=ctrp.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE)),
             activityName='drug AUC Clustering Summary')
    
  }
}