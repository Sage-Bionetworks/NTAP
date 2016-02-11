###
## Basic analysis of single agent drug screens
##
###
script.dir <- dirname(sys.frame(1)$ofile)

#' Get clusters of drugs from AUC values
#' @param aucMat matrix of AUC values - rows are
#' @param h Height at which to cut cluster
#' @return data frame of Drugs and Cluster to which they are assigned
getDrugClusters<-function(aucMat,h=4,doZScore=TRUE){
                                        #now zscore them
    if(doZScore){
        zsMat<-apply(aucMat,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
        nz.zs.mat=zsMat
                                        #reset NA values
    }else{
        nz.zs.mat=aucMat
    }

    nz.zs.mat[which(is.na(nz.zs.mat),arr.ind=T)]<-0.0

    ##there is a pretty blue cluster there, can we do any enrichment?
    drug.dists<-dist(nz.zs.mat)
    hc=hclust(drug.dists)

    ##now cut the clustering'
    drug.clusters<-cutree(hc,h=h)
                                        # hist(sapply(drug.clusters,length))
    return(data.frame(Drug=names(drug.clusters),Cluster=drug.clusters))
}

#' Get clusters of cells from AUC values (transpose of getDrugClusters)
#' @param aucMat matrix of AUC values - rows are
#' @param h Height at which to cut cluster
#' @return data frame of Drugs and Cluster to which they are assigned
getCellClusters<-function(aucMat,h=4,doZScore=TRUE){
    aucMat=t(aucMat)
                                        #now zscore them
    if(doZScore){
        zsMat<-apply(aucMat,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
        nz.zs.mat=zsMat
                                        #reset NA values
    }else{
        nz.zs.mat=aucMat
    }

    nz.zs.mat[which(is.na(nz.zs.mat),arr.ind=T)]<-0.0

    ##there is a pretty blue cluster there, can we do any enrichment?
    drug.dists<-dist(nz.zs.mat)
    hc=hclust(drug.dists)

    ##now cut the clustering'
    drug.clusters<-cutree(hc,h=h)
                                        # hist(sapply(drug.clusters,length))
    return(data.frame(Cell=names(drug.clusters),Cluster=drug.clusters))
}

#' Cluster enrichment analysis - according to ACME paper
#' @param clusters Data frame from getDrugClusters or getCellClusters function
#' @param clusterFeatures Data frame of drug/cell features to ask about
#' @param byDrug True if clustering by drug, otherwise will look for 'cell' feature
#' @param feature Column name of feature to query in clusterFeatures data frame
#' @return ??
clusterEnrichment<-function(clusters,clusterFeatures,byDrug=TRUE,feature='Target'){
    require(dplyr)

    ##join dataset
    ndf<-clusters
    if(byDrug)
        idx=match(ndf$Drug,clusterFeatures$Drug)
    else
        idx=match(ndf$Cell,clusterFeatures$Cell)

    ndf$Target=as.character(clusterFeatures[,feature])[idx]

    ndf=ndf[which(!is.na(ndf$Target)),]
    ndf=ndf[which(ndf$Target!=""),]
    print(paste('Reducing cluster set from',nrow(clusters),'to',nrow(ndf),'after removing NA and blank values'))
    drug.targs= ndf %>% group_by(Target,Cluster) %>% summarise(TimesTargetInCluster=n())
    tot.targs<-ndf %>% group_by(Target) %>% summarize(Total=n())
    clust.size<-ndf %>% group_by(Cluster) %>% summarize(Total=n())

    drug.targs$NumDrugsWithTarget=tot.targs$Total[match(drug.targs$Target,tot.targs$Target)]
    drug.targs$ClusterSize=clust.size$Total[match(drug.targs$Cluster,clust.size$Cluster)]
    targ.size=drug.targs %>% group_by(Cluster)%>% summarize(NumTargets=n())
    drug.targs$UniqueTargets=targ.size$NumTargets[match(drug.targs$Cluster,targ.size$Cluster)]

    drug.targs$Drug=ndf$Drug[match(drug.targs$Cluster,ndf$Cluster)]
    pvals=apply(drug.targs,1,function(x){
        over=as.numeric(x[['TimesTargetInCluster']])
        tt=as.numeric(x[['NumDrugsWithTarget']])
        cs=as.numeric(x[['ClusterSize']])
        mat=matrix(c(over,tt-over,cs-over,nrow(ndf)-tt-cs+over),nrow=2)
        return(fisher.test(mat,alt='g')$p.value)
    })

    drug.targs$Pvalue=pvals
    drug.targs$FDR=p.adjust(pvals,method='fdr')
    return(drug.targs)
}
