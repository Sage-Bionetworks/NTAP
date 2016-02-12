###
## Basic analysis of single agent drug screens
##
###
script.dir <- dirname(sys.frame(1)$ofile)

#' Get clusters of drugs from AUC values
#' @param aucMat matrix of AUC values - rows are
#' @param h Height at which to cut cluster
#' @return data frame of Drugs and Cluster to which they are assigned
getDrugClusters<-function(aucMat,h=4,doubleSigAlpha=NA,doZScore=TRUE){
                                        #now zscore them
  
  ##THIS FUNCTION IS DEPRACATED
  if(!is.na(doubleSigAlpha))
    return( getClusters(aucMat,h,byCol=FALSE,doubleSigAlpha))
  
  if(doZScore){
        zsMat<-apply(aucMat,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
        nz.zs.mat=zsMat
                                        #reset NA values
    }else{
        nz.zs.mat=aucMat
    }

    nz.zs.mat[which(is.na(nz.zs.mat),arr.ind=T)]<-0.0

    ##there is a pretty blue cluster there, can we do any enrichment?
    drug.dists<-as.dist(1-cor(t(nz.zs.mat)))
    hc=hclust(drug.dists)

    ##now cut the clustering'
    drug.clusters<-cutree(hc,h=h)
                                        # hist(sapply(drug.clusters,length))
    return(data.frame(Drug=names(drug.clusters),Cluster=drug.clusters))
}

#' Get clusters of AUC values by row or column
#' @param aucMat matrix of AUC values - rows are
#' @param h Height at which to cut cluster
#' @param byCol set to true to analyze by column (cells) otherwise will cluster rows (drugs)
#' @param doubleSigAlpha - alpha value to use for double sigmoid. if NA, no transform is done
#' @return data frame of Drugs and Cluster to which they are assigned
getClusters<-function(aucMat,h=4,byCol=FALSE,doubleSigAlpha=NA){
    
    #aucMat=t(aucMat)
                                        #now zscore them
    #  if(doZScore){
    #      zsMat<-apply(aucMat,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
    #       nz.zs.mat=zsMat
                                       #reset NA values
    #  }else{
    #      nz.zs.mat=aucMat
    # }
  
    #always cluster by column, so transpose if necessary
    if(!byCol)
      aucMat=t(aucMat)
    
    ##now establish correlations
    if(is.na(doubleSigAlpha)){
      corvals=cor(aucMat,use='pairwise.complete.obs')    
    }else{
      fz.mat=apply(aucMat,2,function(x)
          apply(aucMat,2,function(y)
            fzCor(x,y)))
      rownames(fz.mat)<-colnames(fz.mat)<-colnames(aucMat)
      corvals=dsTransform(fz.mat,alpha=doubleSigAlpha)
    }
    

    corvals[which(is.na(corvals),arr.ind=T)]<-0.0

    ##there is a pretty blue cluster there, can we do any enrichment?
    drug.dists<-as.dist(1-corvals)
    hc=hclust(drug.dists)

    ##now cut the clustering'
    drug.clusters<-cutree(hc,h=h)
                                        # hist(sapply(drug.clusters,length))
    return(data.frame(Cell=names(drug.clusters),Cluster=drug.clusters))
}

#' Adjusted correlation - does standard pearson with fisher z transform
#' and a double signmoid transform to enable clustering
#' @param vec1: first vector
#' @param vec2: second vector
#' 
fzCor<-function(vec1,vec2){
  ##first collect all measured values
  m.vals=intersect(which(!is.nan(vec1)),which(!is.nan(vec2)))
  if(length(m.vals)<4)
    return(NA)
  #compute r
  r=cor(vec1[m.vals],vec2[m.vals])
  if(r==1.0)
    r=0.999999
  #now do fisher z transform
  fz = (log((1+r)/(1-r))/2)*sqrt(length(m.vals)-3)
  fz
}

#'Double Sigmoid Transform
#'Takes two parameters
#'@param k = 3 in both CD and JBMS papers
#'@param alpha = 2.5ish in JBMS papers
dSigTransform<-function(mat,k=3,alpha=2.5){
  #plim=10^((-0.95)*(-log10(min(pnorm(mat),na.rm=TRUE))+log10(0.05))-log10(0.05))
  #zlim=quantile(mat,plim,na.rm=T)
  dval=((mat/alpha)^k)/sqrt(1+(mat/alpha)^(2*k))
  dval
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
    if(byDrug){
        idx=lapply(as.character(ndf$Drug),function(x) which(clusterFeatures$Drug==x))
        names(idx)<-ndf$Drug
        ndf=do.call('rbind',apply(ndf,1,function(x) {
          i=idx[[x[['Drug']]]]
          data.frame(Drug=rep(x[['Drug']],length(i)),
                     Cluster=rep(x[['Cluster']],length(i)),
                     Target=clusterFeatures$Target[i])}))
        
    }else{
        idx=lapply(ndf$Cell,function(x) which(clusterFeatures$Cell==x))
        names(idx)<-ndf$Cell
        ndf=do.call('rbind',apply(ndf,1,function(x) {
          i=idx[[x[['Cell']]]]
          data.frame(Drug=rep(x[['Cell']],length(i)),
                     Cluster=rep(x[['Cluster']],length(i)),
                     Target=clusterFeatures$Target[i])}))
        
    }
#    ndf$Target=as.character(clusterFeatures[,feature])[idx]
       
    ndf=ndf[which(!is.na(ndf$Target)),]
    ndf=ndf[which(ndf$Target!=""),]
    print(paste('Altering cluster set from',nrow(clusters),'to',nrow(ndf),'after removing NA and blank values and accounting for multiple drug targets'))
    drug.targs= ndf %>% group_by(Target,Cluster) %>% summarise(TimesTargetInCluster=n())
    tot.targs<-ndf %>% group_by(Target) %>% summarize(Total=n())
    clust.size<-ndf %>% group_by(Cluster) %>% summarize(Total=n())

    drug.targs$NumDrugsWithTarget=tot.targs$Total[match(drug.targs$Target,tot.targs$Target)]
    drug.targs$ClusterSize=clust.size$Total[match(drug.targs$Cluster,clust.size$Cluster)]
    targ.size=drug.targs %>% group_by(Cluster)%>% summarize(NumTargets=n())
    drug.targs$UniqueTargets=targ.size$NumTargets[match(drug.targs$Cluster,targ.size$Cluster)]

    drug.targs$DrugsWithTarget=apply(drug.targs,1,function(x){
        tc=which(ndf$Target==x[['Target']])
        cc=which(ndf$Cluster==x[['Cluster']])
        paste(ndf$Drug[intersect(tc,cc)],collapse=';')})
    
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
