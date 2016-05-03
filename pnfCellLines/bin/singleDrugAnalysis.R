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
getClusters<-function(aucMat,h=4,k=NA,byCol=FALSE,doubleSigAlpha=NA){

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
      corvals=dSigTransform(fz.mat,alpha=doubleSigAlpha)
    }


    corvals[which(is.na(corvals),arr.ind=T)]<-0.0

    ##there is a pretty blue cluster there, can we do any enrichment?
    drug.dists<-as.dist(1-corvals)
    hc=hclust(drug.dists)

    ##now cut the clustering'
    if(is.na(h)){
      drug.clusters<-cutree(hc,k=k)
      print(paste('Found',length(unique(drug.clusters)),'clusters at k',k))
    }
    else{
      drug.clusters<-cutree(hc,h=h)
      print(paste('Found',length(unique(drug.clusters)),'clusters at height',h))
      }                                   # hist(sapply(drug.clusters,length))
    if(byCol)
      return(data.frame(Cell=names(drug.clusters),Cluster=drug.clusters))
    else
      return(data.frame(Drug=names(drug.clusters),Cluster=drug.clusters))

}

#' Adjusted correlation - does standard pearson with fisher z transform
#' and a double signmoid transform to enable clustering
#' @param vec1: first vector
#' @param vec2: second vector
#'
fzCor<-function(vec1,vec2){
  ##first collect all measured values
  m.vals=intersect(which(!is.nan(vec1)),which(!is.nan(vec2)))
  m.vals=intersect(m.vals,intersect(which(!is.na(vec1)),which(!is.na(vec2))))
  if(length(m.vals)<4)
    return(NA)
 # print(length(m.vals))
  #compute r
  r=cor(vec1[m.vals],vec2[m.vals])
  if(is.na(r))
    return(NA)
  else if(r==1.0)
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
                     Target=clusterFeatures[[feature]][i])}))

    }else{
        idx=lapply(ndf$Cell,function(x) which(clusterFeatures$Cell==x))
        names(idx)<-ndf$Cell
        ndf=do.call('rbind',apply(ndf,1,function(x) {
          i=idx[[x[['Cell']]]]
          data.frame(Drug=rep(x[['Cell']],length(i)),
                     Cluster=rep(x[['Cluster']],length(i)),
                     Target=clusterFeatures[[feature]][i])}))

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

#' test normalization parameters to ensure distributions are good
#' @param matrix of interest to cluster by column or row
#' @byCol if true will evaluate columns not rows
#' @alphas defaults to c(1,5,10,100)
#' @return
testNormalizationParameters<-function(mat,byCol=FALSE,alphas=c(1,5,10,100),prefix=''){
  if(!byCol)
    mat<-t(mat)
  fzMat=apply(mat,2,function(x) apply(mat,2,function(y) fzCor(x,y)))
  rownames(fzMat)<-colnames(fzMat)<-colnames(mat)

  cormat<-stats::cor(mat,use='pairwise.complete.obs')
  pars=alphas
  dsNorms<-lapply(pars,function(x) dSigTransform(fzMat,alpha=x))
  names(dsNorms)<-as.character(pars)


  ##now compute normalization for all values, and compare
  all.inds=do.call("rbind",lapply(1:nrow(cormat),function(x) cbind(rep(x,nrow(cormat)),1:nrow(cormat))))

  all.diffs=t(apply(all.inds,1,function(x){
    #print(paste(rownames(cormat)[x[1]],colnames(cormat)[x[2]]))
    cval=cormat[x[1],x[2]]
    fzval=fzMat[x[1],x[2]]
    #ds=dsNorm[x[1],x[2]]
    #print(paste("Cor:",cval))
    #print(paste("FZ:",fzval))
    #now let's look at the number of samples
    na1=which(!is.na(mat[,x[1]]))
    na2=which(!is.na(mat[,x[2]]))
    #print(paste("Found",length(intersect(na1,na2)),'overlapping samples!'))
    if(!is.na(fzval) && x[1]!=x[2] && fzval>500)
      print(paste("Check:",rownames(cormat)[x[1]],colnames(cormat)[x[2]]))
    retvec=c(Cor=cval,FZ=fzval,Overlap=length(intersect(na1,na2)),Ind1=x[1],Ind2=x[2])
    ds=unlist(lapply(dsNorms,function(y) y[x[1],x[2]]))
    names(ds)<-paste('doubleSig',pars,sep='_')
    return(c(retvec,ds))}))

  df=data.frame(all.diffs)

  ##now plot the results
  library(ggplot2)
  p<-ggplot(df,aes(x=Cor,y=FZ+10))+geom_point(aes(colour=Overlap))+scale_y_log10()+scale_colour_gradientn(colours=rainbow(5))
  png(paste(prefix,'pearsonVsFisherZ.png',sep='_'))
  print(p)
  dev.off()

  ##ok, so this thing really is weeding out spurious correlations
  ##let's plot the double sigmoid transfer as well
  for(a in alphas){
    ddf=df[,c("Cor","Overlap",paste('doubleSig',a,sep='_'))]
    colnames(ddf)<-c("R","SharedMeasurements","doubleSigTransform")
    p<-ggplot(ddf,aes(x=R,y=doubleSigTransform))+geom_point(aes(colour=SharedMeasurements))+scale_colour_gradientn(colours=rainbow(5))
    png(paste(prefix,'_pearsonVsDoubSigTransAlpha',a,'.png',sep=''))
    print(p)
    dev.off()
  }
}


aucPlsr<-function(aucMat,cellClasses,otherFactors=NULL){
    require(pls)
  facts=factor(unlist(cellClasses))
  levs=levels(facts)
  
  
}

#' Check differential auc values between subsets of cells
#' @param aucMat  - matrix of AUC values to be used as input
#' @param cellClasses - list of cell classes - two at most, names should match colnames of aucMat
#'
aucDifferentialResponse<-function(aucMat,cellClasses,otherFactors=NULL){
    require(limma)
  #  all.cells=unlist(cellClasses)

    facts=factor(unlist(cellClasses))
    levs=levels(facts)
    design= model.matrix(~facts)
    
    if(!is.null(otherFactors)){
      facts=cbind(Class=facts,otherFactors)
      design= model.matrix(~.,facts)
      
    }
    colnames(design)[1:2]=levs
    fit <- lmFit(aucMat, design)
    fit <- eBayes(fit)
    tab <- topTable(fit, coef=levs[length(levs)],number=Inf)
    return(tab)
}
