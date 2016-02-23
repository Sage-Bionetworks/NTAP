source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")

aucMat=getCtrpScreensAsMatrix()

fzMat=apply(aucMat,1,function(x) apply(aucMat,1,function(y) fzCor(x,y)))
rownames(fzMat)<-colnames(fzMat)<-rownames(aucMat)
cormat<-stats::cor(t(aucMat),use='pairwise.complete.obs')
#dsNorm<-dSigTransform(fzMat)
pars=c(1,4,5,10,100,1000)
dsNorms<-lapply(pars,function(x) dSigTransform(fzMat,alpha=x))
names(dsNorms)<-as.character(pars)

##now how to we map that back to between -1 and 1?
do.test=FALSE
if(do.test){
  lt=which(fzMat<3,arr.ind=T)
  ct=which(cormat>0.9,arr.ind=T)
  ctm=apply(ct,1,function(x) intersect(which(lt[,1]==x[[1]]),which(lt[,2]==x[[2]]))>0)

##example of high correlation with low-ish fz score:
inds=ct[which(unlist(ctm)),]
no=which(cormat[inds]<0.9999)
inds=ct[no,]
diffs=t(apply(inds,1,function(x){
  print(paste(rownames(cormat)[x[1]],colnames(cormat)[x[2]]))
  cval=cormat[x[1],x[2]]
  fzval=fzMat[x[1],x[2]]
  print(paste("Cor:",cval))
  print(paste("FZ:",fzval))
  #now let's look at the number of samples
  na1=which(!is.na(aucMat[x[1],]))
  na2=which(!is.na(aucMat[x[2],]))
  print(paste("Found",length(intersect(na1,na2)),'overlapping samples!'))
  return(c(Cor=cval,FZ=fzval,Overlap=length(intersect(na1,na2))))}))
}

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
  na1=which(!is.na(aucMat[x[1],]))
  na2=which(!is.na(aucMat[x[2],]))
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
p<-ggplot(df,aes(x=Cor,y=FZ+10))+geom_point(aes(colour=Overlap))+scale_y_log10()
png('spearmanVsFisherZ.png')
print(p)
dev.off()

##ok, so this thing really is weeding out spurious correlations
##let's plot the double sigmoid transfer as well
p<-ggplot(df,aes(x=Cor,y=doubleSig_1))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA1.png')
print(p)
dev.off()

p<-ggplot(df,aes(x=Cor,y=doubleSig_4))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA4.png')
print(p)
dev.off()

p<-ggplot(df,aes(x=Cor,y=doubleSig_5))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA5.png')
print(p)
dev.off()

p<-ggplot(df,aes(x=Cor,y=doubleSig_10))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA10.png')
print(p)
dev.off()

p<-ggplot(df,aes(x=Cor,y=doubleSig_100))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA100.png')
print(p)
dev.off()

p<-ggplot(df,aes(x=Cor,y=doubleSig_1000))+geom_point(aes(colour=Overlap))+scale_colour_gradientn(colours=rainbow(5))
png('spearmanVsDoubSigTransA1000.png')
print(p)
dev.off()

##so, alpha 10 looks ok
##now get clusters at alpha of 10

##now lets store these results
afiles=list.files('.')
csvs=afiles[grep('png',afiles)]
cluster.dir='syn5674273'
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-12/testNormalization.R'
ncats.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ncatsSingleAgentScreens.R'
ctrp.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ctrpSingleAgentScreens.R'
analysis.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/singleDrugAnalysis.R'
for(csv in csvs){
  #first check to see if we're using CTRP or ncats, and original vs. rescored
  if(length(grep('ncats',csv,ignore.case=TRUE))>0){
    
    if(length(grep('rescored',csv,ignore.case=TRUE))>0){
      uf='syn5637634'
    }else{
      uf='syn5522627'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ncats.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='cluster parameter tests')
    
  }else if(length(grep('ctrp',csv,ignore.case=TRUE))>0){
    if(length(grep('rescored',csv,ignore.case=TRUE))>0){
      uf='syn5622708'
    }else{
      uf='syn5632189'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ctrp.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='cluster parameter test')
    
  }else{
    print(csv)
    
  }
}

