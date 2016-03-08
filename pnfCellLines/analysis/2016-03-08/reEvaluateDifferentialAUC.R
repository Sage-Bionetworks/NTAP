##repeat differential analysis of CTRP data with tissue of origin

source("../../bin/ctrpSingleAgentScreens.R")
source('../../bin/singleDrugAnalysis.R')

ccle.calls=getMSSMprocessedMutationCalls()
tissue.calls=getMSSMTissueOrigin()

auc.mat<-getCtrpScreensAsMatrix()

overlap=intersect(colnames(auc.mat),colnames(ccle.calls))
c.calls=ccle.calls[,overlap]
gt=rep('+',length(overlap))
gt[which(c.calls['NF1',]==1)]<-'-'
names(gt)<-overlap
origin=unlist(tissue.calls)[overlap]

auc.mat[is.nan(auc.mat)]<-NA
auc.mat[is.na(auc.mat)]<-mean(auc.mat,na.rm=T)

tab<-aucDifferentialResponse(auc.mat[,overlap],gt,otherFactors=data.frame(Origin=origin))

sigs=rownames(tab)[which(tab$P.Value<0.005)]

#drug.targs<-ctrpDrugTargets()
#all.targs<-as.character(drug.targs$Target)
#names(all.targs)<-drug.targs$Drug

#write to file
write.table(tab,'aucValuesChangingAcrossNF1Status.tab')
pheatmap(auc.mat[sigs,overlap],#annotation_row=data.frame(Target=all.targs),
         clustering_distance_rows='correlation',
         annotation_col=data.frame(Genotype=gt,Tissue=origin),cellwidth=10,cellheight=10,file='aucValuesDistinctAcrossNF1Status_p005.png')

##now separate out by individual tissue type
sapply(unique(origin),function(x){
    cells=names(which(origin==x))
    if(length(cells)<3)
      return(NULL)
    am=auc.mat[,cells]
    gtvals=gt[cells]
    if(length(unique(gtvals))==1)
      return(NULL)
    tab<-aucDifferentialResponse(am,gtvals)    
    sigs=rownames(tab)[which(tab$P.Value<0.005)]
    write.table(tab,paste('aucValuesChangingAcrossNF1StatusIn',x,'cellLines.tab',sep=''))
    if(length(sigs)>2)
      pheatmap(am[sigs,],#annotation_row=data.frame(Target=all.targs),
             clustering_distance_rows='correlation',
             annotation_col=data.frame(Genotype=gt[cells]),cellwidth=10,cellheight=10,file=paste('aucValuesDistinctAcrossNF1StatusIn',x,'cellLines_p005.png',sep=''))
    
})