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


##lastly, do the analysis on the combined matrix
auc.mat<-read.table(synGet('syn5714965')@filePath)
pnf.gt=dfiles$entity.sampleGenotype[match(colnames(auc.mat),dfiles$entity.sampleName)][1:8]
pnf.gt[is.na(pnf.gt)]<-'--'
pnf.gt=sapply(pnf.gt,function(x) if(x=='--'||x=='+-') return('-') else return('+'))
names(pnf.gt)<-colnames(auc.mat)[1:8]
full.gt=gt[colnames(auc.mat)]

over=intersect(c(names(pnf.gt),names(gt)),colnames(auc.mat))
comb.gt=c(pnf.gt,gt)[over]
am=auc.mat[,over]

tab<-aucDifferentialResponse(am,comb.gt)

sigs=rownames(tab)[which(tab$P.Value<0.005)]
write.table(tab,'aucValuesChangingAcrossNF1CombinedMatrixStatus.tab')
nam=as.matrix(am)
nam[which(is.na(nam))]<-mean(nam,na.rm=T)
pheatmap(nam[sigs,],#annotation_row=data.frame(Target=all.targs),
         clustering_distance_rows='correlation',
         annotation_col=data.frame(Genotype=comb.gt),cellwidth=10,cellheight=10,file='aucValuesDistinctAcrossNF1CombinedMatrixStatus_p005.png')

