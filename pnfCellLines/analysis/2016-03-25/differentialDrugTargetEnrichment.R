##testlimmaon AUC data

source("../../bin/ncatsSingleAgentScreens.R")
source('../../bin/singleDrugAnalysis.R')

pid='syn5838770'
all.sigs<-list()
het.sigs<-list()

targetEnrichment<-function(drugList){
  drug.targets<-ncatsDrugTargets()
  all.targs<-setdiff(as.character(drug.targets$Target)[match(drugList,drug.targets$Drug)],NA)
  pvals<-sapply(all.targs,function(a){
    all.targ.drugs<-which(as.character(drug.targets$Target)==a)
    drugs.selected=which(as.character(drug.targets$Drug)%in%drugList)
    drugs.with.target<-which(!is.na(drug.targets$Target))
    fm<-matrix(c(length(intersect(drugs.selected,all.targ.drugs)),length(setdiff(drugs.selected,all.targ.drugs)),
                 length(setdiff(all.targ.drugs,drugs.selected)),length(setdiff(drugs.with.target,union(drugs.selected,all.targ.drugs)))),nrow=2)
    list(NumDrugs=length(intersect(drugs.selected,all.targ.drugs)),Pval=fisher.test(fm,alt='g')$p.value)
          })
  res=data.frame(t(pvals))
  res=res[order(as.numeric(res$Pval),decreasing=F),]
  res$BHCorrectedP=p.adjust(as.numeric(res$Pval),method='BH')
  return(apply(res,2,unlist))
}

filelist=c()

##first get distributiuon of all scores
pdf('drugResponseDistribution.pdf')

for(val in c("TAUC","MAXR","LAC50")){
  
  aucMat<-getValueForAllCells(val)
  hist(aucMat,main=paste('Distribution of',val,'across cells'),xlab=val)
}
dev.off()

#now figure out how they compare
maxrMat<-getValueForAllCells("MAXR")
aucMat<-getValueForAllCells("TAUC")
lac50Mat<-getValueForAllCells("LAC50")

zscore<-function(x){
  (x-mean(x,na.rm=T))/sd(x,na.rm=T)
}

zscoreMat=zscore(maxrMat)+zscore(aucMat)+zscore(lac50Mat)
zscoreMat=zscoreMat/3


require(ggplot2)

auc.v.max=tidyr::gather(data.frame(cor(aucMat,maxrMat,use='pairwise.complete.obs')),"Cell","R",1:8)
auc.v.lac50=tidyr::gather(data.frame(cor(aucMat,lac50Mat,use='pairwise.complete.obs')),"Cell","R",1:8)
max.v.lac50=tidyr::gather(data.frame(cor(lac50Mat,maxrMat,use='pairwise.complete.obs')),"Cell","R",1:8)

auc.v.neglac50=tidyr::gather(data.frame(cor(aucMat,-1*lac50Mat,use='pairwise.complete.obs')),"Cell","R",1:8)
max.v.neglac50=tidyr::gather(data.frame(cor(-1*lac50Mat,maxrMat,use='pairwise.complete.obs')),"Cell","R",1:8)


full.df<-rbind(cbind(Comparison=rep('AUCvsMAXR',nrow(auc.v.max)),auc.v.max),
               cbind(Comparison=rep('AUCvsLAC50',nrow(auc.v.lac50)),auc.v.lac50),
               cbind(Comparison=rep('LAC50vsMAXR',nrow(max.v.lac50)),max.v.lac50),
               cbind(Comparison=rep('AUCvsnegLAC50',nrow(auc.v.neglac50)),auc.v.neglac50),
               cbind(Comparison=rep('negLAC50vsMAXR',nrow(max.v.neglac50)),max.v.neglac50))


p<-ggplot(full.df)+geom_boxplot(aes(y=R,fill=Comparison,x=Cell))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
png('correlationOfScores.png')
print(p)
dev.off()

for(val in c("ZSCORE","TAUC","MAXR","LAC50")){
  if(val=='ZSCORE'){
    aucMat<-zscoreMat
  }else{
  aucMat<-getValueForAllCells(val)

  }
  full.gt=dfiles$entity.sampleGenotype[match(colnames(aucMat),dfiles$entity.sampleName)]
  names(full.gt)<-colnames(aucMat)
  
  targs=ncatsDrugTargets()
  all.targs=targs$Target
  names(all.targs)<-targs$Drug
  
  ##first, remove NA values, set to mean of matrix...
  alt.auc=aucMat
  alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)
  
  ##remove het for this round
  alt.auc=alt.auc[,-6]
  alt.gt=full.gt[colnames(alt.auc)]
  
  #compute differential
  tab<-aucDifferentialResponse(alt.auc,alt.gt)
  write.table(tab,paste(val,'ValuesNoHetDistinctAcrossGenotype.tab',sep=''))
  
  
  ##next do a one-allele check
  ##first, remove NA values, set to mean of matrix...
  alt.auc=aucMat
  alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)
  
  ##now re-assign gts
  oas=list(`++`='+',`+-`='+',`--`='-')
  alt.gt=oas[full.gt]
  names(alt.gt)=names(full.gt)
  #compute differential
  alt.tab<-aucDifferentialResponse(alt.auc,alt.gt)
  write.table(alt.tab,paste(val,'ValuesWithHetDistinctAcrossGenotype.tab',sep=''))
  

  
  for(pval in c(0.01,0.005)){

    sigs=rownames(tab)[which(tab$P.Value<pval)]
    if(length(sigs)>2){
    sig.enrich=targetEnrichment(sigs)
    all.sigs[[val]]<-list(sigs)
    write.table(sig.enrich,file=paste(val,'EnrichedTargetsFrom_NoHetDistinctAcrossGenotype_p',pval,'.txt',sep=''),row.names=T,col.names=T)
    
    pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
           annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file=paste(val,'ValuesNoHetDistinctAcrossGenotype_p',pval,'.png',sep=''))
  }
  
    alt.sigs=rownames(tab)[which(alt.tab$P.Value<pval)]
   if(length(alt.sigs)>2){
   alt.sig.enrich=targetEnrichment(sigs)
    write.table(alt.sig.enrich,file=paste(val,'EnrichedTargetsFrom_WithHetDistinctAcrossGenotype_p',pval,'.txt',sep=''))
    het.sigs[[val]]<-list(alt.sigs)
    #write to file
    pheatmap(aucMat[alt.sigs,],annotation_row=data.frame(Target=all.targs[alt.sigs]),
             annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file=paste(val,'ValuesWithHetDistinctAcrossGenotype_p',pval,'.png',sep=''))
  }
  }
}
  
filelist=list.files('.')

filelist=filelist[-grep('differentialDrugTarget',filelist)]
for(f in filelist){
  synStore(File(f,parentId=pid),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-03-25/differentialDrugTargetEnrichment.R',wasExecuted=T)))
  
}
  ####now do it all again for re-calculated matrix
#########
# #########
# aucMat<-getRecalculatedAUCMatrix()
# full.gt=dfiles$entity.sampleGenotype[match(colnames(aucMat),dfiles$entity.sampleName)]
# names(full.gt)<-colnames(aucMat)
# 
# targs=ncatsDrugTargets()
# all.targs=targs$Target
# names(all.targs)<-targs$Drug
# 
# ##first, remove NA values, set to mean of matrix...
# alt.auc=aucMat
# alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)
# 
# ##remove het for this round
# alt.auc=alt.auc[,-6]
# alt.gt=full.gt[colnames(alt.auc)]
# 
# #compute differential
# tab<-aucDifferentialResponse(alt.auc,alt.gt)
# sigs=rownames(tab)[which(tab$P.Value<0.005)]
# 
# #write to file
# write.table(tab,'recalculatedValuesNoHetDistinctAcrossGenotype.tab')
# pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
#          annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='recalculatedValuesNoHetDistinctAcrossGenotype_p005.png')
# 
# ##next do a one-allele check
# ##first, remove NA values, set to mean of matrix...
# alt.auc=aucMat
# alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)
# 
# ##now re-assign gts
# oas=list(`++`='+',`+-`='-',`--`='-')
# alt.gt=oas[full.gt]
# names(alt.gt)=names(full.gt)
# #compute differential
# tab<-aucDifferentialResponse(alt.auc,alt.gt)
# sigs=rownames(tab)[which(tab$P.Value<0.005)]
# 
# #write to file
# write.table(tab,'recalculatedValuesWithHetDistinctAcrossGenotype.tab')
# pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
#          annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='recalculatedValuesWithHetDistinctAcrossGenotype_p005.png')
# 
# 


