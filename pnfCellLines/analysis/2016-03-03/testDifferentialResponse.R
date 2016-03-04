##testlimmaon AUC data

source("../../bin/ncatsSingleAgentScreens.R")
source('../../bin/singleDrugAnalysis.R')



aucMat<-getValueForAllCells('FAUC')
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
sigs=rownames(tab)[which(tab$P.Value<0.005)]

#write to file
write.table(tab,'faucValuesNoHetDistinctAcrossGenotype.tab')
pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
         annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='faucValuesNoHetDistinctAcrossGenotype_p005.png')

##next do a one-allele check
##first, remove NA values, set to mean of matrix...
alt.auc=aucMat
alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)

##now re-assign gts
oas=list(`++`='+',`+-`='-',`--`='-')
alt.gt=oas[full.gt]
names(alt.gt)=names(full.gt)
#compute differential
tab<-aucDifferentialResponse(alt.auc,alt.gt)
sigs=rownames(tab)[which(tab$P.Value<0.005)]

#write to file
write.table(tab,'faucValuesWithHetDistinctAcrossGenotype.tab')
pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
         annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='faucValuesWithHetDistinctAcrossGenotype_p005.png')

####now do it all again for re-calculated matrix
#########
#########
aucMat<-getRecalculatedAUCMatrix()
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
sigs=rownames(tab)[which(tab$P.Value<0.005)]

#write to file
write.table(tab,'recalculatedValuesNoHetDistinctAcrossGenotype.tab')
pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
         annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='recalculatedValuesNoHetDistinctAcrossGenotype_p005.png')

##next do a one-allele check
##first, remove NA values, set to mean of matrix...
alt.auc=aucMat
alt.auc[is.na(aucMat)]<-mean(aucMat,na.rm=T)

##now re-assign gts
oas=list(`++`='+',`+-`='-',`--`='-')
alt.gt=oas[full.gt]
names(alt.gt)=names(full.gt)
#compute differential
tab<-aucDifferentialResponse(alt.auc,alt.gt)
sigs=rownames(tab)[which(tab$P.Value<0.005)]

#write to file
write.table(tab,'recalculatedValuesWithHetDistinctAcrossGenotype.tab')
pheatmap(aucMat[sigs,],annotation_row=data.frame(Target=all.targs[sigs]),
         annotation_col=data.frame(Genotype=full.gt),cellwidth=10,cellheight=10,file='recalculatedValuesWithHetDistinctAcrossGenotype_p005.png')




