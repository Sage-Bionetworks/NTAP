##compare drug data

require(synapseClient)
require(reshape2)
synapseLogin()
library(pheatmap)
ncats.recalc=read.table(synGet('syn5637634')@filePath,header=T)
ctp.recalc=read.table(synGet('syn5622708')@filePath,header=T)

##now do the matching again
nt.drugs<-unique(ncats.recalc$Drug)
# 
# drug_data='../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt'
# drug_dat=read.table(drug_data,sep='\t',header=T,fill=T,quote='"')
# 
# ##cell line infor
# ccl_data='../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt'
# ccl_dat=read.table(ccl_data,header=T,sep='\t')
# 
# ccl_metadata<-'../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt'
# ccl_metadat<-read.table(ccl_metadata,sep='\t',header=T,as.is=T)
# 
# ccl=data.frame(Exp=ccl_dat$experiment_id,
#                CCL_Id=ccl_dat$master_ccl_id,
#                CCL_Name=ccl_metadat$ccl_name[match(ccl_dat$master_ccl_id,ccl_metadat$master_ccl_id)])

require(reshape2)
ctp.recalc$DrugName=drug_dat$cpd_name[match(ctp.recalc$master_cpd_id,drug_dat$master_cpd_id)]
ctp.recalc$CellLine=ccl$CCL_Name[match(ctp.recalc$experiment_id,ccl$Exp)]

drugmat=acast(ctp.recalc,DrugName~CellLine,value.var='unlist.res.auc.',fun.aggregate = mean)
ncats.mat<-acast(ncats.recalc,Drug~Cell,value.var='AUC',fun.aggregat = mean)

matched.drug.names=sapply(rownames(drugmat),function(x){
  mval=NA
  mval=match(x,nt.drugs)
  if(!is.na(mval))
    return(mval)
  
  y=paste('^',x,'$',sep='')
  mval=grep(y,nt.drugs,ignore.case=T)
  if(length(mval)==1)
    return(mval)
  
  mval=grep(gsub('-','',x),nt.drugs,ignore.case=T)
  if(length(mval)==1)
    return(mval)
  
  return(NULL)
})

drug.matches=unlist(matched.drug.names)
print(paste("Found",length(drug.matches),'drugs that are found in both NCATS and CCL screens'))
#create a combined matrix

#ncmat<-t(sapply(nt.drugs,function(x) colSums(ncats.mat[which(ncats.recalc$Drug==x),-1])))
#rownames(ncmat)<-nt.drugs
ncmat<-ncats.mat[nt.drugs,]

comb.mat=cbind(ncmat[drug.matches,],
               drugmat[match(names(drug.matches),rownames(drugmat)),])


missing=which(apply(comb.mat,1,function(x) length(which(is.nan(x))))>100)
print(paste("Removing",length(missing),'drugs with too many missing values'))
comb.mat=comb.mat[-missing,]

comb.mat[which(is.na(comb.mat),arr.ind=T)]<-0.0
pheatmap(comb.mat,filename='recalcAucClusters.png',cellwidth=10,cellheight = 10)

##that's a big cluster, let's try some basic filtering
cormat=cor(comb.mat)
top50=order(apply(cormat[1:8,],2,mean),decreasing=T)[1:50]

pheatmap(comb.mat[,top50],filename='recalcAucClustersTop50.png',
         cellwidth=10,cellheight = 10,
         clustering_distance_rows='correlation',
         clustering_distnace_cols='correlation')

##now do some spearman
cormat=cor(comb.mat,method='spearman')
top50=order(apply(cormat[1:8,],2,mean),decreasing=T)[1:50]

pheatmap(comb.mat[,top50],filename='recalcAucClustersTop50spearman.png',
         cellwidth=10,cellheight = 10,
         clustering_distance_rows='correlation',
         clustering_distnace_cols='correlation')

cormat=cor(comb.mat)
top25=order(apply(cormat[1:8,],2,mean),decreasing=T)[1:25]

pheatmap(comb.mat[,top50],filename='recalcAucClustersTop25.png',
         cellwidth=10,cellheight = 10,
         clustering_distance_rows='correlation',
         clustering_distnace_cols='correlation')

##now do some spearman
cormat=cor(comb.mat,method='spearman')
top25=order(apply(cormat[1:8,],2,mean),decreasing=T)[1:25]

pheatmap(comb.mat[,top25],filename='recalcAucClustersTop25spearman.png',
         cellwidth=10,cellheight = 10,
         clustering_distance_rows='correlation',
         clustering_distnace_cols='correlation')

pngs=list.files('.')[grep('png',list.files('.'))]
for(p in pngs){
  sf=File(p,parentId='syn5763353')
  synStore(sf,used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-09/clusterRecalcAUCs.R',wasExecuted=TRUE)))
}
