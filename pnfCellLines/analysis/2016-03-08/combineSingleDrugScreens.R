##combine ncats and CTRP recalculated matrices
#source("../../bin/ncatsSingleAgentScreens.R")
#source("../../bin/ctrpSingleAgentScreens.R")

pid='syn5714960'

ncats.recalc=read.table(synGet('syn5637634')@filePath,header=T)
ctp.recalc=read.table(synGet('syn5622708')@filePath,header=T)

##now do the matching again
nt.drugs<-unique(ncats.recalc$Drug)
drug_data=synGet('syn5632193')@filePath
drug_dat=read.table(drug_data,sep='\t',header=T,fill=T,quote='"')

##cell line infor
ccl_data=synGet('syn5632194')@filePath # '../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt'
ccl_dat=read.table(ccl_data,header=T,sep='\t')

ccl_metadata<-synGet('syn5632192')@filePath #../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt'
ccl_metadat<-read.table(ccl_metadata,sep='\t',header=T,as.is=T)

ccl=data.frame(Exp=ccl_dat$experiment_id,
               CCL_Id=ccl_dat$master_ccl_id,
               CCL_Name=ccl_metadat$ccl_name[match(ccl_dat$master_ccl_id,ccl_metadat$master_ccl_id)])

require(reshape2)
ctp.recalc$DrugName=drug_dat$cpd_name[match(ctp.recalc$master_cpd_id,drug_dat$master_cpd_id)]
ctp.recalc$CellLine=ccl$CCL_Name[match(ctp.recalc$experiment_id,ccl$Exp)]

drugmat=acast(ctp.recalc,DrugName~CellLine,value.var='unlist.res.auc.',fun.aggregate = mean)

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

ncmat<-acast(ncats.recalc,Drug~Cell,value.var="AUC",fun.aggregate = function(x) mean(x,na.rm=T))#t(sapply(nt.drugs,function(x) colSums(ncats.recalc[which(ncats.recalc$Drug==x),-1])))
rownames(ncmat)<-nt.drugs

comb.mat=cbind(ncmat[drug.matches,],
               drugmat[match(names(drug.matches),rownames(drugmat)),])

write.table(comb.mat,"ncatsCtrpCombinedMatrix_nplr_recalculated.txt")
synStore(File("ncatsCtrpCombinedMatrix_nplr_recalculated.txt",parentId=pid),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-03-08/combineSingleDrugScreens.R',wasExecuted=TRUE)))

