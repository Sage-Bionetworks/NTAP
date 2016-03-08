##this file loads up the CTRP single agent screens downloaded from
## ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/Broad/CTRPv1.0_2013_pub_Cell_154_1151/
##uploded to synapse at https://www.synapse.org/#!Synapse:syn5622707
library(synapseClient)
synapseLogin()
##AUC data
auc_data<-synGet('syn5622711')@filePath
auc_dat=read.table(auc_data,sep='\t',header=T,quote='"')

##Cell line data per experiment and meta data
ccl_data=synGet('syn5632194')@filePath
ccl_dat=read.table(ccl_data,header=T,sep='\t')
ccl_metadata<-synGet('syn5632192')@filePath
ccl_metadat<-read.table(ccl_metadata,sep='\t',header=T,as.is=T)

##Drug meta data
drug_data=synGet("syn5632193")@filePath #'../../../CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt'
drug_dat=read.table(drug_data,sep='\t',header=T,fill=T,quote='"')


#now match it all up into a single data frame
new.df<-data.frame(AUC=as.numeric(auc_dat$area_under_curve),
                   Drug=drug_dat$cpd_name[match(auc_dat$master_cpd_id,drug_dat$master_cpd_id)],
                   CCL_id=ccl_dat$master_ccl_id[match(auc_dat$experiment_id,ccl_dat$experiment_id)])
new.df$CCL=ccl_metadat$ccl_name[match(new.df$CCL_id,ccl_metadat$master_ccl_id)]

#' Collect CTRP data as matrix
#' @return drug by cell line AUC value matrix
getCtrpScreensAsMatrix<-function(){
    require(reshape2)
    drugmat<-acast(new.df,Drug~CCL,value.var="AUC",fun.aggregate=function(x) mean(x,na.rm=T))#,fill=NA)
    drugmat
}

#' Plot most variable values across all cells
#' @param valname
#' @param ft Format of file
#' @return drug values
plotMostVariableAUCs<-function(ft='png'){
    drug.values<-getCtrpScreensAsMatrix()
    mv=drug.values[order(apply(drug.values,1,var,na.rm=T),decreasing=T)[1:50],]
    mv[which(is.nan(mv),arr.ind=T)]<-0.0
    #instead of genotype, let's get the cell origin
    gt<-ccl_metadat$ccle_primary_site[match(colnames(drug.values),ccl_metadat$ccl_name)]
    names(gt)<-colnames(drug.values)

    pheatmap(t(mv),annotation_row=data.frame(CellType=gt),
             cellwidth=10,cellheight=10,file=paste('drugsWithMostVariableAUCAcrossCellLines.',ft,sep=''))
    return(drug.values)

}

#' Recompute the dose response curve from the data points using NPLR
#' @param recalculate Set to true to recalculate, otherwise will pull from synapse
#' @return Data frame of each cell, drug and value
ctrpDoseResponseCurve<-function(recalculate=TRUE,as.matrix=FALSE){
    require(nplr)
    require(plyr)
    if(recalculate){
        ##get original data points per well
        cpd_data<-synGet('syn5632189')@filePath
        cpd_dat=read.table(cpd_data,sep='\t',header=T,quote='"')


        res.auc=dlply(cpd_dat,c("experiment_id","master_cpd_id"),function(dat){
            res=NA
            try(res<-nplr(x=dat$cpd_conc_umol,y=2^dat$bsub_value_log2/max(2^dat$bsub_value_log2))@AUC[[1]])
            return(AUC=res)
        })

        res.auc=data.frame(attr(res.auc,'split_labels'),unlist(res.auc))
      }else{
        res.auc<-read.table(synGet('syn5622708')@filePath,header=T)
    }
  ##now add cell and drug names
  new.df<-data.frame(AUC=as.numeric(res.auc$unlist.res.auc),
                     Drug=drug_dat$cpd_name[match(res.auc$master_cpd_id,drug_dat$master_cpd_id)],
                     CCL_id=ccl_dat$master_ccl_id[match(res.auc$experiment_id,ccl_dat$experiment_id)])
  new.df$Cell=ccl_metadat$ccl_name[match(new.df$CCL_id,ccl_metadat$master_ccl_id)]
  
  if(as.matrix)
    return(acast(new.df,Drug~Cell,value.var='AUC',fun.aggregate=function(x) mean(x,na.rm=T)))
  else
    return(new.df)
}

#' Get CTRP- predicted target of drug
#' @return data frame of drug and target
ctrpDrugTargets<-function(){
    tlist=apply(drug_dat,1,function(x){
        targs=unlist(strsplit(x[['gene_symbol_of_protein_target']],split=';'))
        drug=rep(x[['cpd_name']],length(targs))
        data.frame(Drug=drug,Target=targs)
    })
    tdat=do.call("rbind",tlist)
    return(tdat)
}


getMSSMprocessedMutationCalls<-function(){
  ##downlods file from synapse, then returns matrix
  tab=read.table(synGet('syn5706505')@filePath,header=T,comment='',as.is=T)
  mat=tab[-c(1:2),-c(1:3)]
  mat<-apply(mat,2,function(x) as.numeric(unlist(x)))
  rownames(mat)<-tab[-c(1:2),1]
  return(mat)
  
}

getMSSMTissueOrigin<-function(){
  tab=read.table(synGet('syn5706505')@filePath,header=T,comment='',as.is=T)
  mat=tab[1,-c(1:3)]
  #names(mat)<-colnames(mat)
  #mat<-apply(mat,2,function(x) as.numeric(unlist(x)))
  #rownames(mat)<-tab[-c(1:2),1]
  return(mat)
  
}
