##now we want to test diffex
library('sleuth')
library('synapseClient')
synapseLogin()


buildModels<-function(){
allfiles<-synapseQuery("select name,id from entity where parentId=='syn5579785'")
h5files=allfiles[grep('*h5',allfiles[,1]),]

h5s=lapply(h5files[,2],function(x) synGet(x))
snames=lapply(h5s,function(x) x@annotations$sampleName)
sgens=lapply(h5s,function(x) x@annotations$sampleGenotype)
filepaths=lapply(h5s,function(x) x@filePath)

##build the contrast matrix from the annotations
df=data.frame(sample=unlist(snames),path=unlist(filepaths),Genotype=factor(unlist(sgens),levels=c("++","+-","--")))
cul=rep('primary',nrow(df))
cul[grep("^i",df$sample)]<-'immortalized'
df$Culture=cul
df$path=as.character(df$path)

sf<-sleuth_prep(df,~ Genotype)

##now get gene names
t2g<-t(sapply(unique(sf$obs_raw$target_id),function(x) {
  res=unlist(strsplit(x,split='|',fixed=T))
  list(transcript=res[1],gene=res[6])}))
t2g<-data.frame(target_id=unique(sf$obs_raw$target_id),t2g)

inc.het=FALSE
if(inc.het){
  sf<-sleuth_prep(df,~ Genotype+Culture,target_mapping=t2g)
  ffit<-sleuth_fit(sf)
  fp <- sleuth_wt(ffit, 'Genotype--')

 res.tab<-sleuth_results(fp,'Genotype--')
 rtab<-apply(res.tab,2,unlist)
 write.table(rtab,file='sleuth_diffex_KO_vs_WT_inc_culture_het.csv',sep=',',row.names=F)
 write.table(summary(sf),file='kal_quant_summary.csv',sep=',',row.names=F,quote=F)
}

##is this any different?  YES - more genes diffex without heet
posneg=df[which(df$Genotype!='+-'),]
posneg$Genotype=factor(as.character(posneg$Genotype),levels=c("++",'--'))

so<-sleuth_prep(posneg,~ Genotype+Culture,target_mapping=t2g)
#write.table(summary(so),file='kal_posnegquant_summary.csv',sep=',',row.names=F,quote=F)

sfit<-sleuth_fit(so)
sp <- sleuth_wt(sfit, 'Genotype--')

res.tab.nohet <- sleuth_results(sp, 'Genotype--')
rtab.nh<-apply(res.tab.nohet,2,unlist)
write.table(rtab.nh,file='sleuth_diffex_KO_vs_WT_inc_culture.csv',sep=',',row.names=F)
}

plotGenesInSamples<-function(obj,transcripts,units="tpm", genes=NULL,annotes=NULL,fname='selectedGenesInData.png'){
  tabd_df <- obj$obs_norm[obj$obs_norm$target_id %in% transcripts,]
  ##select unit to plo
  if (units == "tpm") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             tpm)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "tpm")
  }
  else if (units == "est_counts") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             est_counts)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "est_counts")
  }
  else {
    stop("Didn't recognize the following unit: ", units)
  }
  
  dm=design_matrix(obj)
  genotype=rep("++",nrow(dm))
  names(genotype)<-rownames(dm)
  genotype[which(dm[,'Genotype+-']==1)]<-'+-'
  genotype[which(dm[,'Genotype--']==1)]<-'--'
  
  
  if(is.null(genes))
    genes=transcripts
  rownames(tabd_df) <- genes[tabd_df[,1]]
  names(annotes)<-genes[tabd_df[,1]]
  require(pheatmap)
  if(is.null(annotes))
    pheatmap(t(log2(tabd_df[,-1]+0.01)),cellheight=10,cellwidth=10,annotation_row=data.frame(Genotype=genotype),clustering_distance_rows='correlation',clustering_distance_cols='correlation',filename=fname)
  else
    pheatmap(t(log2(tabd_df[,-1]+0.01)),cellheight=10,cellwidth=10,annotation_row=data.frame(Genotype=genotype),clustering_distance_rows='correlation',clustering_distance_cols='correlation',annotation_col=data.frame(Type=annotes),filename=fname)
  
}

plotVals<-function(data.obj,qval){
  res=sleuth_results(data.obj, 'Genotype--')
  sel=which(res.tab.nohet$qval<qval)
  print(paste("Found",length(sel),'diff ex genes at q=',qval))
  targs=as.character(res.tab.nohet$target_id[sel])
  trans.type=sapply(as.character(targs),function(x) {
    arr=unlist(strsplit(x,split='|',fixed=T))
    arr[length(arr)]
  })
  gene=res$gene[sel]
  trans=res$transcript[sel]
  tnames=paste(genes,trans,sep='_')
  names(tnames)<-targs
  names(trans.type)<-tnames
  plotGenesInSamples(data.obj,targs,'tpm',tnames,trans.type,fname=paste('diffexGenesq',qval,'.png',sep=''))
  
  
}

plotVals(fp,0.05)
plotVals(fp,0.01)
plotVals(fp,0.1)