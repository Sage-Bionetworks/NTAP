##now we want to test diffex
library('sleuth')
library('synapseClient')
synapseLogin()


source("../../bin/runDiffEx.R")

fulldf<-buildDiffExDf()



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
