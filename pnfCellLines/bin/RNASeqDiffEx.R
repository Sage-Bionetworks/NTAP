##now we want to test diffex
library('sleuth')
library('synapseClient')
synapseLogin()


##this function will download and construct the data frame used to build
##the sleuth object
##includes all the variables needed
buildDiffExDf<-function(){

  allfiles<-synapseQuery("select name,id from entity where parentId=='syn5579785'")
  h5files=allfiles[grep('*h5',allfiles[,1]),]

  h5s=lapply(h5files[,2],function(x) synGet(x))
  snames=lapply(h5s,function(x) x@annotations$sampleName)
  sgens=lapply(h5s,function(x) x@annotations$sampleGenotype)
  sors=lapply(h5s,function(x) x@annotations$sampleOrigin)
  sex=rep('Male',length(sgens))
  sex[which(sgens=='++')]<-'Female'
  filepaths=lapply(h5s,function(x) x@filePath)


  ##build the contrast matrix from the annotations
  df=data.frame(sample=unlist(snames),path=unlist(filepaths),Sex=factor(sex),Origin=factor(unlist(sors),levels=rev(unique(unlist(sors)))),Genotype=factor(unlist(sgens),levels=c("++","+-","--")))
  one.all=rep('+',length(sgens))
  one.all[which(sgens=='--')]<-'-'
  df$OneAllele=factor(one.all,levels=c('-','+'))
  cul=rep('primary',nrow(df))
  cul[grep("^i",df$sample)]<-'immortalized'
  df$Culture=cul
  df$path=as.character(df$path)
  return(df)
}

buildGencodeTargetMap<-function(df){

  sf<-sleuth_prep(df,~ OneAllele + Culture + Origin)

  ##now get gene names
  t2g<-t(sapply(unique(sf$obs_raw$target_id),function(x) {
    res=unlist(strsplit(x,split='|',fixed=T))
    list(transcript=res[1],gene=res[6])}))
  t2g<-data.frame(target_id=unique(sf$obs_raw$target_id),t2g)
  t2g<-apply(t2g,2,unlist)
  write.table(t2g,file='../../data/gencodeGeneTranscriptMap.csv',sep=',',col.names=T,row.names=F)
}

##now we need to build the sleuth object
buildSleuthModel<-function(df,inc=c('Sex','Culture'),test='OneAllele+'){

  t2g2<-read.table('../../data/gencodeGeneTranscriptMap.csv',sep=',',header=T)
  ##first build sleuth object
  formstring='~ OneAllele'
  ex=paste(inc,collapse=' + ')
  if(ex!="")
    formstring=paste(formstring,ex,sep=' + ')

  sf<-sleuth_prep(df,as.formula(formstring),target_mapping=t2g2)

  #now do fit
  ffit<-sleuth_fit(sf)

  ##now ask question
  fp <- sleuth_wt(ffit, test)

  return(fp)
}

getSleuthTable<-function(sleuthfit){
  #write results to table
  res.tab<-sleuth_results(sleuthfit,'OneAllele+')
  rtab<-apply(res.tab,2,unlist)
  rtab
}

getGOList<-function(sleuthfit){


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
  genotype=rep("-",nrow(dm))
  names(genotype)<-rownames(dm)
  #genotype[which(dm[,'Genotype+-']==1)]<-'+-'
  genotype[which(dm[,'OneAllele+']==1)]<-'+'


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

plotVals<-function(data.obj,qval,ttype=c(),prefix=''){
  res=sleuth_results(data.obj, 'OneAllele+')
  sel=which(res$qval<qval)
  print(paste("Found",length(sel),'diff ex transcripts at q=',qval))

  targs=as.character(res$target_id)#[sel])
  trans.type=sapply(as.character(targs),function(x) {
    arr=unlist(strsplit(x,split='|',fixed=T))
    arr[length(arr)]
  })

  if(length(ttype)>0){
    sel=intersect(sel,which(trans.type%in%ttype))
    print(length(sel))
  }
  targs=targs[sel]
  gene=res$gene[sel]
  trans=res$transcript[sel]
  tnames=paste(gene,trans,sep='_')

  names(tnames)<-targs
  names(trans.type)<-tnames
  plotGenesInSamples(data.obj,targs,'tpm',tnames,trans.type[sel],fname=paste(prefix,'diffex',paste(ttype,collapse='_'),'TranscriptsQ',qval,'.png',sep=''))


}
