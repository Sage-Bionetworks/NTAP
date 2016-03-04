##now we want to test diffex
library('sleuth')
library('synapseClient')
synapseLogin()


##this function will download and construct the data frame used to build
##the sleuth object
##includes all the variables needed
#' @param outliers - list of cell line names to expcloe
buildDiffExDf<-function(outliers=c()){

  allfiles<-synapseQuery("select name,id,sampleName from entity where parentId=='syn5579785'")
  h5files=allfiles[grep('*h5',allfiles[,1]),]
  
  ol=which(h5files$entity.sampleName%in%outliers)
  if(length(ol)>0)
    h5files=h5files[-ol,]

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
buildSleuthModel<-function(df,inc=c('Sex','Culture'),test='OneAllele',alt='+',prefix=''){

  t2g2<-read.table('../../data/gencodeGeneTranscriptMap.csv',sep=',',header=T)
  ##first build sleuth object
  formstring=paste('~',test)
  ex=paste(inc,collapse=' + ')
  if(ex!="")
    formstring=paste(formstring,ex,sep=' + ')

  print(paste('Building sleuth model for',formstring))
  sf<-sleuth_prep(df,as.formula(formstring),target_mapping=t2g2)

  #now do fit
  ffit<-sleuth_fit(sf)

  ##now ask question
  fp <- sleuth_wt(ffit, paste(test,alt,sep=''))

  return(fp)
}

getSleuthTable<-function(sleuthfit,test,alt){
  #write results to table
  res.tab<-sleuth_results(sleuthfit,paste(test,alt,sep=''))
  rtab<-apply(res.tab,2,unlist)
  rtab
}

getGOList<-function(sleuthfit){


}


plotGenesInSamples<-function(obj,transcripts,units="tpm",
                             genes=NULL,annotes=NULL,fname='selectedGenesInData.pdf',
                             test='OneAllele',alt='+',collapseByGene=FALSE){

    jm = obj$obs_norm
    jm$gene = obj$target_mapping$gene[match(jm$target_id,obj$target_mapping$target_id)]

    tabd_df <- jm[jm$target_id %in% transcripts,]

  ##select unit to plo
  if (units == "tpm") {
    tabd_df <- dplyr::select(tabd_df, gene,target_id, sample,
                             tpm)
    if(collapseByGene)
        tabd_df <- reshape2::dcast(tabd_df, gene ~ sample,
                                   value.var = "tpm", fun.aggregate=sum)
   else
        tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample,
                                   value.var = "tpm")
  }
  else if (units == "est_counts") {
    tabd_df <- dplyr::select(tabd_df, gene,target_id, sample,
                             est_counts)
    if(collapseByGene)
        tabd_df <- reshape2::dcast(tabd_df, gene ~ sample,
                                   value.var = "est_counts", fun.aggregate=sum)
    else
        tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample,
                                   value.var = "est_counts")
  }
  else {
    stop("Didn't recognize the following unit: ", units)
  }

  dm=design_matrix(obj)
  cul=rep('Immortalized',nrow(dm))
  names(cul)<-rownames(dm)
  if('Cultureprimary'%in%colnames(dm))
    cul[which(dm[,'Cultureprimary']==1)]<-'Primary'
  else
    cul=NULL

  if(test=='OneAllele'){
        oa=rep("-",nrow(dm))
        names(oa)<-rownames(dm)
        oa[which(dm[,'OneAllele+']==1)]<-'+'
        adf=data.frame(OneNF1=oa)
  }else if(test%in%c('Genotype','Culture')){

      genotype=rep("++",nrow(dm))
      names(genotype)<-rownames(dm)
      genotype[which(dm[,'Genotype+-']==1)]<-'+-'
      genotype[which(dm[,'Genotype--']==1)]<-'--'
      adf=data.frame(Genotype=genotype)
  }
  if(!is.null(cul))
    adf$Culture=cul




    if(!is.null(genes))
        rownames(tabd_df) <- genes[tabd_df[,1]]
    else
        rownames(tabd_df)=tabd_df[,1]

    if(!is.null(annotes))
        names(annotes)<-rownames(tabd_df)


  require(pheatmap)
  if(is.null(annotes))
      pheatmap(t(log2(tabd_df[,-1]+0.01)),cellheight=10,cellwidth=10,
               annotation_row=adf,clustering_distance_rows='correlation',
               clustering_distance_cols='correlation',filename=fname)
  else
      pheatmap(t(log2(tabd_df[,-1]+0.01)),cellheight=10,cellwidth=10,
               annotation_row=adf,clustering_distance_rows='correlation',
               clustering_distance_cols='correlation',
               annotation_col=data.frame(Type=annotes),filename=fname)

}


plotSingleGene<-function(obj,gene,prefix=''){
    library(ggplot2)
    library(plyr)
    jm <-obj$obs_norm
    jm$gene = obj$target_mapping$gene[match(jm$target_id,obj$target_mapping$target_id)]


    tabd_df <- jm[which(jm$gene==gene),]

    tsum=daply(tabd_df,'sample', function(td){sum(td$tpm)})
    newdf<-data.frame(TPM=tsum,Sample=names(tsum))
    newdf$Genotype<-obj$sample_to_covariates$Genotype[match(names(tsum),obj$sample_to_covariates$sample)]
    newdf$Culture<-obj$sample_to_covariates$Culture[match(names(tsum),obj$sample_to_covariates$sample)]
    newdf$Origin<-obj$sample_to_covariates$Origin[match(names(tsum),obj$sample_to_covariates$sample)]

    pdf(paste(prefix,gene,'expression.pdf',sep=''))
    p=ggplot(newdf)+geom_bar(aes(Origin,TPM,fill=Culture),stat='identity',position='dodge')
    print(p)
    dev.off()
}

plotVals<-function(data.obj,qval,ttype=c(),prefix='',test='OneAllele',alt='+'){
  res=sleuth_results(data.obj, paste(test,alt,sep=''))
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

#general analysis
makeTablesAndPlots<-function(model, test,alt,prefix=''){
  fname=paste(test,'variable',alt,'test',prefix,sep='_')
  ##1
  tab<-getSleuthTable(model,test,alt)
  #store table
  write.table(tab,file=paste(fname,'SleuthResults.csv',sep=''),sep=',',quote=F)
  sigs<-which(as.numeric(tab[,'qval'])<0.1)
  
  pc=sigs[grep('protein_coding',tab[sigs,'target_id'])]
  
  pcg=unique(tab[pc,'gene'])
  print(paste("Found",length(sigs),'differentially expressed transcripts, of which',length(pc),'are protein coding representing',length(pcg),'unique genes'))
  
  pdf(paste(fname,'VolcanoPlot.pdf',sep=''))
  res=plot_volcano(model,paste(test,alt,sep=''),point_alpha=0.8)
  print(res)
  dev.off()
  
  pdf(paste(fname,'PCAPlot.pdf',sep=''))
  res=plot_pca(model,text_labels=TRUE,color_by=test,point_size=8)
  print(res)
  dev.off()
  
  gvals=tab[,'gene']
  names(gvals)<-tab[,'target_id']
  plotGenesInSamples(model,tab[pc,'target_id'],units="tpm",
                     genes=NULL,annotes=NULL,collapseByGene=TRUE,
                     fname=paste(fname,'sigGenesByTpmInHeatmap.pdf'),test=test,alt=alt)
  
  
}
