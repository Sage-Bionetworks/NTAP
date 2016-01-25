library(CNTools)
library(ggbiplot)
library(synapseClient)
library(pheatmap)

#first log into synapse
synapseLogin()


source("../../bin/CNVData.R")

##read in segmentation data
if(!exists("segdata"))
    segdata <- tier1_segmentedData(FALSE)
if(!exists("segdata2"))
    segdata2 <- tier1_segmentedData(TRUE)

ann.vals<-synTableQuery('SELECT "Sample Name","Sample ID","Sample Genotype" FROM syn5014742')@values

idDiffVals<-function(segdat,byval='gene',metric='median',thresh=-.5,byChrom=FALSE){
  require(dplyr)
  cnseg <- CNSeg(segdat)
  rdseg <- getRS(cnseg, by = byval,geneMap=geneInfo, imput = FALSE, XY = FALSE, what =metric)
  segM <- rs(rdseg)
  nzvals<-which(apply(segM,1,function(x) any(as.numeric(x[-c(1:5)])<thresh)))
  nzM<-segM[nzvals,]
  print(paste('found',length(nzvals),byval,'values with logR values less than',thresh))
   
  vals_by_gene=as.data.frame(dlply(nzM,'genename',function(x) apply(x[,-c(1:5)],2,median)))
  ##lets get annotations
  gts=ann.vals[match(rownames(vals_by_gene),ann.vals[,2]),3]
  rownames(vals_by_gene)<-ann.vals[match(rownames(vals_by_gene),ann.vals[,2]),1]
  names(gts)<-rownames(vals_by_gene)
  
  nongenes=c(grep("RNA",colnames(vals_by_gene)),grep("MIR",colnames(vals_by_gene)),
             grep("SNO",colnames(vals_by_gene)),grep("HV",colnames(vals_by_gene)),
             grep("C[0-9]+orf",colnames(vals_by_gene)),grep("BC0",colnames(vals_by_gene)),
             grep("LOC[0-9]+",colnames(vals_by_gene)))
  if(length(nongenes)>0)
    vals_by_gene=vals_by_gene[,-nongenes]
  print(paste(ncol(vals_by_gene),'unique gene-specific values'))
  ##genes by chromosome
  res=nzM %>% group_by(chrom) %>% summarise(genes=n_distinct(genename))

  ##heatmaps by chromosome
  nm=apply(vals_by_gene,2,as.numeric)
  rownames(nm)<-rownames(vals_by_gene)
  chroms=nzM$chrom
  names(chroms)<-nzM$genename
  chroms=chroms[unique(names(chroms))]
  nm=t(nm)[order(chroms[colnames(nm)]),]
  if(byChrom){
    sapply(chroms,function(y){
      gvals=intersect(names(chroms)[which(chroms==y)],rownames(nm))
      pheatmap(nm[gvals,],cellheight=10,cellwidth=10,cluster_rows=TRUE,
               annotation_col=data.frame(Genotype=gts),annotation_row=data.frame(Chromosome=chroms),
               file=paste(metric,byval,'valuesLessThan',thresh,'OnChrom',y,'.png',sep=''))
    })
  }else{
    pheatmap(nm,cellheight=10,cellwidth=10,cluster_rows=FALSE,
             annotation_col=data.frame(Genotype=gts),annotation_row=data.frame(Chromosome=chroms),
             file=paste(metric,byval,'valuesLessThan',thresh,'.png',sep=''))
    }
  return(nm)
}


####Plot clusters - assumes we're in analysis directory!!
if(!exists("geneInfo"))
    geneInfo<-read.table('../../data/hg19_geneInfo.txt')

plotClusteredSegmentData<-function(segdat,byval='gene',metric='median',topGenes=100,prefix=''){

    print(paste("Preparing to analyze segmented CNV data by",metric,'focusing on top',topGenes,byval))
    ##first get data matrix of regions by gene
    cnseg <- CNSeg(segdat)
    if(byval=='gene')
        rdseg <- getRS(cnseg, by = byval,geneMap=geneInfo, imput = FALSE, XY = FALSE, what =metric)
    else
        rdseg <- getRS(cnseg, by = byval, imput = FALSE, XY = FALSE, what = metric)

    segM <- rs(rdseg)

    col.inds=4:ncol(segM)
    if(byval=='gene')
        col.inds=6:ncol(segM)
    ##create data matrix from segM output
    M <- t(do.call("rbind", lapply(col.inds, function(i) as.numeric(as.character(segM[,i])))))
#    idxs <- match(colnames(segM)[col.inds], sample.names)
    colnames(M)<-colnames(segM)[col.inds]
    if(byval=='gene'){
        rownames(M)<-segM[,5]
        ##reduce matrix by individual genes
        print("Reducing matrix to individual genes")
        redM<-sapply(unique(rownames(M)),function(x){
            mm<-which(rownames(M)==x)
            if(length(mm)==1)
                return(M[mm,])
            else
                return(colMeans(M[mm,]))})
        M<-t(redM)

    }
                                            #filter out zero variance
    allvars<-apply(M,1,var,na.rm=T)
    if(length(which(allvars==0))>0){
        M<-M[-which(allvars==0),]
        allvars=allvars[-which(allvars==0)]
    }


    ##ANalysis 0: cluster ALL gene/regions
    fname=paste('all',byval,'by',metric,'logRRatios.pdf',sep='_')
    pdf(fname)
    tm<-M
    pc=prcomp(t(tm),scale.=T,center.=T)

    gb<-ggbiplot(pc,groups=as.factor(origin[colnames(M)]),choices=1:2,obs.scale=1,var.scale=1,ellipse=TRUE,var.axes=FALSE) +
        scale_color_discrete(name = '') +
            theme(legend.direction = 'horizontal', legend.position = 'top')
    print(gb)

    dev.off()
    colnames(tm) <- clnames[colnames(tm)]

    fname=paste('all',byval,'by',metric,'logRRatios_dendro.png',sep='_')
    png(fname)
    plot(hclust(dist(t(tm)),method="ward.D2"),main="Clustering of samples based on distance")
    dev.off()



    ##Analysis 1: how do top 100 most variable genes cluster?

    tm=M[order(allvars,decreasing=T)[1:100],]
                                        #print(head(M))
    if(byval=='gene')
        nfcor=apply(tm,1,function(x) cor(x,M['NF1',]))
    else
        nfcor=rep(0,nrow(M))

    fname=paste('most_variable',topGenes,byval,'by',metric,'logRRatios.pdf',sep='_')
    pdf(fname)
    pc=prcomp(t(tm),scale.=T,center.=T)
    gb<-ggbiplot(pc,groups=as.factor(origin[colnames(M)]),choices=1:2,obs.scale=1,var.scale=1,ellipse=TRUE,var.axes=FALSE) +
        scale_color_discrete(name = '') +
            theme(legend.direction = 'horizontal', legend.position = 'top')
    print(gb)
   # colnames(tm) <- clnames[colnames(tm)]

    dev.off()
    fname=paste('most_variable',topGenes,byval,'by',metric,'logRRatios_dendro.png',sep='_')
    png(fname)
    plot(hclust(dist(t(tm)),method="ward.D2"),main="Clustering of samples based on distance")

    dev.off()


    fname=paste('most_variable',topGenes,byval,'by',metric,'logRRatios_heatmap.png',sep='_')
    if(byval=='gene')
        pheatmap(t(tm),annotation_col=data.frame(NF1Corr=nfcor),
                 annotation_row=data.frame(Genotype=genotype),cellwidth=10,cellheight=10,file=fname)
    else
        pheatmap(t(tm),annotation_row=data.frame(Genotype=genotype),cellwidth=10,cellheight=10,file=fname)
    ## Analysis 2
    ##let's do supervised clustering, look for cnv values that correlate with genotype
    gcors=apply(M,1,function(x) cor(x,as.numeric(as.factor(genotype[colnames(tm)]))))

                                        #or regression p-values
                                        #let's rank all genes by absolute correlation for future functional enrichment!!!!
    write(names(gcors)[order(gcors,decreasing=T)],file='geneRankByAbsCorGenotype.txt')
                                        #sf=File(filename,parentId='syn5014748')
                                        #used(sf)<-synGet('syn5015035')
                                        #synStore(sf)

                                        #now let's take top 100 genes and plot/cluster again
    tm=M[order(gcors,decreasing=T)[1:100],]
    colnames(tm) <- clnames[colnames(tm)]
    

    fname=paste('most_gt_correlated',topGenes,byval,'by',metric,'logRRatios.pdf',sep='_')
    pdf(fname)
    pc=prcomp(t(tm),scale.=T,center.=T)
    ggbiplot(pc,groups=as.factor(origin[colnames(tm)]),choices=1:2,obs.scale=1,var.scale=1,ellipse=TRUE,var.axes=FALSE) #+
    scale_color_discrete(name = '') #+
    theme(legend.direction = 'horizontal', legend.position = 'top')
    dev.off()

    fname=paste('most_gt_correlated',topGenes,byval,'by',metric,'logRRatios_dendro.png',sep='_')
    png(fname)
    colnames(tm) <- clnames[colnames(tm)]
    plot(hclust(dist(t(tm))),main="Clustering of samples based on distance")

    dev.off()
    fname=paste('most_gt_correlated',topGenes,byval,'by',metric,'logRRatios_heatmap.png',sep='_')

    if(byval=='gene')
        pheatmap(t(tm),annotation_col=data.frame(NF1Corr=nfcor),annotation_row=data.frame(Genotype=genotype),
                 cellwidth=10,cellheight=10,file=fname)
    else
        pheatmap(t(tm),annotation_row=data.frame(Genotype=genotype),
                 cellwidth=10,cellheight=10,file=fname)
}

#here are the primary tasks we want to do
main<-function(){
    plotClusteredSegmentData(segdata2,'gene','median',100)
  #  plotClusteredSegmentData(segdata2,'gene','median',1000)

    plotClusteredSegmentData(segdata2,'region','median',100)
  #  plotClusteredSegmentData(segdata2,'region','median',1000)

}
