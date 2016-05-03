##
## Drug sensitivity data files for NCATS single agent screens
##
##
library(synapseClient)
library(data.table)
library(pheatmap)
synapseLogin()
fileparent='syn5522627'

##download files
headerfile='syn5522652'
qr<-synapseQuery(paste("select * from entity where parentId=='",fileparent,"'",sep=''))
hind=which(qr$entity.id==headerfile)
header=read.table(synGet(qr[hind,'entity.id'])@filePath,sep='-',fill=T)
dfiles<-qr[grep('csv',qr$entity.name),]

##first read in all files and update curve class
allfiles<-lapply(dfiles$entity.id,function(x) {
  #print(x)
  res=read.table(synGet(x)@filePath,sep=',',header=T,fill=T,quote='"')
  #print(ncol(res))
  if('CRC'%in%colnames(res))
    cl=res$CRC
  else
    cl=res$CCLASS2
  res$CurveClass<-rep("Inconclusive",length(cl))
  res$CurveClass[which(cl%in%c(1.1, 1.2, 2.1, 2.2))]<-'Enhances'
  res$CurveClass[which(cl%in%c(-1.1, -1.2, -2.1, -2.2))]<-'Inhibits'
  res$CurveClass[which(cl==4)]<-'Inactive'
  return(res)
})

names(allfiles)<-dfiles$entity.sampleName

valsOfInterest<<-c(names(allfiles[[1]])[4:13],c('target','CurveClass'))

#' For each drug, collect a value and put into data frame
#' @param valname Name of value to extract, must be one of the collected values
#' @return A data frame containing a value for each NCAT-screened sample for each drug.
getValueForAllCells<-function(valname){
  #read in data
  if(!valname%in%valsOfInterest){
    print(paste("Value should be one of",paste(valsOfInterest,collapse=',')))
    return(NULL)
  }else{
    print(paste('Computing',valname,'for all cells'))
  }

  #get all drugs
  #ssalldrugs<-unique(allfiles[[1]]$name)
  drug.values<-sapply(dfiles$entity.sampleName,function(x){
      if(valname=='CRC' && !valname%in%names(allfiles[[x]]))
        valname='CCLASS2'
      vals<-allfiles[[x]][[valname]]
      drugs<-allfiles[[x]]$name
      names(vals)<-drugs
      vals
  })
  
  rownames(drug.values)[which(is.na(rownames(drug.values)))]<-'NA'
  dupes=rownames(drug.values)[which(duplicated(rownames(drug.values)))]
  rownames(drug.values)[which(duplicated(rownames(drug.values)))]<-paste(dupes,'alt',sep='_')
  return(drug.values)
}


#' Plots PCA plot for a particular cell of interest
#' Factors in specific features
#' @param cellname Name of cell type to plot
#' @param as.categ If true, use numeric curve class as categorical variable
#' @param use.disc If true, filter curve class further to be 'enhances', 'inhibits','inactive' or 'inconclusive'
plotOneCell<-function(cellname,as.categ=FALSE,use.disc=FALSE){
  td<-allfiles[[cellname]]
  #ntd<-td[,c(5:13,15:36)]
  ntd<-td[,5:13]
  zv<-which(apply(ntd,1,function(x) any(is.na(x))))
  nztd<-ntd[-zv,]
  pc=prcomp(nztd,scale.=TRUE)
  png(paste('CRC',ifelse(as.categ,'cat',''),
            ifelse(use.disc,'discrete',''),'ValuesByDrugFor',cellname,'.png',sep=''))

  if('CRC' %in% names(td))
    cl=td$CRC[-zv]
  else
    cl=td$CCLASS2[-zv]
  if(as.categ){
    cl=as.factor(cl)
  }else if(use.disc){
    cl=as.factor(td$CurveClass[-zv])
  }

  require(ggbiplot)
  p<-ggbiplot(pc,groups=cl)+ggtitle(paste('Drug response panel for',cellname))
  print(p)
  dev.off()
}

#' cluster drugs by response across cell lines
#' DEPRACATED - use
#' @param metric Metric to use to cluster drugs
#' @param h Height at which to divide clusters
#' @return data frame of each drug and each cluster
clusterDrugsByResponse<-function(metric='FAUC',h=4,doZscore=TRUE){
  #get fauc
    fauc=getValueForAllCells(metric)
    script.dir <- dirname(sys.frame(1)$ofile)
    source(paste(script.dir,'singleDrugAnalysis.R',sep='/'))
    return(getDrugClusters(fauc,h,doZScore))

}


#' Plot most variable values across all cells
#' @param valname
#' @param ft Format of file
#' @return drug values
plotMostVariableVals<-function(valname,ft='png'){
    drug.values<-getValueForAllCells(valname)
   if(valname=='CurveClass'){
      av<-drug.values[which(apply(drug.values,1,function(x) length(which(x=='Inhibits'))==5)),]
      levs<-as.factor(av)$levels
      mv<-apply(av,2,function(x) as.numeric(factor(x,levels=levs)))
      colnames(mv)<-colnames(av)
      rownames(mv)<-rownames(av)
    }else{
      drug.var<-apply(drug.values,1,var)
      mv<-drug.values[order(drug.var,decreasing=T)[1:50],]
     }

    gt<-dfiles$entity.sampleGenotype
    names(gt)<-dfiles$entity.sampleName
    pheatmap(t(mv),annotation_row=data.frame(Genotype=gt),cellwidth=10,cellheight=10,file=paste('drugsWithMostVariable',valname,'AcrossCellLines.',ft,sep=''))
    return(drug.values)

}



#' Re-compute drug response curves using the nplr package
#' @param cell Name of cell type to select
#' @param drug Name of drug to plot. If NA, will compute fits for all drugs
#' @param doPlot Will produce a plot of the fitted curve from NPLR
#' @param recalculate Set to true to recalculate (when drug is NA)
#' @return Data frame of drug and AUC value
doseResponseCurve<-function(cell,drug=NA,doPlot=TRUE){
    cell.resp=allfiles[[cell]]

    if(is.na(drug)){
            all.aucs=apply(cell.resp,1,function(drug.dat){
                dvals<-as.numeric(drug.dat[grep('DATA[0-9+]',names(drug.dat))])

                cvals<-as.numeric(drug.dat[grep('^C[0-9+]',names(drug.dat))])
                fit<-NA
            try(
                fit<-nplr(cvals,dvals/max(dvals),useLog=TRUE)@AUC[[1]]
            )

            return(fit)
            })
            names(all.aucs)<-cell.resp$name
            df=data.frame(Cell=rep(cell,length(all.aucs)),
                Drug=cell.resp$name,AUC=all.aucs)

    }
    else{
        drug.dat<-cell.resp[match(drug,cell.resp$name),]
        ##actual points
        dvals<-unlist(drug.dat[grep('DATA[0-9+]',names(drug.dat))])
        cvals<-unlist(drug.dat[grep('^C[0-9+]',names(drug.dat))])

  ##try to fit new model
        require(nplr)
        fit=NA
        try(fit<-nplr(cvals,dvals/max(dvals),useLog=TRUE)@AUC)
        if(doPlot){
            pdf(paste(drug,'doseResponseCurveIn',cell,'.pdf',sep=''))
            plot(fit)
            dev.off()
        }
        all.aucs=fit
        names(all.aucs)<-drug
        df=data.frame(Cell=cell,Drug=drug,AUC=all.aucs)
    }
  return(df)

}

getRecalculatedAUCTab<-function(){
  tab<-read.table(synGet('syn5637634')@filePath,header=T)
  return(tab)
}

getRecalculatedAUCMatrix<-function(){
    tab<-read.table(synGet('syn5637634')@filePath,header=T)
    require(reshape2)
    dmat=acast(tab,Drug~Cell,value.var='AUC',fun.aggregate=function(x) mean(x,na.rm=T))
    
    return(dmat)
}

#' Get NCATS- predicted target of drug
#' @return data frame of drug and target
ncatsDrugTargets<-function(){
    drugTargets<-getValueForAllCells('target')
    drugs=names(drugTargets[,1])
    targs=lapply(as.character(drugTargets[,1]),function(x) unlist(strsplit(x,split=', ')))
    names(targs)<-drugs
    
    fdf=do.call("rbind",lapply(drugs,function(x) data.frame(Drug=rep(x,length(targs[[x]])),Target=targs[[x]])))
    return(unique(fdf))
    
}
