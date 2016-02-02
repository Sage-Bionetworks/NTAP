##
## Drug sensitivity data file
##
##
##
library(synapseClient)
library(data.table)
library(ggbiplot)
library(pheatmap)
synapseLogin()
fileparent='syn5522627'

##download files
headerfile='syn5522652'
qr<-synapseQuery(paste("select * from entity where parentId=='",fileparent,"'",sep=''))
hind=which(qr$entity.id==headerfile)
header=read.table(synGet(qr[hind,'entity.id'])@filePath,sep='-',fill=T)
dfiles<-qr[-hind,]

##first read in all files and update curve class
allfiles<-lapply(dfiles$entity.id,function(x) {
  #print(x)
  res=as.data.frame(fread(synGet(x)@filePath,sep=',',header=T))
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

##now for each drug, collect a value and put into data frame
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
  return(drug.values)
}

plotOneCell<-function(cellname,as.categ=FALSE,use.disc=FALSE){
  td<-allfiles[[cellname]]
  #ntd<-td[,c(5:13,15:36)]
  ntd<-td[,5:13]
  zv<-which(apply(ntd,1,function(x) any(is.na(x))))
  nztd<-ntd[-zv,]
  pc=prcomp(nztd,scale.=TRUE)
  png(paste('CRC',ifelse(as.categ,'cat',''),ifelse(use.disc,'discrete',''),'ValuesByDrugFor',cellname,'.png',sep=''))
  if('CRC' %in% names(td))
    cl=td$CRC[-zv]
  else
    cl=td$CCLASS2[-zv]
  if(as.categ){
    cl=as.factor(cl)
  }else if(use.disc){
    cl=as.factor(td$CurveClass[-zv])
  }
  
  p<-ggbiplot(pc,groups=cl)+ggtitle(paste('Drug response panel for',cellname))
  print(p)
  dev.off()
}

##cluster drugs by response across cell lines
clusterDrugsByResponse<-function(metric='FAUC',h=4){
  #get fauc
  fauc=getValueForAllCells(metric)
  #now zscore them
  zsfauc<-apply(fauc,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
  nz.zs.fauc=zsfauc
  #reset NA values
  nz.zs.fauc[which(is.na(nz.zs.fauc),arr.ind=T)]<-0.0
  #pheatmap(nz.zs.fauc,cellheight=10,cellwidth=10,file='Zscore_FAUC_for_all_cells.pdf')
  
  ##there is a pretty blue cluster there, can we do any enrichment? 
  drug.dists<-dist(nz.zs.fauc)
  hc=hclust(drug.dists)
  
  ##now cut the clustering'
  drug.clusters<-sapply(unique(cutree(hc,h=h)),function(x) names(which(cutree(hc,h=h)==x)))
 # hist(sapply(drug.clusters,length))
  return(drug.clusters)
}

plotMostVariableVals<-function(valname,ft='png'){
    drug.values<-getValueForAllCells(valname)
    
    #navals<-which(apply(drug.values,1,function(x) any(is.na(x))))
    #ddv<-drug.values[-navals,]
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
library(ggplot2)

##need to show dose response curves.
doseResponseCurve<-function(cell,drug){
  cell.resp=allfiles[[cell]]
  drug.dat<-cell.resp[match(drug,cell.resp$name),]
  ##actual points
  dvals<-unlist(drug.dat[grep('DATA[0-9+]',names(drug.dat))])
  cvals<-unlist(drug.dat[grep('^C[0-9+]',names(drug.dat))])
  #drug.dat<-data.frame(substance=rep(drug,length(dvals)),dose=dvals,response=cvals,unit=rep("uM",length(dvals)))
#  fitvals=sapply(cvals,function(x) drug.dat$ZERO+((drug.dat$INF-drug.dat$ZERO)/(1+(log10(x)/drug.dat$LAC50))^(1*drug.dat$HILL)))
  ##try to fit new model
  require(nplr)
  fit=nplr(cvals,dvals/max(dvals),useLog=TRUE)
  pdf(paste(drug,'doseResponseCurveIn',cell,'.pdf',sep=''))
  plot(fit)
  dev.off()

}