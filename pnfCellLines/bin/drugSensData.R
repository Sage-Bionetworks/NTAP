##
## Drug sensitivity data file
##
##
##
library(synapseClient)
library(data.table)
synapseLogin()
fileparent='syn5522627'

##download files
headerfile='syn5522652'
qr<-synapseQuery(paste("select * from entity where parentId=='",fileparent,"'",sep=''))
hind=which(qr$entity.id==headerfile)
header=read.table(synGet(qr[hind,'entity.id'])@filePath,sep='-',fill=T)
dfiles<-qr[-hind,]

valsOfInterest<-names(allfiles[[1]])[4:13]

##now for each drug, collect a value and put into data frame
getValueForAllCells<-function(valname){
  if(!valname%in%valsOfInterest){
    print(paste("Value should be one of",paste(valsOfInterest,collapse=',')))
    return(NULL)
  }else{
    print(paste('Computing',valname,'for all cells'))
  }
    
  #read in data
  allfiles<-lapply(dfiles$entity.id,function(x) as.data.frame(fread(synGet(x)@filePath,sep=',',header=T)))
  names(allfiles)<-dfiles$entity.sampleName
  
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

plotOneCell<-function(cellname){
  td<-allfiles[[cellname]]
  #ntd<-td[,c(5:13,15:36)]
  ntd<-td[,5:13]
  zv<-which(apply(ntd,1,function(x) any(is.na(x))))
  nztd<-ntd[-zv,]
  pc=prcomp(nztd,scale.=TRUE)
  png(paste('CRCValuesByDrugFor',cellname,'.png',sep=''))
  if('CRC' %in% names(td))
    cl=td$CRC[-zv]
  else
    cl=td$CCLASS2[-zv]
    p<-ggbiplot(pc,groups=cl)+ggtitle(paste('Drug response panel for',cellname))
  print(p)
  dev.off()
}

plotMostVariableVals<-function(valname){
    drug.values<-getValueForAllCells(valname)
    
    #navals<-which(apply(drug.values,1,function(x) any(is.na(x))))
    #ddv<-drug.values[-navals,]
    
    drug.var<-apply(drug.values,1,var)
    mv<-drug.values[order(drug.var,decreasing=T)[1:50],]
    gt<-dfiles$entity.sampleGenotype
    names(gt)<-dfiles$entity.sampleName
    pheatmap(t(mv),annotation_row=data.frame(Genotype=gt),cellwidth=10,cellheight=10,file=paste('drugsWithMostVariable',valname,'AcrossCellLines.png',sep=''))
  return(drug.values)  
}