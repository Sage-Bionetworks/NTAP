source("../../bin/ncatsCombinationScreens.R")

##test out plots by parameters -box plots
for(combo in c('6x6','10x10')){
  files<-getFileForCombo(combo,'CTG')
  for(val in c('Beta','DBSumNeg','DBSumPos','Gamma')){
    res<-plotValsAcrossCells(files,paste(combo,'CTG',sep=''),val)
    lms<-doLinearModel(res,paste(combo,'CTG',val,'Values',sep=''))
  }
}
combo='6x6'
files<-getFileForCombo(combo,'CTG')
row.targs<-unique(unlist(lapply(files,function(x) return(x$RowTarget))))
col.targs<-unique(unlist(lapply(files,function(x) return(x$ColTarget))))
write.table(data.frame(Row=row.targs,Col=col.targs),file='NCATS_targets_6x6.tsv',sep='\t',row.names=F,col.names=T)

combo='10x10'
files<-getFileForCombo(combo,'CTG')
row.targs<-unique(unlist(lapply(files,function(x) return(x$RowTarget))))
col.targs<-unique(unlist(lapply(files,function(x) return(x$ColTarget))))
write.table(data.frame(Row=row.targs,Col=col.targs),file='NCATS_targets_10x10.tsv',sep='\t',row.names=F,col.names=T)


combo='6x6'
files<-getFileForCombo(combo,'CTG')
row.targs<-unique(unlist(lapply(files,function(x) return(x$RowName))))
col.targs<-unique(unlist(lapply(files,function(x) return(x$ColName))))
write.table(data.frame(Row=row.targs,Col=col.targs),file='NCATS_drugs_6x6.tsv',sep='\t',row.names=F,col.names=T)

combo='10x10'
files<-getFileForCombo(combo,'CTG')
row.targs<-unique(unlist(lapply(files,function(x) return(x$RowName))))
col.targs<-unique(unlist(lapply(files,function(x) return(x$ColName))))
write.table(data.frame(Row=row.targs,Col=col.targs),file='NCATS_drugs_10x10.tsv',sep='\t',row.names=F,col.names=T)


scripturl='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-05-09/testComboScreens.R'
combo.screen.sid='syn6042830'
fileparent=screendirs[['6x6']]  
cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
six.list<-lapply(cells[,1],function(x){
  synid=cells[match(x,cells[,1]),2]
  allf<-synapseQuery(paste("select id from entity where parentId=='",synid,"'",sep=''))
#  allf<-allf[grep(measure,allf[,1]),]
  return(allf)})

six.files=lapply(unique(unlist(six.list)),function(x) return(list(entity=x)))

fileparent=screendirs[['10x10']]  
cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
ten.list<-lapply(cells[,1],function(x){
  synid=cells[match(x,cells[,1]),2]
  allf<-synapseQuery(paste("select id from entity where parentId=='",synid,"'",sep=''))
  allf<-allf[grep('CTG',allf[,1]),]
  return(allf)})

ten.files=lapply(unique(unlist(six.list)),function(x) return(list(entity=x)))


allfiles<-list.files('./')
sixes<-allfiles[grep('6x6',allfiles)]
tens<-allfiles[grep('10x10',allfiles)]


ten.list=c()

for(f in sixes){
  synStore(File(f,parentId=combo.screen.sid),used=six.files,executed=list(list(url=scripturl)))
}


for(f in tens){
  synStore(File(f,parentId=combo.screen.sid),used=ten.files,executed=list(list(url=scripturl)))
}

##lastly get 6x6 and 10x10 drug combos


