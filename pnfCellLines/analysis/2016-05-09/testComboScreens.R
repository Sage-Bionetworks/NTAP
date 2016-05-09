source("../../bin/ncatsCombinationScreens.R")

##test out plots by parameters -box plots
for(combo in c('6x6','10x10')){
  files<-getFileForCombo(combo,'CTG')
  for(val in c('Beta','DBSumNeg','DBSumPos','Gamma')){
    res<-plotValsAcrossCells(files,paste(combo,'CTG',sep=''),val)
    lms<-doLinearModel(res,paste(combo,'CTG',val,'Values',sep=''))
  }
}

scripturl=''
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
  synStore(File(f,parentId=combo.screen.sid),used=list(six.list),executed=list(list(url=scripturl)))
}


for(f in tens){
  synStore(File(f,parentId=combo.screen.sid),used=list(ten.list),executed=list(list(url=scripturl)))
}
