##exome seq data

require(synapseClient)
synapseLogin()

getAnnovarFiles<-function(reStore=FALSE){
  if(reStore){
  all.files<-synQuery('select name,id from entity where parentId=="syn6086887"')
  all.files$SeqId<-sapply(all.files$entity.name,function(x) unlist(strsplit(x,split='_'))[1])
  ann.files<-all.files[grep('txt',all.files$entity.name),]
  cell.line.map<-synTableQuery('SELECT "Sample Name","Exome-Seq Identifier","Sample Genotype" FROM syn5014742 where "Exome-Seq Identifier" is not NULL')@values
  
  ann.files$CellLine=cell.line.map$`Sample Name`[match(ann.files$SeqId,cell.line.map$`Exome-Seq Identifier`)]
  ann.files$Genotype=cell.line.map$`Sample Genotype`[match(ann.files$SeqId,cell.line.map$`Exome-Seq Identifier`)]
  
  all.dat<-do.call('rbind',lapply(ann.files$entity.id,function(x) {
    tab<-read.table(synGet(x)@filePath,sep='\t',header=T)
    tab$CellLine<-rep(ann.files$CellLine[match(x,ann.files$entity.id)],nrow(tab))
    tab$NF1Genotype<-rep(ann.files$Genotype[match(x,ann.files$entity.id)],nrow(tab))
    
    return(tab)}))
  write.table(all.dat,file='annovarSamplesWithCellLineNameGenotype.tsv',sep='\t',row.names=F)
  synStore(File('annovarSamplesWithCellLineNameGenotype.tsv',parentId='syn6086887'),executed=list(list(url="https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/exomeSeqProcessing.R")))
    }
  else{
    all.dat<-read.table(synGet('syn6174638')@filePath,sep='\t',header=T)
  }
  return(all.dat)
  
  }