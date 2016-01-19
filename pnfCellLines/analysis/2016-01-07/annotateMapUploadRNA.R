###get all RNA-Seq files from belltown directory where kallisto was run
##map ENS transcripts to gene id
##create new tab-delimited file with transcript AND gene id
##annotate with appropriate cell line and upload



require(synapseClient)
library(data.table)
synapseLogin()
belltown.dir='../2016-01-05/'
filedirs=list.files(belltown.dir)
fqc=c(grep('fastqc',filedirs),grep('.out',filedirs))
if(length(fqc)>0)
    filedirs=filedirs[-fqc]
                                        #each directory is a different lane/run.
filedirs=filedirs[grep("C85",filedirs)]
all.dat<-lapply(filedirs,function(x) as.data.frame(fread(paste(belltown.dir,x,'abundance.tsv',sep='/'),sep='\t')))
names(all.dat)<-filedirs


##first get all gene names, add to files
all.genes=all.dat[[1]][,1]
newfiles<-lapply(all.dat,function(x){
    dat<-x[,-1]
    ids<-t(sapply(x[,1],function(x) unlist(strsplit(x,split='|',fixed=T))))
    colnames(ids)<-c('EnsGene','EnsTrans','OttGene','OttTrans','HugoTrans','HugoSymbol','Length','Descrip')
    new.dat<-data.frame(ids,dat)
    return(new.dat)
})

#require(biomaRt)
#ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host='www.ensembl.org')
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)

#epep="ensembl_transcript_id"
#egene='hgnc_symbol'
#gene.mapping<-getBM(attributes=c(epep,egene),filters=c(epep),values=as.list(all.genes),mart=ensembl)

#now create/write new files
#newfiles<-lapply(all.dat,function(x){
#  HGNCSymbol=gene.mapping[match(x[,1],gene.mapping[,1]),2]
#  newf=data.frame(HGNCSymbol,x)
 # return(newf)
#})
#names(newfiles)=names(all.dat)

##then create annotations
samp.mapping=read.table(synGet('syn5562006')@filePath,sep=',',header=T)

tab.q=synTableQuery('SELECT distinct "Sample Name","Sample ID","Sample Origin","Sample Genotype" FROM syn5014742')@values

samp.to.cell<-lapply(filedirs,function(x){
  cellLine=samp.mapping$CellLine[match(paste(x,'0',sep='_'),samp.mapping$File_Name)]
  res=c(tab.q[match(cellLine,tab.q$`Sample Name`),],FlowCellLane=x)
  names(res)<-c('sampleName','sampleID','sampleOrigin','sampleGenotype','flowCellLane')
  res
})
names(samp.to.cell)=filedirs

##write files
file.res<-lapply(filedirs,function(x){
  filename=paste(x,'_RNASeq_Kallisto_gencodev24_quants.tsv',sep='')
  write.table(newfiles[[x]],file=filename,sep='\t',row.names=F)
  newf=File(path=filename,parentId='syn5579785')
  synSetAnnotations(newf)<-samp.to.cell[[x]]
  newf=synStore(newf,executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-01-07/annotateMapUploadRNA.R')
})
##store on synapse
