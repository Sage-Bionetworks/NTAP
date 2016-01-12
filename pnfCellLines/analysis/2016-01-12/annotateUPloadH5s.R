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
all.dat<-lapply(filedirs,function(x) paste(belltown.dir,x,'abundance.h5',sep='/'))
names(all.dat)<-filedirs


##first get all gene names, add to files


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
  filename=paste(x,'_RNASeq_Kallisto_gencodev24_quants.h5',sep='')
  file.copy(all.dat[[x]],filename) 
 #write.table(newfiles[[x]],file=filename,sep='\t',row.names=F)
  newf=File(path=filename,parentId='syn5579785')
  synSetAnnotations(newf)<-samp.to.cell[[x]]
  newf=synStore(newf,executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2015-01-12/annotateUPloadH5s.R')
})
##store on synapse
