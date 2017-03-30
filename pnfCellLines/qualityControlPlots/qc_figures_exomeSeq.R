library(data.table)
tab<-data.frame(fread('exomeSeq/allSamples.tab',sep='\t',header=T,skip=98))
colnames(tab)<-c("CN","Discordance","NumberOfSites","AvMinDepth","Sample_i","Sample_j")

##now get annotations
require(synapseClient)
synapseLogin()

samp.mappings<-synTableQuery('SELECT "Sample Name", "Sample Genotype" FROM syn5014742')@values
colnames(samp.mappings) <- c("sampleName", "genotype")
tab$Genotype_i = samp.mappings$genotype[match(tab$Sample_i,samp.mappings$sampleName)]
tab$Genotype_j = samp.mappings$genotype[match(tab$Sample_j,samp.mappings$sampleName)]

require(reshape2)
require(pheatmap)

atab<-tab
atab$Sample_i=tab$Sample_j
atab$Sample_j=tab$Sample_i
atab$Genotype_i=tab$Genotype_i
atab$Genotype_j=tab$Genotype_j
newtab=rbind(atab,tab)

dmat<-acast(newtab,Sample_j~Sample_i,value.var="Discordance")

itab<-unique(tab[,c(5,7)])
ivals=data.frame(Genotype=itab$Genotype_i)
rownames(ivals)<-itab$Sample_i

jtab<-unique(tab[,c(6,8)])
jvals=data.frame(Genotype=jtab$Genotype_j)
rownames(jvals)<-jtab$Sample_j

pheatmap(dmat,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampleDiscordance.png')
pheatmap(dmat,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampleDiscordance.pdf')

numsites=acast(newtab,Sample_j~Sample_i,value.var='NumberOfSites')
pheatmap(numsites,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampNumOfsites.png')
pheatmap(numsites,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampNumOfsites.pdf')


for (file in c('sampleDiscordance','sampNumOfsites')){
  for(suff in c('.pdf','.png')){  
    fname=paste(file,suff,sep='')
    synStore(File(fname,parentId='syn8498855'),
             used=list(list(entity='syn6674843',wasExecuted=FALSE),
                       list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/qualityControlPlots/qc_figures_exomeSeq.R',wasExecuted=TRUE)))
  }
}
