# Quality Control Figures for SNP data
# Xindi Guo

library(ggplot2)
library(synapseClient)
source("../dataAccess/CNVData.R")

# Get a table of files using updated annotations
snpfiles <- synapseQuery('SELECT id,sampleIdentifier,nf1Genotype FROM entity where parentId=="syn4988794"')
snpfiles <- snpfiles[which(snpfiles$entity.sampleIdentifier != "Not Applicable"),]

# Get the sample data
all.cnv <- tier0_rawData(annot=annotes.snp)

lrr <- do.call("cbind", lapply(all.cnv, function(x) x$"Log R Ratio"))
baf <- do.call("cbind", lapply(all.cnv, function(x) x$"B Allele Freq"))

list.of.lists<-c(Log.R.Ratio=c(),B.Allele.Freq=c(),sampleIdentifier=c())
for(i in names(all.cnv)){
  list.of.lists$Log.R.Ratio=c(list.of.lists$Log.R.Ratio,lrr[,i])
  list.of.lists$B.Allele.Freq=c(list.of.lists$B.Allele.Freq,baf[,i])
  list.of.lists$sampleIdentifier=c(list.of.lists$sampleIdentifier,rep(i,nrow(lrr)))
}
df<-data.frame(list.of.lists)

names(snpfiles) <- c("nf1Genotype","sampleIdentifier","synapseId")
df$nf1Genotype <- snpfiles$nf1Genotype[match(df$sampleIdentifier, snpfiles$sampleIdentifier)]

# Plots
pl=ggplot(df,aes(y=Log.R.Ratio,x=sampleIdentifier))+geom_violin(aes(colour=nf1Genotype)) + coord_flip()
pdf('rotated_violinLrrPlot.pdf')
print(pl)
dev.off()

pb=ggplot(df,aes(y=B.Allele.Freq,x=sampleIdentifier))+geom_violin(aes(colour=nf1Genotype))+coord_flip()
pdf('rotated_violinBafPlot.pdf')
print(pb)
dev.off()

#### upload files
snpqc='syn8498849'
scripturl='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/qualityControlPlots/qc_figures_snp.R'

for(file in c('rotated_violinBafPlot.pdf','rotated_violinLrrPlot.pdf')){
  synStore(File(file,parentId=snpqc),used=list(list(entity="syn4988794",wasExecuted=FALSE),
  	list(url=scripturl,wasExecuted=TRUE)))
}