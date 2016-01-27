##rnaseq boxplots

source("../../bin/RNASeqData.R")

tpm.mat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
library(dplyr)

alldat=lapply(colnames(tpm.mat),function(x) data.frame(CellLine=rep(x,nrow(tpm.mat)),TPM=tpm.mat[,x]))
df<-data.frame(CellLine=unlist(lapply(alldat,function(x) x$CellLine)),TPM=unlist(lapply(alldat,function(x) x$TPM)))

require(ggplot2)
p=ggplot(df)+geom_boxplot(aes(x=CellLine,y=TPM))+scale_y_log10()+ggtitle("Transcripts per million across cell lines")+theme(axis.text.x=element_text(angle = -90, hjust = 0))
png('tpmValsInBoxplot.png')
print(p)
dev.off()


###now do PCA without outlier

source("../../bin/RNASeqDiffEx.R")
df<-buildDiffExDf()
red.df<-df[-which(df$sample=='pNF95.11bC'),]

genotype.mod<-buildSleuthModel(red.df,inc=c('Culture'),test='Genotype',alt='--')
res=makeTablesAndPlots(genotype.mod,test='Genotype',alt='--',prefix='outlierRemoved')

primary.only<-subset(red.df,Culture=='primary')
prim.mod<-buildSleuthModel(red.df,inc=c(),test='OneAllele',alt='+')
prim.res=makeTablesAndPlots(prim.mod,test='OneAllele',alt='+',prefix='outlierRemovedPrimaryOnly')




