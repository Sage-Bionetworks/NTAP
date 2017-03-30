# Quality Control Figures for RNA Seq data
# Xindi Guo

library(ggplot2)
library(dplyr)
library(synapseClient)
source("../dataAccess/RNASeqData.R")

outputPCA <- function(mat,annotes,fileName){
  pca_res <- prcomp(t(mat), center=F, scale=F)
  df <- data.frame(pca_res$x[,c(1:5)])
  
  percent_variation <- pca_res$sdev^2/sum(pca_res$sdev^2) * 100
  
  df<- merge(df,annotes,by.x="row.names",by.y="synapseId")
  
  p1 <- ggplot(data=df, aes(x=PC1,y=PC2,color=sampleName)) + geom_point() + theme_bw(base_size = 14)
  p1 <- p1 + xlab(paste0("PC1",' - (', round(percent_variation[1],2), '%)' ))  + ylab(paste0("PC2",' - ( ', round(percent_variation[2],2), '%)' ))
  
  p2 <- ggplot(data=df, aes(x=PC2,y=PC3,color=sampleName)) + geom_point() + theme_bw(base_size = 14)
  p2 <- p2 + xlab(paste0("PC2",' - (', round(percent_variation[2],2), '%)' ))  + ylab(paste0("PC3",' - ( ', round(percent_variation[3],2), '%)' ))
  
  pdf(fileName)
  print(p1)
  print(p2)
  dev.off()
}

minCount <- 2
minTpm <- 0.1

#### grch38 data
mat.tpm <- rnaKallistoMatrix()
mat.count <- rnaKallistoMatrix(metric = "est_counts")

annotes.grch38 <- samp.mappings[,c("Sample Name","RNA-Seq Data", "Sample Genotype")]
annotes.grch38 <- annotes.grch38[complete.cases(annotes.grch38),]
colnames(annotes.grch38) <- c("sampleName","synapseId","genotype")

# tpm
res <- tidyr::gather(as.data.frame(mat.tpm),key=synapseId,value=TPM)
res <- merge(res,annotes.grch38,by="synapseId")
res$Gene <- rep(rownames(mat.tpm),ncol(mat.tpm))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.tpm <- subset(res, TPM >= minTpm)
p<-ggplot(res.tpm,aes(y=TPM+1,x=genotype))+geom_boxplot(aes(fill=sampleName)) + scale_y_log10()
pdf('boxPlotOfTpmMin0.1.pdf')
print(p)
dev.off()

#PCA plot
outputPCA(mat=mat.tpm,annotes=annotes.grch38,fileName = 'pcaOfTpm.pdf')


## count
res <- tidyr::gather(as.data.frame(mat.count),key=synapseId,value=count)
res <- merge(res,annotes.grch38,by="synapseId")
res$Gene <- rep(rownames(mat.count),ncol(mat.count))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.count <- subset(res, count >= minCount)
p<-ggplot(res.count,aes(y=count,x=genotype))+geom_boxplot(aes(fill=sampleName))+ scale_y_log10()
pdf('boxPlotOfCountMin2.pdf')
print(p)
dev.off()

#PCA plots
outputPCA(mat=mat.count,annotes=annotes.grch38,fileName = 'pcaOfCount.pdf')


#### gencode data
mat.tpm.gencode <- rnaGencodeKallistoMatrix()
mat.count.gencode <- rnaGencodeKallistoMatrix(metric = "est_counts")

annotes.gencode <- samp.mappings[,c("Sample Name","RNA-Seq Data (Gencode)", "Sample Genotype")]
annotes.gencode <- annotes.gencode[complete.cases(annotes.gencode),]
colnames(annotes.gencode) <- c("sampleName","synapseId","genotype")

# tpm
res <- tidyr::gather(as.data.frame(mat.tpm.gencode),key=synapseId,value=TPM)
res <- merge(res,annotes.gencode,by="synapseId")
res$Gene <- rep(rownames(mat.tpm.gencode),ncol(mat.tpm.gencode))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.tpm.gencode <- subset(res, TPM >= minTpm)
p<-ggplot(res.tpm.gencode,aes(y=TPM+1,x=genotype))+geom_boxplot(aes(fill=sampleName)) + scale_y_log10()
pdf('boxPlotOfTpmMin0.1_gencode.pdf')
print(p)
dev.off()

#PCA plot
outputPCA(mat=mat.tpm.gencode,annotes=annotes.gencode,fileName = 'pcaOfTpm_gencode.pdf')


## count
res <- tidyr::gather(as.data.frame(mat.count.gencode),key=synapseId,value=count)
res <- merge(res,annotes.gencode,by="synapseId")
res$Gene <- rep(rownames(mat.count.gencode),ncol(mat.count.gencode))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.count.gencode <- subset(res, count >= minCount)
p<-ggplot(res.count.gencode,aes(y=count,x=genotype))+geom_boxplot(aes(fill=sampleName))+ scale_y_log10()
pdf('boxPlotOfCountMin2_gencode.pdf')
print(p)
dev.off() 

#PCA plot
outputPCA(mat=mat.count.gencode,annotes=annotes.gencode,fileName = 'pcaOfCount_gencode.pdf')


#### upload files
rnaqc <- 'syn8498844'
scripturl <- 'https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/qualityControlPlots/qc_figures_rnaSeq.R'

uploadFile2Synapse <- function(file,parentId,usedEnt,usedScript){
  synStore(File(path = file,parentId = parentId),
           used=list(list(entity=usedEnt,wasExecuted=FALSE),
                     list(url=usedScript,wasExecuted=TRUE)))
}
#the count files
uploadFile2Synapse(file = "boxPlotOfCountMin2.pdf",parentId = rnaqc,usedEnt = "syn5562376", usedScript = scripturl)
uploadFile2Synapse(file = "pcaOfCount.pdf", parentId = rnaqc,usedEnt = "syn5562376", usedScript = scripturl)
uploadFile2Synapse(file = "boxPlotOfCountMin2_gencode.pdf", parentId = rnaqc, usedEnt = "syn5580378", usedScript = scripturl)
uploadFile2Synapse(file = "pcaOfCount_gencode.pdf", parentId = rnaqc,usedEnt = "syn5580378", usedScript = scripturl)

#the tpm files
uploadFile2Synapse(file = "boxPlotOfTpmMin0.1.pdf", parentId = rnaqc, usedEnt = "syn5562378", usedScript = scripturl)
uploadFile2Synapse(file = "pcaOfTpm.pdf", parentId = rnaqc, usedEnt = "syn5562378", usedScript = scripturl)
uploadFile2Synapse(file = "boxPlotOfTpmMin0.1_gencode.pdf", parentId = rnaqc, usedEnt = "syn5580347", usedScript = scripturl)
uploadFile2Synapse(file = "pcaOfTpm_gencode.pdf", parentId = rnaqc, usedEnt = "syn5580347", usedScript = scripturl)

