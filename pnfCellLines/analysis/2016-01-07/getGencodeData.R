source('../../bin/RNASeqData.R')

tpm.mat<-rnaGencodeKallistoMatrix(buildFromFiles=TRUE)
ecounts.mat<-rnaGencodeKallistoMatrix(buildFromFiles=TRUE,metric='est_counts')

##now get some plots going!

plotTranscriptsOfGene("NF1",tpm.mat,'tpm')
plotTranscriptsOfGene("NF1",ecounts.mat,'eCounts')
plotTranscriptsOfGene("NF1",tpm.mat,'tpm',dolog=TRUE)
plotTranscriptsOfGene("NF1",ecounts.mat,'eCounts',dolog=TRUE)

plotTranscriptsOfGene("NF1",tpm.mat,'tpm',dolog=TRUE,ttype=c("protein_coding"))
plotTranscriptsOfGene("NF1",ecounts.mat,'eCounts',dolog=TRUE,ttype=c('protein_coding'))

plotPCA(ecounts.mat,metric='est_counts')
plotPCA(tpm.mat)

plotPCA(ecounts.mat,metric='est_counts',ttype=c('protein_coding'))
plotPCA(tpm.mat,ttype=c('protein_coding'))