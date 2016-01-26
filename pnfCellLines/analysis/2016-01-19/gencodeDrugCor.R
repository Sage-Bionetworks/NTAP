##now use kallisto matrix to see if we have any correlates in drug data
source("../../bin/crossDataComps.R")

#first let's re-assess nf1 expression
tpm.mat<-rnaGencodeKallistoMatrix()
ecounts.mat<-rnaGencodeKallistoMatrix(metric='est_counts')



nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat)
nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat,dolog=T)

nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts')
nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts',dolog=T)

###collapse by protein-coding
pcvals<-grep('protein_coding',rownames(tpm.mat))
prot_coding.tpm=tpm.mat[pcvals,]

#expression of transcripts
nf1.expr<-plotTranscriptsOfGene("NF1",prot_coding.tpm,metric='protCodingTpm')
nf1.expr<-plotTranscriptsOfGene("NF1",prot_coding.tpm,metric='protCodingTpm',dolog=T)

##then also collapse by gene
#all.genes<-unique(sapply(rownames(tpm.mat),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))
#pc.genes<-unique(sapply(rownames(prot_coding.tpm),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))

##add in option to collapse by gene

for(co in c(TRUE,FALSE)){
  for(pc in c(TRUE,FALSE)){
  nf1cor=drugRna(gene='NF1',useGencode=TRUE,collapseAllCounts=co,proteinCoding=pc)
  krascor=drugRna(gene='KRAS',useGencode=TRUE,collapseAllCounts=co,proteinCoding=pc)
  egfrcor=drugRna(gene='EGFR',useGencode=TRUE,collapseAllCounts=co,proteinCoding=pc)
  }
}


nf1cor.t=drugRna(valname='TAUC',gene='NF1',useGencode=TRUE)
krascor.t=drugRna(valname='TAUC',gene='KRAS',useGencode=TRUE)
egfrcor.t=drugRna(valname='TAUC',gene='EGFR',useGencode=TRUE)


nf1cor=drugRna(gene='NF1',useGencode=TRUE,doLog=TRUE)
krascor=drugRna(gene='KRAS',useGencode=TRUE,doLog=TRUE)
egfrcor=drugRna(gene='EGFR',useGencode=TRUE,doLog=TRUE)

nf1cor.t=drugRna(valname='TAUC',gene='NF1',useGencode=TRUE,doLog=TRUE)
krascor.t=drugRna(valname='TAUC',gene='KRAS',useGencode=TRUE,doLog=TRUE)
egfrcor.t=drugRna(valname='TAUC',gene='EGFR',useGencode=TRUE,doLog=TRUE)

all=plotPCA(tpm.mat,metric='tpm',ttype=c())
protcod=plotPCA(tpm.mat,metric='tpm',ttype=c('protein_coding'))
  

##plot most variable genes
plotMostVariable(tpm.mat,metric='tpm',ttype=c(),doLog=F,top=100)
plotMostVariable(tpm.mat,metric='tpm',ttype=c(),doLog=T,top=100)
plotMostVariable(ecounts.mat,metric='est_counts',ttype=c(),doLog=F,top=100)
plotMostVariable(ecounts.mat,metric='est_counts',ttype=c(),doLog=T,top=100)
plotMostVariable(tpm.mat,metric='tpm',doLog=T,top=100,ttype='protein_coding')


