##now use kallisto matrix to see if we have any correlates in drug data
source("../../bin/crossDataComps.R")

#first let's re-assess nf1 expression
tpm.mat<-rnaGencodeKallistoMatrix()
ecounts.mat<-rnaGencodeKallistoMatrix(metric='est_counts')

nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat)
nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat,dolog=T)

nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts')
nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts',dolog=T)

nf1cor=drugRna(gene='NF1',useGencode=TRUE)
krascor=drugRna(gene='KRAS',useGencode=TRUE)
egfrcor=drugRna(gene='EGFR',useGencode=TRUE)

nf1cor.t=drugRna(valname='FAUC',gene='NF1',useGencode=TRUE,qthresh=0.01)
krascor.t=drugRna(valname='FAUC',gene='KRAS',useGencode=TRUE,qthresh=0.01)
egfrcor.t=drugRna(valname='FAUC',gene='EGFR',useGencode=TRUE,qthresh=0.01)


nf1cor=drugRna(gene='NF1',useGencode=TRUE,doLog=TRUE)
krascor=drugRna(gene='KRAS',useGencode=TRUE,doLog=TRUE)
egfrcor=drugRna(gene='EGFR',useGencode=TRUE,doLog=TRUE)

nf1cor.t=drugRna(valname='FAUC',gene='NF1',useGencode=TRUE,doLog=TRUE,qthresh=0.01)
krascor.t=drugRna(valname='FAUC',gene='KRAS',useGencode=TRUE,doLog=TRUE,qthresh=0.01)
egfrcor.t=drugRna(valname='FAUC',gene='EGFR',useGencode=TRUE,doLog=TRUE,qthresh=0.01)

all=plotPCA(tpm.mat,metric='tpm',ttype=c())
protcod=plotPCA(tpm.mat,metric='tpm',ttype=c('protein_coding'))
  

##plot most variable genes
plotMostVariable(tpm.mat,metric='tpm',ttype=c(),doLog=F,top=100)
plotMostVariable(tpm.mat,metric='tpm',ttype=c(),doLog=T,top=100)
plotMostVariable(ecounts.mat,metric='est_counts',ttype=c(),doLog=F,top=100)
plotMostVariable(ecounts.mat,metric='est_counts',ttype=c(),doLog=T,top=100)
plotMostVariable(tpm.mat,metric='tpm',doLog=T,top=100,ttype='protein_coding')

nf1cor.t=drugGene(valname='FAUC',gene='NF1',useGencode=TRUE)
krascor.t=drugGene(valname='FAUC',gene='KRAS',useGencode=TRUE)
egfrcor.t=drugGene(valname='FAUC',gene='EGFR',useGencode=TRUE)


nf1cor.t=drugGene(valname='FAUC',gene='NF1',useGencode=TRUE,by.pval=F)
krascor.t=drugGene(valname='FAUC',gene='KRAS',useGencode=TRUE,by.pval=F)
egfrcor.t=drugGene(valname='FAUC',gene='EGFR',useGencode=TRUE,by.pval=F)
tubb.t=drugGene(valname='FAUC',gene='TUBB',useGencode=TRUE,by.pval=F)

