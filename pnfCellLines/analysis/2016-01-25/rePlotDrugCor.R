##now use kallisto matrix to see if we have any correlates in drug data
source("../../bin/crossDataComps.R")

#first let's re-assess nf1 expression
tpm.mat<-rnaGencodeKallistoMatrix()
ecounts.mat<-rnaGencodeKallistoMatrix(metric='est_counts')



#nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat)
#nf1.expr<-plotTranscriptsOfGene("NF1",tpm.mat,dolog=T)

#nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts')
#nf1.expr<-plotTranscriptsOfGene("NF1",ecounts.mat,metric='est_counts',dolog=T)

###collapse by protein-coding
pcvals<-grep('protein_coding',rownames(tpm.mat))
prot_coding.tpm=tpm.mat[pcvals,]

#expression of transcripts
#nf1.expr<-plotTranscriptsOfGene("NF1",prot_coding.tpm,metric='protCodingTpm')
#nf1.expr<-plotTranscriptsOfGene("NF1",prot_coding.tpm,metric='protCodingTpm',dolog=T)

##then also collapse by gene
#all.genes<-unique(sapply(rownames(tpm.mat),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))
#pc.genes<-unique(sapply(rownames(prot_coding.tpm),function(x) unlist(strsplit(x,split='.ENST',fixed=T))[1]))

##add in option to collapse by gene

for(co in c(TRUE,FALSE)){
  for(pc in c(TRUE,FALSE)){
    for(dl in c(TRUE,FALSE)){
    nf1cor=drugRna(gene='NF1',useGencode=TRUE,doLog=dl,collapseAllCounts=co,proteinCoding=pc)
    krascor=drugRna(gene='HRAS',useGencode=TRUE,doLog=dl,collapseAllCounts=co,proteinCoding=pc)
    egfrcor=drugRna(gene='PLXDC2',useGencode=TRUE,doLog=dl,collapseAllCounts=co,proteinCoding=pc)
    ablcor=drugRna(gene='ABL1',useGencode=TRUE,doLog=dl,collapseAllCounts=co,proteinCoding=pc)
    pdgcor=drugRna(gene='PDGFRA',useGencode=TRUE,doLog=dl,collapseAllCounts=co,proteinCoding=pc)
    
        }
  }
}




