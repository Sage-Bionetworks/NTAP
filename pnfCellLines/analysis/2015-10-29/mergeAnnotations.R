##This is a short script designed to merge the annotation files provided from UCSC (hg19) to be a 'geneInfo' file for CNTools

xref<-read.table('../../data/ucsc_kgXref_hg19_2015_10_29.csv',header=F,as.is=T,sep='\t',quote='"')[,c(1,5)]


genepos<-read.table('../../data/ucsc_knownGene_hg19_2015_10_29.csv',header=F,as.is=T,sep='\t',quote='"')


geneInfo<-cbind(genepos[,c(2,4,5,1)],xref[match(genepos[,1],xref[,1]),2])
geneInfo[,1]=sapply(geneInfo[,1],function(x) gsub("chr",'',x))
colnames(geneInfo)<-c('chrom','start','end','geneid','genename')
write.table(geneInfo,file='../../data/hg19_geneInfo.txt')
