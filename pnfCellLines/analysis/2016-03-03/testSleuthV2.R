##now we want to test diffex
source("../../bin/RNASeqDiffEx.R")


##first get total data frame
full.df<-buildDiffExDf(outliers=c("pNF95.11bC"))

#build target map once
#buildGencodeTargetMap(full.df)

##now compute different models of differential expression

##first do immortalized vs. primary
fullmod<-buildSleuthModel(full.df,inc=c('Sex','OneAllele'),test='Culture',alt='primary')

##
tab<-getSleuthTable(fullmod,'Culture','primary')

##wirte list for go
write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')

pcoding<-tab[grep('protein_coding',tab[,1]),]

fground=unique(tab[which(as.numeric(tab[,'qval'])<0.1),'gene'])
bground=unique(tab[which(as.numeric(tab[,'qval'])>0.1),'gene'])
bground=union(bground,fground)
write(fground,'prim_immor_q1_diffexgenes.txt')
write(bground,'prim_immor_q1_bggenes.txt')


fground=unique(tab[which(as.numeric(tab[,'qval'])<0.01),'gene'])
bground=unique(tab[which(as.numeric(tab[,'qval'])>0.01),'gene'])
bground=union(bground,fground)
write(fground,'prim_immor_q01_diffexgenes.txt')
write(bground,'prim_immor_q01_bggenes.txt')


fground=unique(tab[which(as.numeric(tab[,'pval'])<0.01),'gene'])
bground=unique(tab[which(as.numeric(tab[,'pval'])>0.01),'gene'])
bground=union(bground,fground)
write(fground,'prim_immor_p01_diffexgenes.txt')
write(bground,'prim_immor_p01_bggenes.txt')


allgenes=unique(tab[,'gene'])##pretty sure this preserves order
write(allgenes,file='prim_immor_allgenes_ranked.txt')

##now let's get protein coding
pcgenes<-unique(tab[grep('protein_coding',tab[,'target_id']),'gene'])
write(pcgenes,file='prim_immor_all_protcod_genes_ranked.txt')


##now try to plot
fground.tf=unique(tab[which(as.numeric(tab[,'qval'])<0.01),1])

plotVals(fullmod,qval=0.1,test='Culture',alt='primary')
plotVals(fullmod,qval=0.1,ttype=c('protein_coding'),test = 'Culture',alt='primary')

####
##then do one-allele test again.
fullmod<-buildSleuthModel(full.df)

##
tab<-getSleuthTable(fullmod,test='OneAllele',alt='+')

##wirte list for go
write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')

pcoding<-tab[grep('protein_coding',tab[,1]),]
 
fground=unique(tab[which(as.numeric(tab[,'qval'])<0.1),'gene'])
bground=unique(tab[which(as.numeric(tab[,'qval'])>0.1),'gene'])
bground=union(bground,fground)
write(fground,'oneallele_q1_diffexgenes.txt')
write(bground,'oneallele_q1_bggenes.txt')


fground=unique(tab[which(as.numeric(tab[,'qval'])<0.01),'gene'])
bground=unique(tab[which(as.numeric(tab[,'qval'])>0.01),'gene'])
bground=union(bground,fground)
write(fground,'q01_diffexgenes.txt')
write(bground,'q01_bggenes.txt')


fground=unique(tab[which(as.numeric(tab[,'pval'])<0.01),'gene'])
bground=unique(tab[which(as.numeric(tab[,'pval'])>0.01),'gene'])
bground=union(bground,fground)
write(fground,'one_allele_p01_diffexgenes.txt')
write(bground,'one_allele_p01_bggenes.txt')


allgenes=unique(tab[,'gene'])##pretty sure this preserves order
write(allgenes,file='allgenes_ranked.txt')

##now let's get protein coding
pcgenes<-unique(tab[grep('protein_coding',tab[,'target_id']),'gene'])
write(pcgenes,file='all_protcod_genes_ranked.txt')


##now try to plot
fground.tf=unique(tab[which(as.numeric(tab[,'qval'])<0.01),1])

plotVals(fullmod,qval=0.1)
plotVals(fullmod,qval=0.1,ttype=c('protein_coding'))


pdf=subset(full.df,Culture=='primary')
pmod<-buildSleuthModel(pdf,inc=c("Sex"))

##
ptab<-getSleuthTable(pmod,test='OneAllele',alt='+')
plotVals(pmod,qval=0.1,prefix='primaryOnly')
plotVals(pmod,qval=0.1,ttype=c('protein_coding'),prefix='primaryOnly')

plotVals(pmod,qval=0.01,prefix='primaryOnly')
plotVals(pmod,qval=0.01,ttype=c('protein_coding'),prefix='primaryOnly')

fground=unique(tab[which(as.numeric(ptab[,'qval'])<0.01),'gene'])
bground=unique(tab[which(as.numeric(ptab[,'qval'])>0.01),'gene'])
bground=union(bground,fground)
write(fground,'primOnly_q01_diffexgenes.txt')
write(bground,'primOnly_q01_bggenes.txt')

