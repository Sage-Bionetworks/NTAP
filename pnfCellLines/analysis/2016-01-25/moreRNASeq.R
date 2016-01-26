##additional RNAseq analysis

source("../../bin/RNASeqDiffEx.R")

#first read in old gene lists for functional enrichment
culture.genes<-read.table('../2016-01-20/Culture_variable_primary_testSleuthResults.csv',sep=',')
genotype.genes<-read.table('../2016-01-20/Genotype_variable_--_testSleuthResults.csv',sep=',')
oa.genes<-read.table("../2016-01-20/OneAllele_variable_+_testSleuthResults.csv",sep=',')

write(unique(as.character(culture.genes$gene)),file='rankedCultureExpression.txt')
write(unique(as.character(genotype.genes$gene)),file='rankedGenotypeExpression.txt')
write(unique(as.character(oa.genes$gene)),file='rankedOneAlleleExpression.txt')

write(unique(as.character(culture.genes$gene[which(culture.genes$qval<0.1)])),file='sigDiffExGenesCulture.txt')
write(unique(as.character(genotype.genes$gene[which(genotype.genes$qval<0.1)])),file='sigDiffExGenesGenotype.txt')
write(unique(as.character(oa.genes$gene[which(oa.genes$qval<0.1)])),file='sigDiffExGenesOneAllele.txt')

##now do primary comparison
df<-buildDiffExDf()

prim<-subset(df,Culture=='primary')

genotype.mod<-buildSleuthModel(df,inc=c('Culture'),test='Genotype',alt='--')
culture.mod<-buildSleuthModel(df,inc=c('Genotype'),test='Culture',alt='primary')
oneallele.mod<-buildSleuthModel(df,inc=c('Culture'),test='OneAllele',alt='+')
primonly.mod<-buildSleuthModel(prim,inc=c(),test='Genotype',alt='--')
res=makeTablesAndPlots(primonly.mod,test='Genotype',alt='--',prefix='primaryOnly')
#now do sleuth live and annotate volcano plots. 

primonly.genes<-read.table('Genotype_variable_--_test_primaryOnlySleuthResults.csv',sep=',')
write(unique(as.character(primonly.genes$gene[which(primonly.genes$qval<0.1)])),file='sigDiffExGenesPrimaryOnly.txt')
