##apply gene expression profiles to GSEA/MSIGDB using GSVA
library(GSVA)

source('../../bin/RNASeqData.R')


library(data.table)
library(GSEABase)
library(GSVAdata)
library(pheatmap)
#cut and paste from vignette
#cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,
#                            +                            min.sz=10, max.sz=500, verbose=TRUE)$es.obs,
#        +                            dir=cacheDir, prefix=cachePrefix)

data(c2BroadSets)

#> gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)$es.obs
rna.data<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)

#get entrez ids for RNA count matrix
eids<-as.data.frame(fread('../../data/HugoGIDsToEntrez_DAVID.txt',sep='\t',header=T))

eid.match=eids[match(rownames(rna.data),eids[,1]),2]
eid.rows=which(!is.na(eid.match))
ent.rna<-rna.data[eid.rows,]
rownames(ent.rna)=eid.match[eid.rows]

es<-gsva(ent.rna,c2BroadSets,rnaseq=T,no.bootstraps=1000)
##let's write out the pathways to a tab-delimited file for future analysis


##pathway analysis directory
synid='syn5688828'
path.vals=es$es.obs
path.bootstrap=es$bootstrap$p.vals.sign
write.table(path.vals,'cellLineBroadPathEnrich.tab',sep='\t')
write.table(path.bootstrap,'cellLineBroadPath1000BootstrapPvals.tab',sep='\t')



#let's cluster, see how that goes? 
res=path.vals
plot(hclust(as.dist(1-cor(res)),method='ward.D2'))
colnames(res)[1]="ipNF05.5 (mixed clone)" 
colnames(res)[11]="ipNF05.5 (single clone)"

genotype=samp.mappings[match(colnames(res),samp.mappings[,1]),4]
names(genotype)<-colnames(res)

vars<-apply(res,1,var,na.rm=T)
mostvar<-res[order(vars,decreasing=T)[1:100],]
#names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10, annotation_col=data.frame(Genotype=genotype),
         filename='mostVariablePathways_GSVA.png')


mostvar<-res[order(vars,decreasing=F)[1:100],]
#names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10, 
         annotation_col=data.frame(Genotype=genotype),
         filename='leastVariablePathways_GSVA.png')


require(limma)

design= model.matrix(~factor(genotype))
colnames(design)=c('','+-','++')
fit <- lmFit(res, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="++")
sig=res[rownames(allGenes)[which(allGenes$P.Value<0.05)],]
pheatmap(sig,cellheight=10,cellwidth=10,annotation_col=data.frame(Genotype=genotype),
         filename='sigDiff_HomoZ_Pathways_GSVA.png')



for(file in c('cellLineBroadPathEnrich.tab','sigDiff_HomoZ_Pathways_GSVA.png','cellLineBroadPath1000BootstrapPvals.tab','mostVariablePathways_GSVA.png','leastVariablePathways_GSVA.png')){
  sf=File(file,parentId=synid)
  synStore(sf,activityName='GSVA enrichment analysis',
           used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-25/evalRnaUsingGSVA.R',wasExecuted=TRUE)))
}
