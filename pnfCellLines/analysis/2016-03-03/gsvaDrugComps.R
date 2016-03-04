
##compare drug response to gene expression pathway mappings

source("../../bin/ncatsSingleAgentScreens.R")

##lastly let's look at GSVA clustrers
gsva=as.matrix(read.table(synGet('syn5689231')@filePath))
colnames(gsva)[1]='ipNF05.5 (mixed clone)'
colnames(gsva)[11]='ipNF05.5 (single clone)'
gsva.cells=intersect(colnames(gsva),rownames(drug.by.cell))
gsva.drugcor=cor(t(gsva[,gsva.cells]),drug.by.cell[gsva.cells])