source("../../bin/crossDataComps.R")

drugRna('MAXR','NF1')
drugRna('MAXR','KRAS')
drugRna('MAXR','BRAF')
drugRna('MAXR','ENAH')
drugRna('MAXR','NRAS')
drugRna('TAUC','EGFR')

drugRna('TAUC','NF1')
drugRna('TAUC','KRAS')
drugRna('TAUC','BRAF')

drugRna('TAUC','NRAS')
drugRna('TAUC','ENAH')

drugGene('TAUC','KRAS')
drugGene('TAUC','NF1')
drugGene('TAUC','BRAF')
drugGene('TAUC','EGFR')

drugGene('FAUC','KRAS')
drugGene('FAUC','NF1')
drugGene('FAUC','BRAF')
drugGene('FAUC','EGFR')

