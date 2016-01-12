##now we want to test diffex
library('sleuth')
library('synapseClient')

h5files<-synapseQuery("select name,id from entity where parentId=='syn5579785'")
h5files=h5files[grep('*h5',h5files[,1]),]

h5s=lapply(h5files[,2],function(x) synGet(x))
snames=lapply(h5s,function(x) x@annotations$sampleName)
sgens=lapply(h5s,function(x) x@annotations$sampleGenotype)
filepaths=lapply(h5s,function(x) x@filePath)

##build the contrast matrix from the annotations
df=data.frame(sample=unlist(snames),path=unlist(filepaths),Genotype=unlist(sgens))
df$path=as.character(df$path)
so<-sleuth_prep(df,~ Genotype)
so<-sleuth_fit(so)
so <- sleuth_wt(so, 'condition++')
