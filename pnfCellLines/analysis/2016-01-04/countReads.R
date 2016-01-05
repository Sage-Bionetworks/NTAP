###new attempt to count reads
alldirs=list.dirs('../2015-12-22')
library('jsonlite')
for(a in alldirs[-1]){
    abund=paste(a,'abundance.tsv',sep='/')
    json=paste(a,'run_info.json',sep='/')
    tab<-read.table(abund,sep='\t',header=T)
    allcounts=sum(tab[,4])
    tot=jsonlite::fromJSON(json)$n_processed
    print(paste(allcounts,tot,allcounts/tot,a))
}
