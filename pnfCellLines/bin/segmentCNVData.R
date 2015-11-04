###
### Goal of this script is to read in CNV analysis from Synapse and perform the following analyses
### 1- analyze probes in NF1 region
### 2- Do segmentation analysis, upload results and cluster
###

library(data.table)
library(DNAcopy)
library(CNTools)

library(synapseClient)

#first login to synapse
source("../../bin/CNVData.R")
                                        #read in the annotation file. this can be big.
annot<-snp_annotation_data()

#now get raw files
sample.data<-tier0_rawData(annot)

####STEP 1: analyze probes from OMNI SNP data

##get regions that are within our general region: chr17:29000019 to 30427403
chr17.snps=annot[grep('chr17',annot$chrpos),]
all.pos<-annot$pos
all.chr<-annot$chr

pos<-all.pos[which(all.chr=='17')]
chr17.snps=chr17.snps[intersect(which(pos>29000019),which(pos<30427403)),]


is.autosome <- as.character(all.chr) %in% as.character(1:22)
lrr <- do.call("cbind", lapply(sample.data, function(x) x$"Log R Ratio"))
baf <- do.call("cbind", lapply(sample.data, function(x) x$"B Allele Freq"))

##add variable to signify if SNP is in region of interest
in.region<-rep(FALSE,nrow(lrr))
in.region[which(annot$chrpos%in%chr17.snps$chrpos)]<-TRUE
rm(annot)

#how can we plot baf only?
library(ggplot2)
##create ggplot-associated data.frame
list.of.lists<-c(LogRRatio=c(),BAlleleFreq=c(),Sample=c(),SampleType=c(),NF1Region=c(),Position=c(),Chromosome=c(),Origin=c())

for(i in names(sample.data)){
    #first get general population distribution
    list.of.lists$LogRRatio=c(list.of.lists$LogRRatio,lrr[,i])
    list.of.lists$BAlleleFreq=c(list.of.lists$BAlleleFreq,baf[,i])
    list.of.lists$Sample=c(list.of.lists$Sample,rep(i,nrow(lrr)))
    if(i=='CIDR'){
        st<-rep('Control',nrow(lrr))

    }else if (i%in%c('W1','W2','W3','W4','W5')){
        st<-rep('Primary',nrow(lrr))

    }else{
        st<-rep('CellLine',nrow(lrr))

    }
    list.of.lists$SampleType=c(list.of.lists$SampleType,st)
    list.of.lists$NF1Region=c(list.of.lists$NF1Region,in.region)
    list.of.lists$Position=c(list.of.lists$Position,all.pos)
    list.of.lists$Chromosome=c(list.of.lists$Chromosome,all.chr)
    list.of.lists$Origin=c(list.of.lists$Origin,rep(origin[which(names(sample.data)==i)],nrow(lrr)))

}
df<-data.frame(list.of.lists)
rm(list.of.lists)

#now do the plotting
pdf('ntap_cnv_dist_plots.pdf')

m<-ggplot(df,aes(x=LogRRatio,colour=SampleType,linetype=NF1Region))
m<-m + geom_density() + xlim(-2.5,2)
print(m)

m<-ggplot(df,aes(x=BAlleleFreq,colour=SampleType,linetype=NF1Region))
m<-m + geom_density() + xlim(-.1,1.1)
print(m)

dev.off()

sf=File('ntap_cnv_dist_plots.pdf',parentId='syn5014748')
synStore(sf)

##what if we just plot LRR and BAF values within region of interest?
##now get the region around chr17 to get region of interest, then re-plot b-allele frequency and logR
nf1.df<-subset(df,Chromosome=='17')
nf1.df<-subset(nf1.df,Position>28000000&Position<31000000)
pdf('ntap_cnv_chr17_values.pdf')
m<-ggplot(nf1.df)
m<-m +geom_point(aes(x=Position,y=BAlleleFreq,colour=Origin,shape=SampleType))
print(m)

m<-ggplot(nf1.df,aes(x=Position,y=LogRRatio,colour=Origin,shape=SampleType))
m<-m +geom_point()
print(m)

##now do individual subset
for(i in unique(nf1.df$Origin)){
    tdf<-subset(nf1.df,Origin==i)
    m<-ggplot(tdf)
    m<-m +geom_point(aes(x=Position,y=BAlleleFreq,colour=Sample,shape=SampleType))+ggtitle(paste('NF1 region for sample',i))
    print(m)

    m<-ggplot(tdf,aes(x=Position,y=LogRRatio,colour=Sample,shape=SampleType))
    m<-m +geom_point()+ggtitle(paste('NF1 region for sample',i))
    print(m)
}

dev.off()
sf=File('ntap_cnv_chr17_values.pdf',parentId='syn5014748')
synStore(sf)
###STEP 3: do segmentation analysis

cna <- CNA(lrr[is.autosome,], as.character(all.chr)[is.autosome], all.pos[is.autosome], data.type="logratio",names(sample.data))


smoothed.cna <- smooth.CNA(cna)
segment.smoothed.cna <- segment(smoothed.cna, verbose=1)


segment.smoothed.cna.sundo <- segment(smoothed.cna, undo.splits="sdundo",undo.SD=2,verbose=1)
chr17.smoothed.sundo=subset(segment.smoothed.cna.sundo,chromlist=c("17"))


##write the files
write.table(segment.smoothed.cna$output, file="ntap_clp_cbs_noundo.seg",
            sep="\t",quote=FALSE,row.names=F)
write.table(segment.smoothed.cna.sundo$output, file="ntap_clp_cbs_undosd2.seg",
            sep="\t",quote=FALSE,row.names=F)

##put these files in synapse analysis directory
sf=File('ntap_clp_cbs_noundo.seg',parentId='syn5014748')
synStore(sf,used=list(list(name='segmentCNVData.R',
                url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/segmentCNVData.R',wasExecuted=TRUE),
                list(entity='syn5005069',wasExecuted=FALSE),
                list(name='CNVData.R',
                     url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/CNVData.R',wasExecuted=TRUE)),
         activityName='Segmentation analysis of copy alterations')



sf=File('ntap_clp_cbs_undosd2.seg',parentId='syn5014748')
synStore(sf,used=list(list(name='segmentCNVData.R',
                url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/segmentCNVData.R',wasExecuted=TRUE),
                list(entity='syn5005069',wasExecuted=FALSE),
                list(name='CNVData.R',
                     url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/CNVData.R',wasExecuted=TRUE),
                list(name='CNVData.R',
                     url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/CNVData.R',wasExecuted=TRUE)),
         activityName='Segmentation analysis of copy alterations with filtering')
