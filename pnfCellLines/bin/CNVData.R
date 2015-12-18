###
### CNVData
### Primary R script designed to make CNV data from PNF cell lines easily available.
###

library(data.table)
library(DNAcopy)
library(CNTools)

library(synapseClient)
synapseLogin()

##download all meta data from original files
if(!exists('snpfiles'))
    snpfiles<-synapseQuery('SELECT id,sampleName,sampleGenotype,sampleID,sampleOrigin FROM entity where parentId=="syn4988794"')

snpfiles=snpfiles[which(!is.na(snpfiles$entity.sampleID)),]

#now get annotations
origin<-snpfiles$'entity.sampleOrigin'
names(origin)<-snpfiles$entity.sampleID

clnames<-paste(snpfiles$'entity.sampleName',snpfiles$'entity.sampleGenotype')
names(clnames)<-snpfiles$entity.sampleID

genotype<-snpfiles$'entity.sampleGenotype'
names(genotype)<-clnames

cnv_annotation_data<-function(){
    annots=snpfiles
    colnames(annots)<-c('Origin','sample','Genotype','synapseID','Name')
    return(annots)
}

##here is a basic function to get annotation data
snp_annotation_data<-function(){
    ##need to downlod and read in large annotation file as well
    print("Retrieving OMNI Array SNP annotation data from Synapse...")
    anndata<-synGet('syn5005069')
    annot <- as.data.frame(fread(anndata@filePath,sep=",",header=T))
    return(annot)
}

#first function is to get 'raw' data from SNP OMNI arrays under syn4988794
#and return as single data frame
tier0_rawData<-function(annot=NA){
                                        #first login to synapse
    synapseLogin()

    if(is.na(annot))
        annot=snp_annotation_data()

    print('Now retreiving original CNV data from NTAP CLP OMNI arrays...')
                                        #collect sample files and process them into data frame
    sample.data<-lapply(snpfiles$entity.id,function(synid){
        fname=synGet(synid)
        print(paste("Getting sample",snpfiles$entity.sampleID[match(synid,snpfiles$entity.id)]))
        data <- as.data.frame(fread(fname@filePath,sep=",",header=T))
        ad<-data[match(annot$Name,data$'SNP Name'),]
        return(ad)
    })

                                        #update sample names
    names(sample.data) <- snpfiles$entity.sampleID

    return(sample.data)
}


#this function collects the segmented data, which has been processed and analyzed and uploaded to synapse
tier1_segmentedData<-function(filterBySD=TRUE){
    print(paste('Retriving segmented CNV data from Synapse...'))
    if(filterBySD){ #get the segments that differ by at least 2 SDs with neighboring regions
        segdata <- read.table(synGet('syn5015036')@filePath,header=T)
    }else{
        segdata <- read.table(synGet('syn5015035')@filePath,header=T)
    }
    return(segdata)
}




dataSummaryPlots<-function(annot=NA){

    if(is.na(annot))
        annot<-snp_annotation_data


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
    pdf('ntap_cnv_plots.pdf')

    m<-ggplot(df,aes(x=LogRRatio,colour=SampleType,linetype=NF1Region))
    m<-m + geom_density() + xlim(-2.5,2)
    print(m)

    m<-ggplot(df,aes(x=BAlleleFreq,colour=SampleType,linetype=NF1Region))
    m<-m + geom_density() + xlim(-.1,1.1)
    print(m)

    dev.off()

    #sf=File('ntap_cnv_plots.pdf',parentId='syn5014748')
    #synStore(sf)

    ##what if we just plot LRR and BAF values within region of interest?
    ##now get the region around chr17 to get region of interest, then re-plot b-allele frequency and logR
    nf1.df<-subset(df,Chromosome=='17')
    nf1.df<-subset(nf1.df,Position>28000000&Position<31000000)
    pdf('ntap_cnv_chr17_values.pdf')
    m<-ggplot(nf1.df)
    m<-m +geom_point(aes(x=Position,y=BAlleleFreq,colour=Origin,shape=SampleType))
    print(m)

}
