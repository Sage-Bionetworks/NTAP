##now we want to test diffex

source("../../bin/RNASeqDiffEx.R")

fulldf<-buildDiffExDf()


doAnalysis<-function(model, test,alt){
    fname=paste(test,'variable',alt,'test',sep='_')
    ##1
    tab<-getSleuthTable(model,test,alt)
                                        #store table
    write.table(tab,file=paste(fname,'SleuthResults.csv',sep=''),sep=',',quote=F)
    sigs<-which(as.numeric(tab[,'qval'])<0.1)

    pc=sigs[grep('protein_coding',tab[sigs,'target_id'])]

    pcg=unique(tab[pc,'gene'])
    print(paste("Found",length(sigs),'differentially expressed transcripts, of which',length(pc),'are protein coding representing',length(pcg),'unique genes'))

    pdf(paste(fname,'VolcanoPlot.pdf',sep=''))
    p<-plot_volcano(model,paste(test,alt,sep=''),point_alpha=0.8)
    print(p)
    dev.off()

    pdf(paste(fname,'PCAPlot.pdf',sep=''))
    p<-plot_pca(model,text_labels=TRUE,color_by=test,point_size=8)
    print(p)
    dev.off()

    gvals=tab[,'gene']
    names(gvals)<-tab[,'target_id']
    plotGenesInSamples(model,tab[pc,'target_id'],units="tpm",
                       genes=NULL,annotes=NULL,collapseByGene=TRUE,
                       fname=paste(fname,'sigGenesByTpmInHeatmap.pdf'),test=test,alt=alt)


}



#genotype.mod<-buildSleuthModel(fulldf,inc=c('Culture'),test='Genotype',alt='--')
#culture.mod<-buildSleuthModel(fulldf,inc=c('Genotype'),test='Culture',alt='primary')
#oneallele.mod<-buildSleuthModel(fulldf,inc=c('Culture'),test='OneAllele',alt='+')


##now for each model, we need to
#1-create ranked list of genes
#2-do Volcano plotso
#3-do PCA
#4-do Heatmap of diff ex genes
##

doAnalysis(genotype.mod,'Genotype','--')
doAnalysis(culture.mod,'Culture','primary')
doAnalysis(oneallele.mod,'OneAllele','+')
plotSingleGene(genotype.mod,'NF1','tpmValuesFromGenotypeModel.pdf')
