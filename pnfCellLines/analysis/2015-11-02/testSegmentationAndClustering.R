#!/bin/bash

#test segmentation analysis
source("../../bin/segmentCNVData.R") #(better to run on EC2)

##also let's add the annotation here

#test clustering
source('../../bin/clusterCNVBySample.R')
main()
#now upload images for wiki
parentid='syn5014748'

#here are the files used for the analysis performed
filelist=list(list(name='clusterCNVBySample.R',url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/clusterCNVBySample.R',wasExecuted=TRUE),list(entity='syn5015036',wasExecuted=FALSE),list(name='CNVData.R',url='https://raw.githubusercontent.com/sgosline/NTAP/master/pnfCellLines/bin/CNVData.R',wasExecuted=TRUE))

#image 1 - unsupervised clustering
sf=File('all_gene_by_median_logRRatios_dendro.png',parentId=parentid)
#add annotations
#add in activity
synStore(sf,used=filelist,activityName='Plotting of all segmented CNV data',
         activityDescription='To execute source clusterCNVBySample.R and then run main()')

#image 2 - demonstration of gene variability
sf=File('most_variable_100_gene_by_median_logRRatios_heatmap.png',parentId=parentid)

#add annotations
#add in activity

synStore(sf,used=filelist,activityName='Plotting of genes with most variable CNAs',
         activityDescription='To execute source clusterCNVBySample.R and then run main()')


#image 3 - correlated by genotype
sf=File('most_gt_correlated_100_gene_by_median_logRRatios_heatmap.png',parentId=parentid)
#add annotations
#add in activity
synStore(sf,used=filelist,activityName=
         'Plotting of genes with CNAs that correlate with genotype',
         activityDescription
         ='To execute source clusterCNVBySample.R and then run main()')
