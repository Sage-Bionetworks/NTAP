#!/bin/bash

#test segmentation analysis
source("../../bin/segmentCNVData.R")

#test clustering
source('../../bin/clusterCNVBySample.R')

#now upload images for wiki

#image 1 - unsupervised clustering
sf=File('all_gene_by_median_logRRatios_dendro.png')
#add annotations
#add in activity
synStore(sf)

#image 2 - demonstration of gene variability
sf=File('most_variable_100_gene_by_median_logRRatios_heatmap.png')
#add annotations
#add in activity
synStore(sf)


#image 3 - correlated by genotype
sf=File('most_gt_correlated_100_gene_by_median_logRRatios_heatmap.png')
#add annotations
#add in activity
synStore(sf)
