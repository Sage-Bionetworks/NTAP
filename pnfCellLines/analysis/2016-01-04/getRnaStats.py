'''
Script to run picard tools on belltown to analyze bam files output from kallisto
'''

import os,sys,re


##module load picard-tools
picard_home='/gluster/toolbox/picard-tools/1.94/'

#these are the tools to run
rs_met='CollectRnaSeqMetrics.jar'
as_metrics='CollectAlignmentSummaryMetrics.jar'

bamdir='../2015-12-22'

for file in os.listfiles(bamdir):
    fullfile=os.path.join(bamdir,file)
    rs_outfile=re.sub('.bam','_rnaseqMetrics.txt',file)

    as_outfile=re.sub('.bam','_alignMetrics.txt',file)
