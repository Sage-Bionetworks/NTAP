'''
Run kallisto on the deposited FASTQ files
'''

import os,sys,re

fqdir='/scratch/NTAP/rnaseq/hhao_123771/FASTQ/'

allfqfiles=[a for a in os.listdir(fqdir) if 'gz' in a]

allsamps=set()
[allsamps.add('_'.join(a.split('_')[0:2])) for a in allfqfiles]

for a in allsamps:
    f1=os.path.join(fqdir,a+'_0_1.fastq.gz')
    f2=os.path.join(fqdir,a+'_0_2.fastq.gz')
    kcmd='~/kallisto_linux-v0.42.4/kallisto quant -i Homo_sapiens.GRCh38.rel79.idx '+f1+' '+f2+' -o '+a
    print kcmd
    os.system(kcmd)
