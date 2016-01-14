'''
Run kallisto on the deposited FASTQ files
'''

import os,sys,re

fqdir='/scratch/NTAP/rnaseq/hhao_123771/FASTQ/'

allfqfiles=[a for a in os.listdir(fqdir) if 'gz' in a]

allsamps=set()
[allsamps.add('_'.join(a.split('_')[0:2])) for a in allfqfiles]

##first make new index
iname='gencode.v24.transcripts.idx'
icmd='~/kallisto_linux-v0.42.4/kallisto index gencode.v24.transcripts.fa.gz -i '+iname
if not os.path.exists(iname):	
    print icmd
    os.system(icmd)

for a in allsamps:
    f1=os.path.join(fqdir,a+'_0_1.fastq.gz')
    f2=os.path.join(fqdir,a+'_0_2.fastq.gz')
    ofile=a+'.out'
    kcmd='~/kallisto_linux-v0.42.4/kallisto quant --bootstrap-samples=100 --threads=16 -i '+iname+' '+f1+' '+f2+' -o '+a+' 2>'+ofile
    print kcmd
    os.system(kcmd)
    
