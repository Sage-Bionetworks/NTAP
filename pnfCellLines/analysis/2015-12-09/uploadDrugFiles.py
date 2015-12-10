#/usr/bin/python
import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()

fd='NTAP MIPE qHTS Data Dec2015'
allfiles=os.listdir(fd)

for f in allfiles:
    celltype=re.sub("NTAP ","",f.split(' MIPE ')[0])
    #first assign genotype, sampleName, sampleOrigin
    sct=celltype.split()
    if sct[0]=='ipNF02.3':
        gt='++'
        sn='ipn02.3'
        so='pn02.3'
    elif sct[0]=='ipNF95.11C':
        gt='+-'
        sn='ipnNF95.11c'
        so='pnNF95.11c'
    else:
        gt='--'
        if len(sct)==1:
            sn=re.sub('C','c',sct[0])
            so=sn[1:]
        else:
            if sct[1]=='MC':
                sn=sct[0]+' (mixed clone)'
                so=sct[0][1:]
            elif sct[1]=='SC':
                sn=sct[0]+' (single clone)'
                so=sct[0][1:]
            elif sct[1]=='C_T':
                sn=sct[0]+'C'
                so=sn[1:]
            else:
                print 'Need to annotate '+celltype
    print 'File: %s, Genotype: %s, Sample: %s, Origin: %s'%(f,gt,sn,so)
    nf=synapseclient.entity.File(os.path.join(fd,f),parent='syn5522627',annotations={'sampleGenotype':gt,'sampleName':sn,'sampleOrigin':so})
    syn.store(nf)
