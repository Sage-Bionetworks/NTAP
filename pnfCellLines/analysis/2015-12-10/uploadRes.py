'''
upload results
'''
import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()

allplots=os.listdir('./')

rfi='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/drugSensData.R'

pid='syn5522627'
allsi=syn.query("select sampleName,id from entity where parentId=='%s'"%(pid))['results']

ids=[a['entity.id'] for a in allsi if a['entity.sampleName'] is not None]
redsi=[a for a in allsi if a['entity.id'] in ids]

for fi in allplots:
    if 'png' not in fi:
        continue
    elif 'CRCValues' in fi: #need toa ssociate with single sid
        cellname=re.sub(".png",'',fi.split('For')[1])
        sival=''
        sival=[a['entity.id'] for a in redsi if a['entity.sampleName'][0]==cellname]
        #print cellname,sival
        sivals=sival
    elif 'drugsWithMostVariable' in fi:
        sivals=ids
    newf=synapseclient.entity.File(fi,parentId='syn5522801')
    print fi,sivals
    syn.store(newf,used=sivals,executed=rfi)
