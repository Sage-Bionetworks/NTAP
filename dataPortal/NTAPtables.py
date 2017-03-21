import synapseclient
from synapseclient import Table
import pandas as pd
syn = synapseclient.login()

#Table 1:
#Project Entity -- Number of Files uploaded -- Number of contributors  -- Date of latest upload
#> synapseid, number of files total to the project, number of contributors, last modified on date
#https://www.synapse.org/#!Synapse:syn4939478/wiki/411657 
allNTAPprojects = syn.tableQuery('SELECT * FROM syn5867440')
allNTAPprojects_df = allNTAPprojects.asDataFrame()
allNTAPprojects_df = allNTAPprojects_df.drop_duplicates("Synapse_ID")

projectUploadActivitySynId = "syn7804884"
projectTracker = syn.tableQuery('select * from %s' % projectUploadActivitySynId)
projectTrackerDf = projectTracker.asDataFrame()
projectTrackerDf['lateModified'] = projectTrackerDf['lateModified'].fillna(0)
removeSamples = []
for index, synId in enumerate(allNTAPprojects_df['Synapse_ID']):
	temp = syn.chunkedQuery('select id, createdByPrincipalId, modifiedOn from file where projectId == "%s"' % synId)
	modifiedOn = []
	createdBy = []
	count = 0
	for x in temp:
		modifiedOn.append(x['file.modifiedOn'])
		createdBy.append(x['file.createdByPrincipalId'])
		count = count + 1
	if len(modifiedOn) > 0:
		mod = max(modifiedOn)
	else:
		mod = 0
	oldValues = projectTrackerDf[projectTrackerDf['projectEntity'] == synId]
	newValues = [synId, count, len(set(createdBy)), mod, allNTAPprojects_df['Active'][index]]
	if not oldValues.empty:
		if not all([old == new for old, new in zip(oldValues.values[0], newValues)]):
			projectTrackerDf[projectTrackerDf['projectEntity'] == synId] = newValues
		else:
			removeSamples.append(synId)
	else:
		projectTrackerDf = projectTrackerDf.append(pd.DataFrame([newValues],columns=['projectEntity','numberOfFiles','numberOfContributors','lateModified','Active']))

newUploads = projectTrackerDf[~projectTrackerDf['projectEntity'].isin(removeSamples)]
if not newUploads.empty:
	newUploads['lateModified'] = newUploads['lateModified'].apply(int)
	newUploads['numberOfFiles'] = newUploads['numberOfFiles'].apply(int)
	newUploads['numberOfContributors'] = newUploads['numberOfContributors'].apply(int)
	newUploads['lateModified'][newUploads['lateModified'] == 0] = ""
	schema = syn.get(projectUploadActivitySynId)
	tablestore = Table(schema, newUploads, etag=projectTracker.etag)
	tablestore = syn.store(tablestore)
else:
	print("No updates!")


#Table 2: Files by assay type
#Assay Type -- Number of Files -- Number of Cell Lines
#> assay, grab number of unique synapseid, sampleIdentifier
#https://www.synapse.org/#!Synapse:syn4939478/wiki/411658
ntap_generated_data_synId = "syn7805078"
ntap_generated_data = syn.tableQuery('SELECT * FROM %s' % ntap_generated_data_synId)
ntap_generated_data_df = ntap_generated_data.asDataFrame()

annot_synIds = ["syn7506024","syn7805075","syn7992153"]
assaysNumSynId = {}
assaysNumSampleId = {}

for synId in annot_synIds:
	annotations = syn.tableQuery('SELECT * FROM %s' % synId)
	annot_table = annotations.asDataFrame()
	assays = set(annot_table['assay'])
	for i in assays:
		numAssay = len(set(annot_table['synapseId'][annot_table['assay'] == i]))
		numId = len(set(annot_table['sampleIdentifier'][annot_table['assay'] == i]))
		if assaysNumSynId.get(i) is None:
			assaysNumSynId[i] = numAssay
		else:
			assaysNumSynId[i] = assaysNumSynId[i] + numAssay
		if assaysNumSampleId.get(i) is None:
			assaysNumSampleId[i] = numId
		else:
			assaysNumSampleId[i] = assaysNumSampleId[i] + numId

removeAssay = []
for i in assaysNumSynId.keys():
	oldValues = ntap_generated_data_df[ntap_generated_data_df['assay'] == i]
	newValues = {'assay':i, 'numberOfFiles':assaysNumSynId[i], 'numberOfCellLines':assaysNumSampleId[i]}
	if oldValues.empty:
		ntap_generated_data_df = ntap_generated_data_df.append(pd.DataFrame(newValues, index=[0]))
	else:
		if not all([oldValues[old].values[0] == newValues[old] for old, new in zip(oldValues, newValues)]):
			ntap_generated_data_df['numberOfCellLines'][ntap_generated_data_df['assay'] == i] = newValues['numberOfCellLines']
			ntap_generated_data_df['numberOfFiles'][ntap_generated_data_df['assay'] == i] = newValues['numberOfFiles']
		else:
			removeAssay.append(i)
		
newUploads = ntap_generated_data_df[~ntap_generated_data_df['assay'].isin(removeAssay)]
if not newUploads.empty:
	newUploads['numberOfFiles'] = newUploads['numberOfFiles'].apply(int)
	newUploads['numberOfCellLines'] = newUploads['numberOfCellLines'].apply(int)
	schema = syn.get(ntap_generated_data_synId)
	tablestore = Table(schema, newUploads, etag=ntap_generated_data.etag)
	tablestore = syn.store(tablestore)
else:
	print("No updates!")

