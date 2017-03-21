import synapseclient

syn = synapseclient.login()

def writeWholeTable(table):
	rowset = table.asRowSet()
	wikiTable =  "|".join([head['name'] for head in rowset['headers']]) + "\n"
	wikiTable =  wikiTable + "|".join(["---"]*len(rowset['headers'])) + "\n"
	for row in rowset['rows']:
		row['values'] = map(str, row['values'])
		wikiTable  = wikiTable + "|".join(row['values']) + "\n"
	return(wikiTable)


def writeProjectTables(table, getProjectIndex):
	"""
	table: 			  table query object
	getProjectIndex:  index of rowSet where you want to display the project name and synapse id link
	"""
	rowset = table.asRowSet()
	wikiTable =  "|".join([head['name'] for head in rowset['headers']]) + "\n"
	wikiTable =  wikiTable + "|".join(["---"]*len(rowset['headers'])) + "\n"
	for row in rowset['rows']:
		row['values'] = map(str, row['values'])
		try:
			projName =  syn.get(row['values'][getProjectIndex]).name
			row['values'][getProjectIndex] = "[%s](%s)" % (projName, row['values'][getProjectIndex])
		except Exception as e:
			print(e)
		wikiTable  = wikiTable + "|".join(row['values']) + "\n"
	return(wikiTable)

##Projects   syn4939478/wiki/235831

# table = syn.tableQuery('SELECT Name, Project_Title, Tumor_Focus, Synapse_ID FROM syn5867440 Where Active = True')
# firstTable = writeProjectTables(table, -1)

# table = syn.tableQuery('SELECT * FROM syn7239496')
# secondTable = writeProjectTables(table, 1)

# markdown = ("Many NTAP-funded initiatives already have data that are being uploaded to Synapse so they can be downloaded.  "
#  "Below is a table of projects with data.\n\n%s\n"
#  "In addition to the NTAP-funded grants, NTAP is working closely with [Sage Bionetworks](http://www.sagebase.org) "
#  "to carry out systems-level analysis on NF-related datasets.  "
#  "Those projects can be browsed below.\n\n%s\n"
#  "To view a complete list of NTAP projects, go to [NTAP Projects](https://www.synapse.org/#!Synapse:syn5867440/)" ) % (firstTable, secondTable)

# wikipage = syn.getWiki(syn.get('syn4990358'),"411306")
# wikipage.markdown = markdown
# syn.store(wikipage)

### Data upload   syn4939478/wiki/411657

table = syn.tableQuery('select projectEntity, numberOfFiles, numberOfContributors, lateModified from syn7804884 where Active = True order by "lateModified" DESC')
firstTable = writeProjectTables(table, 0)
#table = syn.tableQuery('select * from syn7805078')
#secondTable = writeWholeTable(table)

markdown = "Here is a summary of the latest activity by project:\n%s\n\n" % firstTable


 #"There are numerous types of data available:\n\n%s") % (firstTable, secondTable)

wikipage = syn.getWiki(syn.get('syn4939478'),"411657")
wikipage.markdown = markdown
syn.store(wikipage)