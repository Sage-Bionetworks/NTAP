##get significanly correlated cell lines

require(synapseClient)

synapseLogin()

res=synapseQuery("select * from entity where parentId=='syn5594111'")

for(i in 1:nrow(res)){
  fname=res[[i,'entity.name']]
  synid=res[[i,'entity.id']]
  tab<-read.table(synGet(synid)@filePath)
  
}