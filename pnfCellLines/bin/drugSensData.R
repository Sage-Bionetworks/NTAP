##
## Drug sensitivity data file
##
##
##
library(synapseClient)

fileparent='syn5522627'

headerfile='syn5522652'
qr<-synapseQuery(paste("select * from entity where entity.parentId==",fileparent,sep=''))
