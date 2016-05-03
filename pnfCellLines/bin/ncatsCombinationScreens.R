##begin to evaluate combination screens....

##
## Drug sensitivity data files for NCATS single agent screens
##
##
library(synapseClient)
library(data.table)
library(pheatmap)
synapseLogin()

screendirs=list(`10x10`='syn5611797',`6x6`='syn5611796')

getFileForCombo<-function(comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
    fileparent=screendirs[[comboScreen]]  
    cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
    print(paste('Retrieved',nrow(cells),'cell types for',comboScreen,'screen'))
    cell.dat<-lapply(cells[,1],function(x){
      synid=cells[match(x,cells[,1]),2]
      allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
      meta<-data.frame(fread(synGet(allf[grep('metadata',allf[,1]),2])@filePath))
      calc<-data.frame(fread(synGet(allf[grep('calc',allf[,1]),2])@filePath))
      ##not sure which values to use?  
      
      ##return rowtrug, row target, col drug, col target, and some synergy score....
    })
      
}