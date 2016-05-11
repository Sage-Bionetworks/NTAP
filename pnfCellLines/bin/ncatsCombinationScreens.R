##begin to evaluate combination screens....

##
## Drug sensitivity data files for NCATS single agent screens
##
##
library(synapseClient)
library(data.table)
library(pheatmap)
synapseLogin()
library(ggplot2)
require(reshape2)
screendirs=list(`10x10`='syn5611797',`6x6`='syn5611796')

getFileForCombo<-function(comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
    if(measure=='CCG'&&comboScreen=='6x6'){
      print('CCG data not available for 6x6, evaluating with CTG')
      measure='CTG'
    }
      
    fileparent=screendirs[[comboScreen]]  
    cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
    print(paste('Retrieved',nrow(cells),'cell types for',comboScreen,'screen'))
    cell.dat<-lapply(cells[,1],function(x){
      synid=cells[match(x,cells[,1]),2]
      allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
      allf<-allf[grep(measure,allf[,1]),]
      meta<-data.frame(fread(synGet(allf[grep('metadata',allf[,1]),2])@filePath))
      calc<-data.frame(fread(synGet(allf[grep('calc',allf[,1]),2])@filePath))
      resp<-data.frame(fread(synGet(allf[grep('resp',allf[,1]),2])@filePath))
      ##not sure which values to use?  
      combined<-cbind(meta,calc)
      ##return rowtrug, row target, col drug, col target, and some synergy score....
      return(combined)
    })
    names(cell.dat)<-cells[,1]
    return(cell.dat)
}



##how to plot values acros cells - experiment with the best way to visualize the data in its 
##initial form
plotValsAcrossCells<-function(file.list,prefix='',value='Beta'){
  fres<-do.call('rbind',lapply(names(file.list),function(y){
    x<-file.list[[y]]
    data.frame(Row=x$RowName,Col=x$ColName,
            RowTarg=x$RowTarget,ColTarg=x$ColTarget,
            Value=x[,value],Cell=y)
  }))
  
  
  fname=paste(prefix,value,'Values',sep='')
    
  ##plot values by cell
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Cell))
  ggsave(paste(fname,'ByCell.png',sep=''))
  
  ##plot values by row, column to see if there are outliers 
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Row,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByRow.png',sep=''))
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Col,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByColumn.png',sep=''))
  
  ##then row/column targets
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=RowTarg,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByRowTarget.png',sep=''))
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=ColTarg,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByColTarget.png',sep=''))
  return(fres)
}



doLinearModel<-function(value.df,prefix){
  
  lres<-lm(Value~Row+Col+Cell+Row*Col,value.df)
  dcoeffs<-data.frame(summary(lres)$coefficients)
  #print(head(coeffs))
  dcoeffs$adjustedP<-p.adjust(as.numeric(dcoeffs[,4]))
  dcoeffs<-dcoeffs[order(dcoeffs[,4],decreasing=F),]
  write.table(dcoeffs,paste(prefix,'DrugLinearModelCoefficients.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  
  lres<-lm(Value~RowTarg+ColTarg+Cell+RowTarg*ColTarg,value.df)
  tcoeffs<-data.frame(summary(lres)$coefficients)
  tcoeffs$adjustedP<-p.adjust(as.numeric(tcoeffs[,4]))
  tcoeffs<-tcoeffs[order(tcoeffs[,4],decreasing=F),]
  write.table(tcoeffs,paste(prefix,'DrugTargetLinearModelCoefficients.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  
  return(list(drug=dcoeffs,target=tcoeffs))
  
}

plotLMResults<-function(lmres,value.df,pval=0.05){
    combos<-lmres[grep(':',rownames(lmres)),]
    combos<-cbind(combos,t(sapply(rownames(combos),function(x) unlist(strsplit(gsub("Row|Col",'',x),split=':')))))
    colnames(combos)[6:7]<-c('Row','Col')    
    #pmat<-reshape2::dcast(combos,Row~Col,value.var='Pr...t..')
    apmat<-reshape2::dcast(combos,Row~Col,value.var='adjustedP')
    rownames(apmat)<-apmat$Row
    apmat<-apmat[,-1]
    apmat[which(is.na(apmat),arr.ind=T)]<-1.0
    
    mat<-reshape2::dcast(value.df,Row~Col,value.var='Value',fun.aggregate=mean)
    rownames(mat)<-mat$Row
    mat<-mat[,-1]
    mat[which(is.na(mat),arr.ind=T)]<-0.0
    rord<-order(apply(mat,1,sum))
    cord<-order(apply(mat,2,sum))
    
    par(mfrow=c(2,1))
    pheatmap(mat[rord,cord],cellheight=10,cellwidth=10,cluster_rows = F,cluster_cols = F)
}
    
    
    
    