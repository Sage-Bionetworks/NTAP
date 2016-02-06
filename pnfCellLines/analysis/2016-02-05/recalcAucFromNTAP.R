##compare to CTP dataa
source('../../bin/drugSensData.R')


require(plyr)
require(nplr)


all.aucs=sapply(dfiles$entity.sampleName,function(x) doseResponseCurve(x,NA,FALSE))

##now try to reshape into matrix
require(reshape2)
drugmat<-acast(new.df,Drug~CCL,value.var="AUC",fun.aggregate=function(x) mean(x,na.rm=T))#,fill=NA)

##now try to match to NTAP data
source("../../bin/drugSensData.R")
fauc.vals=getValueForAllCells("FAUC")
nt.drugs<-rownames(fauc.vals)
matched.drug.names=sapply(rownames(drugmat),function(x){
  mval=NA
  mval=match(x,nt.drugs)
  if(!is.na(mval))
    return(mval)

  y=paste('^',x,'$',sep='')
  mval=grep(y,nt.drugs,ignore.case=T)
  if(length(mval)==1)
    return(mval)

  mval=grep(gsub('-','',x),nt.drugs,ignore.case=T)
  if(length(mval)==1)
    return(mval)

  return(NULL)
})

drug.matches=unlist(matched.drug.names)
print(paste("Found",length(drug.matches),'drugs that are found in both NCATS and CCL screens'))
#create a combined matrix

##first do some sort of z-scoring...
comb.mat=cbind(fauc.vals[drug.matches,],
               drugmat[match(names(drug.matches),rownames(drugmat)),])


missing=which(apply(comb.mat,1,function(x) length(which(is.nan(x))))>100)
print(paste("Removing",length(missing),'drugs with too many missing values'))
comb.mat=comb.mat[-missing,]


comb.norm<-apply(comb.mat,2,function(x) (x-mean(x,na.rm=T))/(sd(x,na.rm=T)))
comb.norm[which(is.na(comb.norm),arr.ind=T)]<-0.0
dmat=dist(t(comb.norm))

##now cluster
h=hclust(dmat)

##now do the clustering?
primsite=ccl_metadat$ccle_primary_site[match(colnames(drugmat),ccl_metadat$ccl_name)]
names(primsite)=colnames(drugmat)
primsite=c(primsite,sapply(colnames(fauc.vals),function(x) return("pNFs")))

require(pheatmap)
pheatmap(comb.norm,cellheight=10,cellwidth=10,filename='allAucZscoreByCell.png')

ccors=cor(comb.norm)

most.cor= union(colnames(ccors)[1:8],names(sort(apply(ccors[1:8,],2,median,na.rm=T),decreasing=T)[1:80]))
pheatmap(comb.norm[,most.cor],clustering_method='ward',cellheight=10,cellwidth=10,annotation_col=data.frame(PrimarySite=as.factor(primsite)),
         filename='top80CorrelatedAucZscores.png')

most.cor= union(colnames(ccors)[1:8],names(sort(apply(ccors[1:8,],2,median,na.rm=T),decreasing=T)[1:50]))
pheatmap(comb.norm[,most.cor],clustering_method='ward',cellheight=10,cellwidth=10,annotation_col=data.frame(PrimarySite=as.factor(primsite)),
         filename='top50CorrelatedAucZscores.png')

most.cor= union(colnames(ccors)[1:8],names(sort(apply(ccors[1:8,],2,median,na.rm=T),decreasing=T)[1:20]))
pheatmap(comb.norm[,most.cor],cellheight=10,cellwidth=10,clustering_method='ward.D2',clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',annotation_col=data.frame(PrimarySite=primsite),
         filename='top20CorrelatedAucZscores.png')
