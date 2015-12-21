##test various curve-fitting metrics
library(nplr)
library(drc)
#library(drfit)
source("../../bin/drugSensData.R")

curveCompare<-function(datfile){

  res=apply(datfile,1,function(x){
    resp=as.numeric(x[grep('DATA[0-9]+',colnames(datfile))])
    nresp=convertToProp(resp)#as.numeric(resp/max(resp))
    conc=as.numeric(x[grep('^C[0-9]+',colnames(datfile))])
  
  ##1-- start with newton-raphson method from fred/brian
  n.res=nplr(conc,nresp)
  
  ##2-- DRC (r package cited in paper)
  df=data.frame(dose=conc,response=nresp,substance=rep(x[2],length(resp)),unit=1)
  names(df)[3]='substance'
  d.res=drcfit(df)
##  compare with points from GRID provided by NCATS
})
}

all.res=lapply(allfiles, curveCompare)

#now what? 