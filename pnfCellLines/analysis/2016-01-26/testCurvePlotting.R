##test curve-plotting from parameters


source("../../bin/drugSensData.R")

sapply(names(allfiles),function(x){
  print(x)
  df=allfiles[[x]]
  if("CRC"%in%colnames(df))
    drugs=df[which(df[,"CRC"]%in%c(-1.1,-1.2)),"name"]
  else
    drugs=df[which(df[,"CCLASS2"]%in%c(-1.1,-1.2)),"name"]  
  print(paste("found",length(drugs),'effective compounds for',x))
  if(length(drugs)>5){
    print("sampling only 5")
    drugs=sample(drugs,5)
  }
  sapply(drugs,function(y) doseResponseCurve(x,y))
})
