source("../../bin/singleDrugAnalysis.R")
source("../../bin/crossDataComps.R")

alpha.pars=c(1,4,5,10,100,1000)

###build functions to identify ideal ALPHA parameter for transform
#' find best alpha parameter in CTD2/CCLE data
#' @param alpha.params - list of parameters to sample
#' @return file names to upload to synapse
testCTD2Norm<-function(){
  source("../../bin/ccleData.R")
  source("../../bin/ctrpSingleAgentScreens.R")
  ##ctrp vs. ccle
  ctrpMat=as.data.frame(getCtrpScreensAsMatrix())
  ccle.tpm<-getCCLEDataTPM(removeDupes=TRUE)

  for(ap in alpha.pars){
    res=computeDrugRnaNormalizedCor(drugMat=ctrpMat,rnaMat=ccle.tpm,prefix='ctrpOriginal',sampleCombos=10000,alpha=5)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synapseStore(sf,used=list())
    }
      
  }
  
}

###build functions to identify ideal ALPHA parameter for transform
#' find best alpha parameter in NCATS data
#' @param alpha.params - list of parameters to sample
#' @return file names to upload to synapse
testNcatsNorm<-function(){
  ##ncats vs. cell lines
  genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
  drugMat<-getValueForAllCells("FAUC")
}




