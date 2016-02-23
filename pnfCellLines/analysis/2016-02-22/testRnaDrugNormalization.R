source("../../bin/singleDrugAnalysis.R")
source("../../bin/crossDataComps.R")

alpha.pars=c(1,4,5,10,100,1000)

this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-22/testRnaDrugNormalization.R'

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
    res=computeDrugRnaNormalizedCor(drugMat=ctrpMat,rnaMat=ccle.tpm,prefix='ctrpOriginal',sampleCombos=NA,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synapseStore(sf,used=list(url=this.script,wasExecuted=TRUE))
    }
      
  }
  
}

###build functions to identify ideal ALPHA parameter for transform
#' find best alpha parameter in NCATS data
#' @param alpha.params - list of parameters to sample
#' @return file names to upload to synapse
testNcatsNorm<-function(){
  ##ncats vs. cell lines
  source('../../bin/RNASeqData.R')
  source("../../bin/ncatsSingleAgentScreens.R")
  genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE)
  drugMat<-getValueForAllCells("FAUC")
  
  for(ap in alpha.pars){
    res=computeDrugRnaNormalizedCor(drugMat=drugMat,rnaMat=genCodeMat,prefix='ncatsOriginal',sampleCombos=NA,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synapseStore(sf,used=list(url=this.script,wasExecuted=TRUE))
    }
    
  }
}


testCTD2Norm()
testNcatsNorm()

