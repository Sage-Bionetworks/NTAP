source("../../bin/singleDrugAnalysis.R")
source("../../bin/crossDataComps.R")
require(parallel)
alpha.pars=c(1,4,5,10,100,1000)

this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-22/testRnaDrugNormalization.R'

###build functions to identify ideal ALPHA parameter for transform
#' find best alpha parameter in CTD2/CCLE data
#' @param alpha.params - list of parameters to sample
#' @return file names to upload to synapse
testCTD2Norm<-function(numSamps=NA){
  source("../../bin/ccleData.R")
  source("../../bin/ctrpSingleAgentScreens.R")
  ##ctrp vs. ccle
  ctrpMat=as.data.frame(getCtrpScreensAsMatrix())
  ccle.tpm<-getCCLEDataTPM(removeDupes=TRUE)

  #also get re-normalized CTRP
  reMat<-ctrpDoseResponseCurve(FALSE,TRUE)
 
  #for(ap in alpha.pars){
  r=mclapply(alpha.pars[-1],function(ap){

    res=computeDrugRnaNormalizedCor(drugMat=ctrpMat,rnaMat=ccle.tpm,prefix='ctrpOriginal',sampleCombos=numSamps,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE)))
    }

    res=computeDrugRnaNormalizedCor(drugMat=reMat,rnaMat=ccle.tpm,prefix='ctrpRescored',sampleCombos=numSamps,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE)))
    }

  })

}

###build functions to identify ideal ALPHA parameter for transform
#' find best alpha parameter in NCATS data
#' @param alpha.params - list of parameters to sample
#' @return file names to upload to synapse
testNcatsNorm<-function(numSamps=NA){
  ##ncats vs. cell lines
  source('../../bin/RNASeqData.R')
  source("../../bin/ncatsSingleAgentScreens.R")
  genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)
  drugMat<-getValueForAllCells("FAUC")
  remat=getRecalculatedAUCMatrix()
 # for(ap in alpha.pars){
  r=mclapply(alpha.pars,function(ap){
    res=computeDrugRnaNormalizedCor(drugMat=drugMat,rnaMat=genCodeMat,prefix='ncatsOriginal',sampleCombos=numSamps,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE)))
    }

    res=computeDrugRnaNormalizedCor(drugMat=reMat,rnaMat=genCodeMat,prefix='ncatsRescored',sampleCombos=numSamps,alpha=ap)
    for(f in res){
      sf=File(f,parentId='syn5679539')
      synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE)))
    }
  })
}


testCTD2Norm(500000)
testNcatsNorm(500000)
