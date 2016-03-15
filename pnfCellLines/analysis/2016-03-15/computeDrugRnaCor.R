source("../../bin/singleDrugAnalysis.R")
source("../../bin/crossDataComps.R")

#alpha.pars=c(1,4,5,10,100,1000)

this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-03-15/computeDrugRnaCor.R'

run.all.norm<-function(){
    require(parallel)
    source("../../bin/ccleData.R")
    source("../../bin/ctrpSingleAgentScreens.R")
    ##ctrp vs. ccle
     ctrpMat=as.data.frame(getCtrpScreensAsMatrix())
    ccle.tpm<-getCCLEDataTPM(removeDupes=TRUE)

                                        #also get re-normalized CTRP
    ctrp.reMat<-ctrpDoseResponseCurve(FALSE,TRUE)
    ctrp.alpha=10

    source('../../bin/RNASeqData.R')
    source("../../bin/ncatsSingleAgentScreens.R")
    genCodeMat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)
    ncats.drugMat<-getValueForAllCells("FAUC")
    ncats.remat=getRecalculatedAUCMatrix()

    ##now create tuples
    run.list<-list(ctrp.original=list(exp=ccle.tpm,drug=ctrpMat,alpha=ctrp.alpha,name='ctrpOriginal'),
                   ctrp.rescored=list(exp=ccle.tpm,drug=ctrp.reMat,alpha=ctrp.alpha,name='ctrpRescored'),
                   ncats.original=list(exp=genCodeMat,drug=ncats.drugMat,alpha=1,name='ncatsOriginal'),
                   ncats.rescored=list(exp=genCodeMat,drug=ncats.remat,alpha=1,name='ncatsRescored'))
    all.files=mclapply(run.list,function(x){
        computeDrugRnaNormalizedCor(drugMat=x$drug,rnaMat=x$exp,prefix=x$name,sampleCombos=NA,alpha=x$alpha)
    },mc.cores=4)

    for(f in all.files){
        sf=File(f,parentId='syn5728996')
        synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE)))
  }

}


run.all.norm()
