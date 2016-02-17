
##double checking these scripts to ensure they map back to the BROAD data

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")

#plotMostVariableAUCs()


##first get original auc vals
orig.ncats<-getValueForAllCells("FAUC")
orig.ctrp<-getCtrpScreensAsMatrix()

#rescored.ncats<-getRecalculatedAUCMatrix()

rescored.ncats<-acast(getRecalculatedAUCTab(),Drug~Cell,value.var="AUC",fun.aggregate = mean)
rescored.ctrp<-acast(ctrpDoseResponseCurve(FALSE),Drug~Cell,value.var="AUC",fun.aggregate = mean)

##now test normalization on rows and columns
testNormalizationParameters(orig.ncats,byCol=FALSE,alphas=c(1,5,10,100),prefix='ncatsDrug')
testNormalizationParameters(orig.ncats,byCol=TRUE,alphas=c(1,5,10,100),prefix='ncatsCell')

testNormalizationParameters(rescored.ncats,byCol=FALSE,alphas=c(1,5,10,100),prefix='rescoredNcatsDrug')
testNormalizationParameters(rescored.ncats,byCol=TRUE,alphas=c(1,5,10,100),prefix='rescoredNcatsCell')

testNormalizationParameters(orig.ctrp,byCol=FALSE,alphas=c(1,5,10,100),prefix='ctrpDrug')
testNormalizationParameters(orig.ctrp,byCol=TRUE,alphas=c(1,5,10,100),prefix='ctrpCell')

testNormalizationParameters(rescored.ctrp,byCol=FALSE,alphas=c(1,5,10,100),prefix='rescoredCtrpDrug')
testNormalizationParameters(rescored.ctrp,byCol=TRUE,alphas=c(1,5,10,100),prefix='rescoredCtrpCell')

##now add in more cell parameters
testNormalizationParameters(orig.ncats,byCol=TRUE,alphas=c(20,40,60,80),prefix='ncatsCell')
testNormalizationParameters(rescored.ncats,byCol=TRUE,alphas=c(20,40,60,80),prefix='rescoredNcatsCell')

