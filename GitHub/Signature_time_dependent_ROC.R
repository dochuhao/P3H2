library(survival)
library(survminer)
library(timeROC)

bioROC=function(inputFile=null, rocFile=null){

	rt=read.table(inputFile, header=T, sep="\t", check.names=F)

	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)

	pdf(file=rocFile, width=2.8, height=2.8)

	par(oma=c(0,0,0,0))

	par(mar=c(4,4,0,0) + 0.1)
	plot(ROC_rt,time=1,col="#7F3C8D",title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col="#11A579",add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col="#F2B701",add=TRUE,title=FALSE,lwd=2)

	legend('bottomright',
	        c(paste0('1-year: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('3-year: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('5-year: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("#7F3C8D","#11A579","#F2B701"),lwd=2,bty = 'n')
	dev.off()
}

bioROC(inputFile="riskTrain.txt", rocFile="ROC.train.pdf")
bioROC(inputFile="riskTestGEO1.txt", rocFile="ROC.test1.pdf")
