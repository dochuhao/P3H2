library(survival)
library(survminer)

bioSUR=function(inputFile=null, surFile=null){
  rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
  res.cut<- surv_cutpoint(rt, time="futime", event="fustat",variables = c("riskScore"))
  res.cat<- surv_categorize(res.cut)

  riskScore<-survfit(Surv(futime,fustat)~riskScore,data = res.cat)
  pdf(file=surFile, onefile=FALSE, width=3.2, height=2.8)
  surPlot=ggsurvplot(riskScore,
                     data = res.cat,
                     conf.int=F,
                     risk.table = F,
                     pval=T,
                     pval.size=6,
                     legend.title="Riskscore",
                     legend.labs=levels(factor(res.cat[,"riskScore"])),
                     legend = c(0.7, 0.85),
                     font.legend=15,
                     xlab="Time(years)",
                     palette=c("#EDC948","#59A14F"))
  print(surPlot)
  dev.off()
}

bioSUR(inputFile="riskTrain.txt", surFile="sur.train.pdf")
bioSUR(inputFile="riskTestGEO1.txt", surFile="sur.test1.pdf")
