library(survival)
library(survminer)

bioSUR=function(inputFile=null, surFile=null){
  rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
  res.cut<- surv_cutpoint(rt, time="futime", event="fustat",variables = c("riskScore"))
  res.cat<- surv_categorize(res.cut)

  riskScore<-survfit(Surv(futime,fustat)~riskScore,data = res.cat)
  pdf(file=surFile, onefile=FALSE, width=3, height=2)
  surPlot=ggsurvplot(riskScore,
                     data = res.cat,
                     conf.int=F,
                     risk.table = F,
                     pval=T,
                     pval.size=6,
                     legend.title="Riskscore",
                     legend.labs=levels(factor(res.cat[,"riskScore"])),
                     legend = c(1, 0.5),
                     font.legend=15,
                     xlab=" ",
                     palette=c("#EDC948","#59A14F"))
  print(surPlot)
  dev.off()
}

bioSUR(inputFile="TrainACTA2.txt", surFile="ACTA2.train.pdf")
bioSUR(inputFile="TrainFRAS1.txt", surFile="FRAS1.train.pdf")
bioSUR(inputFile="TrainP3H2.txt", surFile="P3H2.train.pdf")
bioSUR(inputFile="TestACTA2.txt", surFile="ACTA2.test.pdf")
bioSUR(inputFile="TestFRAS1.txt", surFile="FRAS1.test.pdf")
bioSUR(inputFile="TestP3H2.txt", surFile="P3H2.test.pdf")
