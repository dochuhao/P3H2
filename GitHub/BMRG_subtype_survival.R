#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

clusterFile="BMcluster.txt"
cliFile="time.txt"

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

length=length(levels(factor(rt$BMcluster)))
diff=survdiff(Surv(futime, fustat) ~ BMcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ BMcluster, data = rt)
#print(surv_median(fit))

bioCol=c( "#E15759","#4E79A7")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit,
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=5.5,
		           legend.title="BMcluster",
		           legend.labs=levels(factor(rt[,"BMcluster"])),
		           legend = c(0.85, 0.85),
		           font.legend=14,
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette = bioCol,
		           #surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=6,height=5)
print(surPlot)
dev.off()
