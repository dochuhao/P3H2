library(survival)
library(survminer)

clusterFile="DEGcluster.txt"
cliFile="time.txt"

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

length=length(levels(factor(rt$DEGcluster)))
diff=survdiff(Surv(futime, fustat) ~ DEGcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ DEGcluster, data = rt)
#print(surv_median(fit))

bioCol=c("#F28E2B","#59A14F","#E15759","#4E79A7","#76B7B2","#EDC948","#B07AA1")
bioCol=bioCol[1:length(levels(factor(rt[,"DEGcluster"])))]
surPlot=ggsurvplot(fit,
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="DEGcluster",
		           legend.labs=levels(factor(rt[,"DEGcluster"])),
		           legend = c(0.85, 0.5),
		           font.legend=13,
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=F,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf", onefile = FALSE, width=5, height=3.5)
print(surPlot)
dev.off()
