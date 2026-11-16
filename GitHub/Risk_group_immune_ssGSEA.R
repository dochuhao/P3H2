library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="mRNA_tpm.txt"
riskFile="risk.txt"
gmtFile="immune.gmt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=log2(data+1)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

# GSVA >= 1.50 uses method-specific parameter objects.
if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  ssgseaParams = GSVA::ssgseaParam(
    exprData = data,
    geneSets = GSEABase::geneIds(geneSets),
    alpha = 0.25,
    normalize = TRUE
  )
  ssgseaScore = GSVA::gsva(ssgseaParams, verbose = FALSE)
} else {
  ssgseaScore = GSVA::gsva(
    data, geneSets, method = "ssgsea",
    tau = 0.25, ssgsea.norm = TRUE, verbose = FALSE
  )
}

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(risk))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
scorerisk=cbind(ssgseaScore, risk)

data=melt(scorerisk, id.vars=c("risk"))
colnames(data)=c("risk", "Immune", "Fraction")

bioCol=c("#E15759","#4E79A7","#F28E2B","#76B7B2","#59A14F","#EDC948","#B07AA1")
bioCol=bioCol[1:length(levels(factor(data[,"risk"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="risk",
     ylab="Immune infiltration",
     xlab="",
     legend.title="Risk",
     palette=bioCol)
p=p+rotate_x_text(50)

pdf(file="boxplot.pdf", width=8, height=4.5)
p+stat_compare_means(aes(group=risk), method="wilcox.test", hide.ns=FALSE, label.y=max(data$Fraction, na.rm=TRUE)*1.08, vjust=0.5,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()
