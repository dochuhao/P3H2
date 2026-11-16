library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="mRNA_tpm.txt"
clusterFile="DEGcluster.txt"
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=log2(data+1)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data,
                geneSets,
                min.sz=10,
                max.sz=500,
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
#Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
#gsvaCluster=cbind(gsvaCluster, Project)

adj.P.Val.Filter=0.24
allType=as.vector(gsvaCluster$DEGcluster)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){

	treat=gsvaCluster[gsvaCluster$DEGcluster==comp[2,i],]
	con=gsvaCluster[gsvaCluster$DEGcluster==comp[1,i],]
	data=rbind(con, treat)

	Type=as.vector(data$DEGcluster)
	ann=data[,c(ncol(data), (ncol(data)-1))]
	data=t(data[,-c((ncol(data)-1), ncol(data))])
	design=model.matrix(~0+factor(Type))
	colnames(design)=levels(factor(Type))
	fit=lmFit(data, design)
	contrast=paste0(comp[2,i], "-", comp[1,i])
	cont.matrix=makeContrasts(contrast, levels=design)
	fit2=contrasts.fit(fit, cont.matrix)
	fit2=eBayes(fit2)

	allDiff=topTable(fit2,adjust='fdr',number=200000)
	allDiffOut=rbind(id=colnames(allDiff),allDiff)
	write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)

	diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
	diffSigOut=rbind(id=colnames(diffSig),diffSig)
	write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)

	bioCol=c("#11A579", "#F2B701","#7F3C8D","#80BA5A","#3969AC", "#E73F74", "#11A579", "#80BA5A","#7F3C8D",  "#A5AA99")
	ann_colors=list()
	m6aCluCol=bioCol[1:length(levels(factor(allType)))]
	names(m6aCluCol)=levels(factor(allType))
	ann_colors[["DEGcluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]

	termNum=10
	diffTermName=as.vector(rownames(diffSig))
	diffLength=length(diffTermName)
	if(diffLength<termNum){termNum=diffLength}
	hmGene=c("KEGG_CIRCADIAN_RHYTHM_MAMMAL","KEGG_PRIMARY_IMMUNODEFICIENCY",
	         "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
	         "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_ARGININE_AND_PROLINE_METABOLISM",
	         "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS","KEGG_FATTY_ACID_METABOLISM",
	         "KEGG_SULFUR_METABOLISM","KEGG_STEROID_BIOSYNTHESIS")
	#hmGene=diffTermName[1:termNum]
	hmExp=data[hmGene,]
	pdf(file=paste0(contrast,".heatmap.pdf"), width=12, height=3)
	pheatmap(hmExp,
	         annotation=ann,
	         annotation_colors = ann_colors,
	         color = colorRampPalette(c(rep("#59A14F",2), "white", rep("#F28E2B",2)))(50),
	         cluster_cols =F,
	         show_colnames = F,
	         gaps_col=as.vector(cumsum(table(Type))),
	         scale="row",
	         fontsize = 10,
	         fontsize_row=10,
	         fontsize_col=10)
	dev.off()
}
