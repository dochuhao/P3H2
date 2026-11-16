library(limma)
library(ConsensusClusterPlus)
expFile="interGeneExp.txt"
workDir="."

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=1000,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=999999,
              plot="pdf")

clusterNum=3
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("DEGcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$DEGcluster))
cluster$DEGcluster=letter[match(cluster$DEGcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="DEGcluster.txt", sep="\t", quote=F, col.names=F)
