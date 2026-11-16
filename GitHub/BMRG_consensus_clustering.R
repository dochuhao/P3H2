#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
expFile="uniSigGeneExp.txt"
workDir="."

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

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

clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("BMcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$BMcluster))
cluster$BMcluster=letter[match(cluster$BMcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="BMcluster.txt", sep="\t", quote=F, col.names=F)
