options(stringsAsFactors=FALSE)
src <- Sys.getenv('IMMUNE_INPUT_DIR')
stopifnot(nzchar(src), dir.exists(src))
out <- Sys.getenv('IMMUNE_OUTPUT_DIR')
stopifnot(nzchar(out))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
canon<-function(x)sub('^TARGET[-_]','',x)
read_scores<-function(file){
 x<-as.matrix(read.delim(file,check.names=FALSE,row.names=1));storage.mode(x)<-'double'
 colnames(x)<-canon(colnames(x));x
}
groups<-read.delim(file.path(src,'41.ssGSEA/risk.txt'),row.names=1,check.names=FALSE)
rownames(groups)<-canon(rownames(groups));ids<-rownames(groups)
stopifnot(length(ids)==85)
cat('Risk-group counts:');print(table(groups$risk))
comparison<-function(mat,group,panel,test){
 ans<-lapply(seq_len(ncol(mat)),function(j){
  y<-mat[,j]
  pp<-if(sd(y)==0)1 else if(test=='Wilcoxon')wilcox.test(y~group,exact=FALSE)$p.value else kruskal.test(y~group)$p.value
  data.frame(panel=panel,feature=colnames(mat)[j],test=test,n=length(y),P=pp,
    high_median=if(test=='Wilcoxon')median(y[group=='high']) else NA_real_,
    low_median=if(test=='Wilcoxon')median(y[group=='low']) else NA_real_)
 })
 z<-do.call(rbind,ans);z$BH_FDR<-p.adjust(z$P,'BH');z
}
s6b<-read_scores(file.path(src,'41.ssGSEA/ssGSEA.result.txt'))
cat('Figure6B rows:',nrow(s6b),'\n');print(rownames(s6b))
r6b<-comparison(t(s6b[,ids,drop=FALSE]),groups$risk,'Figure6B','Wilcoxon')
xall<-read_scores(file.path(src,'42 Xcell/xCell_mRNA_tpm_xCell_0400081323.txt'))
xtypes<-setdiff(rownames(xall),c('ImmuneScore','StromaScore','MicroenvironmentScore'))
r6c<-comparison(t(xall[xtypes,ids,drop=FALSE]),groups$risk,'Figure6C_all64','Wilcoxon')
xshown<-read_scores(file.path(src,'42 Xcell/xCellsig.txt'))
cat('Figure6C shown rows:',nrow(xshown),'\n');print(rownames(xshown))
r6c$shown_in_Figure6C<-r6c$feature%in%rownames(xshown)
xcom<-read_scores(file.path(src,'42 Xcell/xCellscore.txt'))
r6d<-comparison(t(xcom[,ids,drop=FALSE]),groups$risk,'Figure6D','Wilcoxon')
est<-read.delim(file.path(src,'43.estimate/TMEscores.txt'),row.names=1,check.names=FALSE)
rownames(est)<-canon(rownames(est))
r6e<-comparison(as.matrix(est[ids,c('StromalScore','ImmuneScore','ESTIMATEScore')]),groups$risk,'Figure6E','Wilcoxon')
ss3<-read_scores(file.path(src,'25.ssGSEACor/ssgseaOut.txt'))
cl<-read.delim(file.path(src,'25.ssGSEACor/DEGcluster.txt'),row.names=1,check.names=FALSE)
rownames(cl)<-canon(rownames(cl));g<-factor(cl[ids,1],levels=c('A','B','C'))
cells<-c('aDCs','B_cells','CD8+_T_cells','DCs','iDCs','Macrophages','Mast_cells','Neutrophils','NK_cells','pDCs','T_helper_cells','Tfh','Th1_cells','Th2_cells','TIL','Treg')
functions<-setdiff(rownames(ss3),cells)
r3g<-comparison(t(ss3[cells,ids,drop=FALSE]),g,'Figure3G','Kruskal-Wallis')
r3h<-comparison(t(ss3[functions,ids,drop=FALSE]),g,'Figure3H','Kruskal-Wallis')
results<-list(Figure6B=r6b,Figure6C=r6c,Figure6D=r6d,Figure6E=r6e,Figure3G=r3g,Figure3H=r3h)
for(nm in names(results)){
 write.csv(results[[nm]],file.path(out,paste0(nm,'_nominal_and_BH_tests.csv')),row.names=FALSE)
 cat('\n',nm,'nominal P<.05:',sum(results[[nm]]$P<.05),'BH FDR<.05:',sum(results[[nm]]$BH_FDR<.05),'\n')
 print(results[[nm]][results[[nm]]$P<.05,],row.names=FALSE)
}

# Reproduce Figure6A and adjust its rank correlation for purity and metastasis.
risk<-read.delim(file.path(src,'40.immuneCor/riskTrain.txt'),row.names=1,check.names=FALSE)
rownames(risk)<-canon(rownames(risk))
ci<-read.delim(file.path(src,'39.CIBERSORT/CIBERSORT-Results.txt'),row.names=1,check.names=FALSE)
rownames(ci)<-canon(rownames(ci));ci<-ci[ci[['P-value']]<0.05,,drop=FALSE]
ix<-intersect(rownames(ci),rownames(risk));ci<-ci[ix,,drop=FALSE];risk<-risk[ix,,drop=FALSE]
clinical<-read.delim(file.path(src,'35.Nomo/clinical.txt'),row.names=1,check.names=FALSE)
rownames(clinical)<-canon(rownames(clinical))
x<-risk$riskScore;x<-pmin(x,quantile(x,0.99))
y<-ci[['T cells CD4 memory activated']]
cell_cols<-setdiff(names(ci),c('P-value','Correlation','RMSE'))
cors<-do.call(rbind,lapply(cell_cols,function(nm){
 z<-ci[[nm]]
 if(sd(z)==0)return(data.frame(cell_type=nm,SD=0,rho=NA_real_,P=NA_real_))
 cc<-cor.test(x,z,method='spearman',exact=FALSE)
 data.frame(cell_type=nm,SD=sd(z),rho=unname(cc$estimate),P=cc$p.value)
}))
cors$BH_FDR_all22<-p.adjust(cors$P,'BH')
eligible<-cors$SD>=0.01
cors$BH_FDR_original_SDfilter<-NA_real_
cors$BH_FDR_original_SDfilter[eligible]<-p.adjust(cors$P[eligible],'BH')
write.csv(cors,file.path(out,'Figure6A_CIBERSORT_correlations.csv'),row.names=FALSE)
z<-data.frame(purity=est[ix,'TumorPurity'],metastasis=as.integer(clinical[ix,'Metastasis']=='Metastatic'))
partial<-function(x,y,z){
 xx<-rank(x);yy<-rank(y)
 design<-cbind(1,apply(as.matrix(z),2,rank))
 if(ncol(as.matrix(z))==0)design<-matrix(1,length(x),1)
 ex<-lm.fit(design,xx)$residuals;ey<-lm.fit(design,yy)$residuals
 r<-cor(ex,ey)
 df<-length(x)-ncol(design)-1
 p<-2*pt(-abs(r*sqrt(df/(1-r*r))),df)
 c(rho=r,P=p)
}
partial_all<-t(vapply(cell_cols,function(nm){
 if(sd(ci[[nm]])==0)return(c(rho=NA_real_,P=NA_real_))
 partial(x,ci[[nm]],z)
},numeric(2)))
cors$purity_metastasis_adjusted_rho<-partial_all[,'rho']
cors$purity_metastasis_adjusted_P<-partial_all[,'P']
cors$purity_metastasis_adjusted_BH_FDR<-p.adjust(cors$purity_metastasis_adjusted_P,'BH')
write.csv(cors,file.path(out,'Figure6A_CIBERSORT_correlations.csv'),row.names=FALSE)
cor_rows<-list()
for(nm in c('Unadjusted','Purity','Metastasis','Purity_and_metastasis')){
 zz<-switch(nm,Unadjusted=z[,FALSE,drop=FALSE],Purity=z['purity'],Metastasis=z['metastasis'],Purity_and_metastasis=z)
 if(nm=='Unadjusted'){
  cc<-cor.test(x,y,method='spearman',exact=FALSE);point<-c(rho=unname(cc$estimate),P=cc$p.value)
 }else point<-partial(x,y,zz)
 set.seed(20261113)
 bs<-replicate(2000,{j<-sample.int(length(x),length(x),replace=TRUE)
  if(nm=='Unadjusted')cor(x[j],y[j],method='spearman') else partial(x[j],y[j],zz[j,,drop=FALSE])['rho']})
 cor_rows[[nm]]<-data.frame(model=nm,n=length(x),rho=point['rho'],nominal_P=point['P'],
   lower95=unname(quantile(bs,.025,na.rm=TRUE)),upper95=unname(quantile(bs,.975,na.rm=TRUE)))
}
cor_rows<-do.call(rbind,cor_rows)
write.csv(cor_rows,file.path(out,'Figure6A_purity_metastasis_adjustment.csv'),row.names=FALSE)
write.csv(data.frame(sample=ix,risk_original=risk$riskScore,risk_capped99=x,activated_CD4_memory_fraction=y,z),
 file.path(out,'Figure6A_analysis_data.csv'),row.names=FALSE)
cat('\nFigure6A cell result:\n');print(cors[cors$cell_type=='T cells CD4 memory activated',])
cat('\nFigure6A confounding sensitivity:\n');print(cor_rows)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('\nPanelwise BH recomputation completed. Original figure annotations were nominal, not BH-adjusted.\n')
