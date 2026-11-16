options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
src <- Sys.getenv('P3H2_PROG_INPUTS',unset='P3H2_prognosis_inputs')
out <- Sys.getenv('P3H2_PROG_OUTPUTS',unset='P3H2_prognosis_outputs')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tt <- read.delim(file.path(src,'06.unicox/time.txt'),row.names=1,check.names=FALSE)
bm <- read.delim(file.path(src,'06.unicox/BMexp.txt'),row.names=1,check.names=FALSE)
clin <- read.delim(file.path(src,'35.Nomo/clinical.txt'),row.names=1,check.names=FALSE)
tme <- read.delim(file.path(src,'43.estimate/TMEscores.txt'),row.names=1,check.names=FALSE)
single <- read.delim(file.path(src,'30.sursingle/TrainP3H2.txt'),row.names=1,check.names=FALSE)
multi <- read.delim(file.path(src,'29.multiCox/lassoSigExp.txt'),row.names=1,check.names=FALSE)
canon <- function(x)sub('^TARGET[-_]','',x)
for(nm in c('tt','clin','tme','single','multi')) {
 z<-get(nm);rownames(z)<-canon(rownames(z));assign(nm,z)
}
colnames(bm)<-canon(colnames(bm))
ids <- rownames(tt)
stopifnot(nrow(tt)==85, all(tt$fustat%in%0:1), sum(tt$fustat)==27,
 all(ids%in%colnames(bm)),all(ids%in%rownames(clin)),all(ids%in%rownames(tme)),all(ids%in%rownames(single)))
d<-data.frame(id=ids,time_years=tt$futime/365,event=tt$fustat,P3H2=as.numeric(bm['P3H2',ids]),
 Age=clin[ids,'Age'],Sex=factor(clin[ids,'Gender'],levels=c('Female','Male')),
 Metastasis=factor(clin[ids,'Metastasis'],levels=c('Non-metastatic','Metastatic')),
 tme[ids,,drop=FALSE],check.names=FALSE)
rownames(d)<-ids
d$AgeGroup<-factor(as.integer(d$Age>=15),levels=0:1,labels=c('Under15','15orOlder'))
d$P3H2_SD<-as.numeric(scale(d$P3H2))
d$P3H2_log2_SD<-as.numeric(scale(log2(d$P3H2+1)))
stopifnot(max(abs(single[ids,'riskScore']-d$P3H2))<1e-6,
 all(single[ids,'fustat']==d$event),max(abs(single[ids,'futime']-d$time_years))<1e-8)
audit <- data.frame(check=c('n','events_coded_1','censored_coded_0','P3H2_SD_TPM','single_gene_input_matches_P3H2','single_gene_event_codes_match'),
 value=c(nrow(d),sum(d$event),sum(d$event==0),sd(d$P3H2),TRUE,TRUE))
write.csv(audit,file.path(out,'coding_and_expression_audit.csv'),row.names=FALSE)

# Figure 4H uses the expression itself as a variable named riskScore; high/low labels are expression groups, not the multi-gene risk index.
cut <- surv_cutpoint(data.frame(futime=d$time_years,fustat=d$event,riskScore=d$P3H2),
                     time='futime',event='fustat',variables='riskScore')
threshold<-cut$cutpoint$cutpoint
d$ExpressionGroup<-factor(ifelse(d$P3H2>threshold,'High','Low'),levels=c('Low','High'))
km<-survfit(Surv(time_years,event)~ExpressionGroup,data=d)
lr<-survdiff(Surv(time_years,event)~ExpressionGroup,data=d)
groupfit<-coxph(Surv(time_years,event)~ExpressionGroup,data=d)
gs<-summary(groupfit)
group_summary<-data.frame(cutoff_TPM=threshold,n_low=sum(d$ExpressionGroup=='Low'),n_high=sum(d$ExpressionGroup=='High'),
 high_vs_low_HR=gs$conf.int[1,1],lower95=gs$conf.int[1,3],upper95=gs$conf.int[1,4],
 nominal_logrank_P=pchisq(lr$chisq,1,lower.tail=FALSE))
write.csv(group_summary,file.path(out,'P3H2_expression_group_audit.csv'),row.names=FALSE)

forms <- list(
 Unadjusted_TPM = Surv(time_years,event)~P3H2_SD,
 Clinical_TPM = Surv(time_years,event)~P3H2_SD+AgeGroup+Sex+Metastasis,
 Purity_only_TPM = Surv(time_years,event)~P3H2_SD+TumorPurity,
 Clinical_purity_TPM = Surv(time_years,event)~P3H2_SD+AgeGroup+Sex+Metastasis+TumorPurity,
 Clinical_immune_stromal_TPM = Surv(time_years,event)~P3H2_SD+AgeGroup+Sex+Metastasis+ImmuneScore+StromalScore,
 Unadjusted_log2 = Surv(time_years,event)~P3H2_log2_SD,
 Clinical_log2 = Surv(time_years,event)~P3H2_log2_SD+AgeGroup+Sex+Metastasis,
 Clinical_purity_log2 = Surv(time_years,event)~P3H2_log2_SD+AgeGroup+Sex+Metastasis+TumorPurity,
 Clinical_immune_stromal_log2 = Surv(time_years,event)~P3H2_log2_SD+AgeGroup+Sex+Metastasis+ImmuneScore+StromalScore)
fits<-list();res<-list();phs<-list()
for(nm in names(forms)) {
 f<-coxph(forms[[nm]],data=d,x=TRUE,y=TRUE,model=TRUE)
 fits[[nm]]<-f;s<-summary(f);j<-grep('^P3H2',rownames(s$coefficients))
 res[[nm]]<-data.frame(model=nm,scale=if(grepl('log2',nm))'per SD log2(TPM+1)' else 'per SD untransformed TPM',
   n=f$n,events=f$nevent,coefficient=s$coefficients[j,1],HR=s$conf.int[j,1],
   lower95=s$conf.int[j,3],upper95=s$conf.int[j,4],nominal_P=s$coefficients[j,5])
 pp<-cox.zph(f)
 phs[[nm]]<-data.frame(model=nm,term=rownames(pp$table),pp$table,check.names=FALSE)
}
rr<-do.call(rbind,res);write.csv(rr,file.path(out,'P3H2_adjusted_Cox_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,phs),file.path(out,'P3H2_Cox_PH_tests.csv'),row.names=FALSE)
orig_uni<-coxph(Surv(time_years,event)~P3H2,data=d)
orig_multi<-coxph(Surv(futime,fustat)~ACTA2+FRAS1+P3H2,data=multi)
cat('Original-scale univariate P3H2:\n');print(summary(orig_uni)$coefficients)
cat('Original three-gene model P3H2:\n');print(summary(orig_multi)$coefficients['P3H2',])

correlate<-function(mat,category) {
 rows<-lapply(colnames(mat),function(nm){
  z<-mat[,nm]
  if(!is.numeric(z)||sd(z,na.rm=TRUE)==0)return(NULL)
  co<-cor.test(d$P3H2,z,method='spearman',exact=FALSE)
  data.frame(category=category,feature=nm,n=sum(is.finite(z)),rho=unname(co$estimate),P=co$p.value)
 })
 z<-do.call(rbind,rows);z$BH_FDR<-p.adjust(z$P,'BH');z
}
cor_est<-correlate(as.matrix(d[,c('TumorPurity','ImmuneScore','StromalScore','ESTIMATEScore')]),'ESTIMATE')
ci<-read.delim(file.path(src,'39.CIBERSORT/CIBERSORT-Results.txt'),row.names=1,check.names=FALSE)
rownames(ci)<-canon(rownames(ci))
stopifnot(all(ids%in%rownames(ci)))
ci<-ci[ids,,drop=FALSE]
cell_cols<-setdiff(colnames(ci),c('P-value','Correlation','RMSE'))
cor_ci<-correlate(as.matrix(ci[,cell_cols]),'CIBERSORT_all_samples')
xc<-read.delim(file.path(src,'42 Xcell/xCell_mRNA_tpm_xCell_0400081323.txt'),row.names=1,check.names=FALSE)
if(all(ids%in%canon(colnames(xc)))) {
 colnames(xc)<-canon(colnames(xc));xc<-t(as.matrix(xc[,ids,drop=FALSE]))
 cor_x<-correlate(xc,'xCell')
} else stop('xCell full-score input orientation or sample identifiers did not match')
write.csv(cor_est,file.path(out,'P3H2_ESTIMATE_correlations.csv'),row.names=FALSE)
write.csv(cor_ci[order(cor_ci$BH_FDR),],file.path(out,'P3H2_CIBERSORT_correlations.csv'),row.names=FALSE)
write.csv(cor_x[order(cor_x$BH_FDR),],file.path(out,'P3H2_xCell_correlations.csv'),row.names=FALSE)
keep<-ci[['P-value']]<0.05
valid_ci<-do.call(rbind,lapply(cell_cols,function(nm) {
 co<-cor.test(d$P3H2[keep],ci[keep,nm],method='spearman',exact=FALSE)
 data.frame(feature=nm,n=sum(keep),rho=unname(co$estimate),P=co$p.value)
}))
valid_ci$BH_FDR<-p.adjust(valid_ci$P,'BH')
write.csv(valid_ci[order(valid_ci$BH_FDR),],file.path(out,'P3H2_CIBERSORT_Pbelow005_correlations.csv'),row.names=FALSE)
write.csv(d,file.path(out,'P3H2_cohort_analysis_data.csv'),row.names=FALSE)
saveRDS(list(data=d,fits=fits,km=km,group_audit=group_summary),file.path(out,'P3H2_prognosis_audit.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('\nExpression groups:\n');print(group_summary)
cat('\nAdjusted P3H2 results (exploratory post-selection sensitivity):\n');print(rr)
cat('\nESTIMATE correlations:\n');print(cor_est)
cat('\nTop xCell correlations:\n');print(head(cor_x[order(cor_x$BH_FDR),],10))
cat('\nOriginal Xena-derived endpoint file retained: 85 patients, 27 deaths; 1=death, 0=alive at last follow-up (right-censored). No GDC endpoint substitution.\n')
