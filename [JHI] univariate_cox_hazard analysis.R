source('[JHI] prepare univariate_cox_hazard analysis.R')
library("survival");library("survminer")
####
## Perform the univariate cox regression for total 27471 genes to select the significant harzarous genes on RFS.
univ_cox_analysis<-function(rfs_data){
  covariates<-colnames(rfs_data)[c(-1:-3)]
  univ_cox_formulas<-sapply(covariates,function(x)as.formula(paste('Surv(time,status)~',x)))
  univ_cox_models<-lapply(univ_cox_formulas,function(x){
    coxph(x,data=rfs_data)
  })
  univ_cox_results<-lapply(univ_cox_models,
                           function(x){ 
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2);#coeficient beta
                             HR <-signif(x$coef[2], digits=2);#exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             HR <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                             res<-c(beta, HR, wald.test, p.value)
                             names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                           "p.value")
                             return(res)
                             #return(exp(cbind(coef(x),confint(x))))
                           })
  res<-as.data.frame(t(as.data.frame(univ_cox_results)))
  return(res)
}
###

#Ex)
  con2_rfs_data<-generate_rfs_dataframe(con2_patientdata_of_nte,con2_patientdata_of_not_nte)
  univ_cox_analysis_result<-univ_cox_analysis(con2_rfs_data)
  #write.csv(res,'uni_cox_analysis_results.csv')

