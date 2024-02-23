source('[JHI] luad_nte_and_not_nte_clinical_data.R')
source('[JHI]download_luad_samples_data.R')

###############################Finding DEGs function step####################

find_degs<-function(nte_patient_samples,not_nte_patient_samples){
  library(org.Hs.eg.db)
  nte_num<-ncol(nte_patient_samples)
  not_nte_num<-ncol(not_nte_patient_samples)
  #Shapiro test(p>0.05) to select genes that follow normal-dist
  normaldist_genes_nte<-which(
    apply(nte_patient_samples,1,function(x){
      if (sum(x)!=length(x)*x[1]){
        shapiro.test(x)$p.value
      }
      else{
        return(0)
      }
    })>0.05)
  
  normaldist_genes_not_nte<-which(
    apply(not_nte_patient_samples,1,function(x){
      if (sum(x)!=length(x)*x[1]){
        shapiro.test(x)$p.value
      }
      else{
        return(0)
      }
    })>0.05)
  
  common_normaldist_genes_index<-intersect(normaldist_genes_nte,normaldist_genes_not_nte)
  
  #Merge two group for following statistical analysis
  whole_comparing_samples<-cbind(nte_patient_samples,not_nte_patient_samples)
  test_index<-c(rep(0,nrow(nte_patient_samples)))
  test_index[common_normaldist_genes_index]<-1
  whole_comparing_samples<-as.data.frame(cbind(whole_comparing_samples,test_index))
  
  #t-test for normal-dist genes and wilcoxon test for not norml-dist genes
  stat_test_pval_fold<-t(as.matrix(apply(whole_comparing_samples,1,function(x){
    if (x[nte_num+not_nte_num+1]==1){
      #t.test(x[1:nte_num],x[(nte_num+1):(not_nte_num+nte_num)])$p.value
      c(t.test(x[1:nte_num],x[(nte_num+1):(not_nte_num+nte_num)])$p.value,median(x[(nte_num+1):(not_nte_num+nte_num)])/median(x[1:nte_num]))
    }
    else{
      #wilcox.test(x[1:nte_num],x[(nte_num+1):(not_nte_num+nte_num)])$p.value
      c(t.test(x[1:nte_num],x[(nte_num+1):(not_nte_num+nte_num)])$p.value,median(x[(nte_num+1):(not_nte_num+nte_num)])/median(x[1:nte_num]))
    }
  })))
  
  fold2_stat_test_pval<-as.matrix(sort(stat_test_pval[which(stat_test_pval_fold[,2]>2|stat_test_pval_fold[,2]<0.5),1]))
  
  #Fdr correction
  for (i in 1:length(fold2_stat_test_pval)){
    if (fold2_stat_test_pval[i]>=0.01){
      stop_point<<-i
      break
    }
  }
  DEGs_list<-keys(org.Hs.eg.db,keytype = 'ENSEMBL')[keys(org.Hs.eg.db,keytype = 'ENSEMBL')%in%rownames(fold2_stat_test_pval)[1:stop_point-1]==TRUE]
  
  return (DEGs_list)
}
#####################################################


##nte patient samples and non-nte patient samples 
con1_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con1_patientdata_of_nte$bcr_patient_barcode==TRUE]
con1_not_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con1_patientdata_of_not_nte$bcr_patient_barcode==TRUE]

con2_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con2_patientdata_of_nte$bcr_patient_barcode==TRUE]
con2_not_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con2_patientdata_of_not_nte$bcr_patient_barcode==TRUE]

a<-find_degs(con1_nte_patient_samples,con1_not_nte_patient_samples)
b<-find_degs(con2_nte_patient_samples,con2_not_nte_patient_samples)
