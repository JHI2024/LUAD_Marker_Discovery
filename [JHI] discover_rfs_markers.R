#Download result of cox analysis
library(readxl)
univ_cox_analysis_results<-read_xlsx('uni_cox_analysis_results.xlsx')

#Select the significantly expressed genes for rfs with p-value
library(org.Hs.eg.db)
rfs_marker<-subset(univ_cox_analysis_results,univ_cox_analysis_results$p.value<0.001)
DEGs_list<-keys(org.Hs.eg.db,keytype = 'ENSEMBL')[keys(org.Hs.eg.db,keytype = 'ENSEMBL')%in%rfs_marker$...1==TRUE]
deg_index<-which(keys(org.Hs.eg.db,keytype = 'ENSEMBL')%in%DEGs_list==TRUE)
DEGs_info<-data.frame(ENSEMBL=DEGs_list,Gene_code=keys(org.Hs.eg.db,keytype = 'SYMBOL')[deg_index])
writexl::write_xlsx(DEGs_info,'Univariate cox hazard analysis based RFS markers.xlsx')


