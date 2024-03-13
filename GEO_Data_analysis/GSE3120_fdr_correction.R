gse31210_cox_results<-read.csv('gse31210 uni_cox_analysis_results.csv')
sort_pval<-sort(gse31210_cox_results$p.value)
index<-vector()
for (i in 1:nrow(gse31210_cox_results)){
  index[i]=which(gse31210_cox_results$p.value%in%sort_pval[i]==TRUE)
}

gse31210_cox_results<-gse31210_cox_results[index,]

for (i in 1:nrow(gse31210_cox_results)){
  if (gse31210_cox_results$p.value[i]>=(0.001*i/22600)){
    stop_point<<-i
    break
  }
}

gse31210_cox_results
stop_point

gse31210_cox_results$p.value
sort(rfs_marker_tcga)
sort(rfs_marker_gse37745)
