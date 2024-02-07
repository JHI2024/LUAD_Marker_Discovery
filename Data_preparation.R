setwd("data")
library(TCGAbiolinks)

# 1. TCGA clinical data
target_cancer <- "TCGA-LUAD"
query<-GDCquery(
  project = target_cancer,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  file.type = ".txt",
)

GDCdownload(
  query = query
)

clinical<-GDCprepare(query)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})
