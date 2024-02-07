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

# Removing irrelevant rows from the clinical information object (rows containing column names and column dictionary indices)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})
