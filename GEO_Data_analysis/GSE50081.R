# 14. GSE50081
library(pbapply)
library(dplyr)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(GEOquery)
source('[JHI] univariate_cox_hazard analysis_GEO_ver.R')

data <- getGEO("GSE50081")
clindata <- data$GSE50081_series_matrix.txt.gz@phenoData@data
  #writexl::write_xlsx(clindata,'raw_clinical_GSE50081_data.xlsx')

clinical_GSE50081<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$`Sex:ch1`, FUN=function(x){
                                if(x=="F"){
                                  "female"
                                }else if(x=="M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$`histology:ch1`, FUN=function(x){
                                if(x=="adenocarcinoma"){
                                  "ADC"
                                }else if (x=="squamous cell carcinoma"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking:ch1`, FUN=function(x){
                                if(x=="Current"){
                                  "Ever smoker"
                                }else if(x=="Ex-smoker"){
                                  "Ever smoker"
                                }else if(x=="Never"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`Stage:ch1`, FUN=function(x){
                                if(x == "1A"){
                                  "Stage IA"
                                }else if(x == "1B"){
                                  "Stage IB"
                                }else if(x == "2A"){
                                  "Stage IIA"
                                }else if(x=="2B"){
                                  "Stage IIB"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(clindata$`t-stage:ch1`, FUN=function(x){
                                if(x=="1"){
                                  "T1"
                                }else if(x=="2"){
                                  "T2"
                                }else if(x=="3"){
                                  "T3"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(clindata$`n-stage:ch1`, FUN=function(x){
                                if(x=="0"){
                                  "N0"
                                }else if(x=="1"){
                                  "N1"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`survival time:ch1`, FUN=function(x){
                                as.numeric(x)*12
                              }),
                              vital_status=sapply(clindata$`status:ch1` , FUN=function(x){
                                if(x=="alive"){
                                  "Alive"
                                }else if(x=="dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              disease_free_survival=sapply(clindata$`disease-free survival time:ch1`, FUN=function(x){
                                round(as.numeric(x)*365)
                              }),
                              relapse=sapply(clindata$`recurrence:ch1`, FUN=function(x){
                                if(x=="N"){
                                  1
                                }else if(x=="Y"){
                                  2
                                }else
                                  "[Not Available]"
                              })
)
  #writexl::write_xlsx(clinical_GSE50081,'clinical_GSE50081_data.Xlsx')

data.MA=data$GSE50081_series_matrix.txt.gz@assayData$exprs
data.Platform<-data$GSE50081_series_matrix.txt.gz@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

# # MAD calculation
# if(mad.exclusion!=0){
#   mad.exp<-apply(data.MA, MARGIN=1, mad)
#   data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
# }
# 
# if(qnorm){
#   data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
# }

expression_GSE50081<-t(data.MA)
adc.ind<-which(clinical_GSE50081$histology=="ADC")
clinical_GSE50081<-clinical_GSE50081[adc.ind,]
expression_GSE50081<-expression_GSE50081[adc.ind,]
####


####
#Make the data.frame for RFS analysis
con_clinical_GSE50081<-subset(clinical_GSE50081,
                              clinical_GSE50081$relapse !='[Not Available]'
                              & (clinical_GSE50081$ajcc_pathologic_stage=="Stage IA" |clinical_GSE50081$ajcc_pathologic_stage=="Stage IB" ))

con_expression_GSE50081<-subset(expression_GSE50081,rownames(expression_GSE50081)%in%con_clinical_GSE50081$Accession==TRUE)
colnames(con_expression_GSE50081)<-apply(as.matrix(colnames(con_expression_GSE50081)),1,function(x){paste0('G',x)})


rfs_data_GSE50081<-cbind(data.frame(time=as.numeric(con_clinical_GSE50081[,10]),status=as.numeric(con_clinical_GSE50081[,11])),
                         con_expression_GSE50081)
  #write.csv(rfs_data_GSE50081,'rfs_data_GSE50081.csv')
####


####
#Implement Cox analysis for RFS data of GSE50081
result_cox_gse50081<-univ_cox_analysis(rfs_data_GSE50081)
  write.csv(result_cox_gse50081,'gse50081 uni_cox_analysis_results.csv')
sort(as.numeric(result_cox_gse50081$p.value))
