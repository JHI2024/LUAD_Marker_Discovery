# 11. GSE37745
library(pbapply)
library(dplyr)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(GEOquery)
source('[JHI] univariate_cox_hazard analysis_GEO_ver.R')
####
####
#Prepare RFS data.frame of GSE37745
#GEO version 차이로 인한 문제
data <- getGEO("GSE37745")
clindata <- data$GSE37745_series_matrix.txt.gz@phenoData@data
  #writexl::write_xlsx(clindata,'raw_clinical_GSE37745_data.xlsx')
clinical_GSE37745<-data.frame(Accession=clindata$geo_accession,
                              histology=sapply(clindata$`histology:ch1`, FUN=function(x){
                                if(x=="adeno"){
                                  "ADC"
                                }else if(x=="squamous"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$`gender:ch1`, FUN=function(x){
                                if(x=="female"){
                                  "female"
                                }else if(x=="male"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`tumor stage:ch1`, FUN=function(x){
                                if(x=="1a"){
                                  "Stage IA"
                                }else if(x=="1b"){
                                  "Stage IB"
                                }else if(x=="2a"){
                                  "Stage IIA"
                                }else if(x=="2b"){
                                  "Stage IIB"
                                }else if(x=="3a"){
                                  "Stage IIIA"
                                }else if(x=="3b"){
                                  "Stage IIIB"
                                }else if(x=="4"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              adjuvant_chemo=sapply(clindata$`adjuvant treatment:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="no"){
                                  "no"
                                }else if(x=="yes"){
                                  "yes"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`days to determined death status:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`dead:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="no"){
                                  "Alive"
                                }else if(x=="yes"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`recurrence:ch1`, FUN=function(x){
                                if(x=="no"){
                                  1
                                }else if(x=="yes"){
                                  2
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              disease_free_survival=sapply(clindata$`days to recurrence / to last visit:ch1`, FUN=function(x){
                                as.numeric(x)
                              })
)
  #writexl::write_xlsx(clinical_GSE37745,'clinical_GSE37745_data.Xlsx')

data.MA=data$GSE37745_series_matrix.txt.gz@assayData$exprs
data.Platform<-data$GSE37745_series_matrix.txt.gz@featureData@data

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

expression_GSE37745<-data.MA
adc.ind<-which(clinical_GSE37745$histology=="ADC")
clinical_GSE37745<-clinical_GSE37745[adc.ind,]
expression_GSE37745<-t(expression_GSE37745[, adc.ind])


####
#Make the data.frame for RFS analysis
con_clinical_GSE37745<-subset(clinical_GSE37745,clinical_GSE37745$relapse!="[Not Available]"
                              & clinical_GSE37745$disease_free_survival !='NA'
                              & (clinical_GSE37745$ajcc_pathologic_stage=="Stage IA" |clinical_GSE37745$ajcc_pathologic_stage=="Stage IB" ) )

con_expression_GSE37745<-subset(expression_GSE37745,rownames(expression_GSE37745)%in%con_clinical_GSE37745$Accession==TRUE)
colnames(con_expression_GSE37745)<-apply(as.matrix(colnames(con_expression_GSE37745)),1,function(x){paste0('G',x)})

rfs_data_GSE37745<-cbind(data.frame(time=con_clinical_GSE37745[,9],status=as.numeric(con_clinical_GSE37745[,8])),
                         con_expression_GSE37745)
  #write.csv(rfs_data_GSE37745,'rfs_data_GSE37745.csv')
####


####
#Implement Cox analysis for RFS data of GSE37745
result_cox_gse37745<-univ_cox_analysis(rfs_data_GSE37745)
  #write.csv(result_cox_gse37745,'gse37745 uni_cox_analysis_results.csv')


