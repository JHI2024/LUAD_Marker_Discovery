library(pbapply)
library(dplyr)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(GEOquery)
source('[JHI] univariate_cox_hazard analysis_GEO_ver.R')
####



####
#Prepare RFS data.frame of GSE41271

#GEO version 차이로 인한 문제
data <- getGEO('GSE41271')

#GEO version 차이로 인한 문제
clindata <- data$GSE41271_series_matrix.txt.gz@phenoData@data
  #writexl::write_xlsx(clindata,'raw_clinical_GSE41271_data.xlsx')
clinical_GSE41271<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$characteristics_ch1.1, FUN=function(x){
                                if(x=="gender: F"){
                                  "female"
                                }else if(x=="gender: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              race=sapply(clindata$`race:ch1`, FUN=function(x){
                                if(x=="African American"){
                                  "black or african american"
                                }else if (x=="Asian"){
                                  "asian"
                                }else if(x=="Hispanic"){
                                  "hispanic"
                                }else if(x=="Caucasian"){
                                  "white"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$characteristics_ch1, FUN=function(x){
                                if(x=="histology: Adenocarcinoma"){
                                  "ADC"
                                }else if(x=="histology: Squamous"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`tobacco history:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="Y"){
                                  return("Ever smoker")
                                }else if(x=="N"){
                                  return("Never smoker")
                                }else{
                                  return("[Not Available]")
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x == "final patient stage: IA"){
                                  "Stage IA"
                                }else if(x == "final patient stage: IB"){
                                  "Stage IB"
                                }else if(x == "final patient stage: IIA"){
                                  "Stage IIA"
                                }else if(x == "final patient stage: IIB"){
                                  "Stage IIB"
                                }else if(x == "final patient stage: IIIA"){
                                  "Stage IIIA"
                                }else if(x == "final patient stage: IIIB"){
                                  "Stage IIIB"
                                }else if(x == "final patient stage: IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                surgery.date<-as_date(clindata$`date of surgery:ch1`[k])
                                last.followup.date<-as_date(clindata$`last follow-up survival:ch1`)[k]
                                survival.days<-as.numeric(last.followup.date-surgery.date)
                                if(is.na(survival.days)){
                                  return("[Not Available]")
                                }else{
                                  return(survival.days)  
                                }
                              }),
                              vital_status=sapply(clindata$characteristics_ch1.7 , FUN=function(x){
                                if(x=="vital statistics: A"){
                                  "Alive"
                                }else if(x=="vital statistics: D"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              state=sapply(clindata$characteristics_ch1.9, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="recurrence: N"){
                                  1
                                }else if(x=="recurrence: Y"){
                                  2
                                }else
                                  "[Not Available]"
                              }),
                              time=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                surgery.date<-as_date(clindata$`date of surgery:ch1`[k])
                                last.followup.date<-as_date(clindata$`last follow-up recurrence:ch1`)[k]
                                survival.days<-as.numeric(last.followup.date-surgery.date)
                                if(is.na(survival.days)){
                                  return("[Not Available]")
                                }else{
                                  return(survival.days)  
                                }
                              })
)
  #writexl::write_xlsx(clinical_GSE41271,'clinical_GSE41271_data.Xlsx')

data.MA=data$GSE41271_series_matrix.txt.gz@assayData$exprs
data.Platform<-data$GSE41271_series_matrix.txt.gz@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
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
# if(qnorm){
#   data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
# }

expression_GSE41271<-t(data.MA)
adc.ind<-which(clinical_GSE41271$histology=="ADC")
clinical_GSE41271<-clinical_GSE41271[adc.ind,]
expression_GSE41271<-expression_GSE41271[adc.ind,]
####


####
#Make the data.frame for RFS analysis
con_clinical_GSE41271<-subset(clinical_GSE41271,
                              clinical_GSE41271$state!="[Not Available]"
                              & clinical_GSE41271$time !='[Not Available]'
                              & (clinical_GSE41271$ajcc_pathologic_stage=="Stage IA" |clinical_GSE41271$ajcc_pathologic_stage=="Stage IB" ))
con_expression_GSE41271<-subset(expression_GSE41271,rownames(expression_GSE41271)%in%con_clinical_GSE41271$Accession==TRUE)
colnames(con_expression_GSE41271)<-apply(as.matrix(colnames(con_expression_GSE41271)),1,function(x){paste0('G',x)})


rfs_data_GSE41271<-cbind(data.frame(time=as.numeric(con_clinical_GSE41271[,10]),status=as.numeric(con_clinical_GSE41271[,9])),
                              con_expression_GSE41271)
  #write.csv(rfs_data_GSE41271,'rfs_data_GSE41271.csv')
####


####
#Implement Cox analysis for RFS data of GSE41271
result_cox_gse41271<-univ_cox_analysis(rfs_data_GSE41271)
  #write.csv(result_cox_gse41271,'gse41271 uni_cox_analysis_results.csv')
  