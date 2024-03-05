require(TCGAbiolinks)

# 1. TCGA clinical data
target_cancer <- "TCGA-LUAD"
query<-GDCquery(
  project = target_cancer,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = 'BCR Biotab'
)

GDCdownload(
  query = query
)

clinical<-GDCprepare(query)


# Removing irrelevant rows from the clinical information object (rows containing column names and column dictionary indices)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})


##clinical data category
clinical_patient<-clinical$clinical_patient_luad
clinical_radiation<-clinical$clinical_radiation_luad
clinical_drug<-clinical$clinical_drug_luad
clinical_nte<-clinical$clinical_nte_luad
clinical_followup<-clinical$clinical_follow_up_v1.0_luad


##nte and not nte patient data from follow up data
followupdata_of_nte<-subset(clinical_followup,clinical_followup$new_tumor_event_dx_days_to!='[Not Available]'&clinical_followup$new_tumor_event_dx_days_to!='[Not Applicable]')
followupdata_of_nte<-subset(followupdata_of_nte,duplicated(followupdata_of_nte$bcr_patient_barcode)==FALSE)

followupdata_of_clinical_not_nte<-subset(clinical_followup,clinical_followup$new_tumor_event_dx_days_to=='[Not Applicable]')
followupdata_of_clinical_not_nte<-subset(followupdata_of_clinical_not_nte,duplicated(followupdata_of_clinical_not_nte$bcr_patient_barcode)==FALSE)

##Raw patient clinical data(특별한 통제 변인 설정하지 않은 두 환자 집단)
  #raw_patientdata_of_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_nte$bcr_patient_barcode==TRUE)

  #raw_patientdata_of_not_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_clinical_not_nte$bcr_patient_barcode==TRUE)


#Condition1: Residual tumor=R0이고 neoadjuvant treatmenT=NO에 해당되는 두 환자 집단  
#con1_patientdata_of_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_nte$bcr_patient_barcode==TRUE 
                           #& clinical_patient$history_neoadjuvant_treatment=='No'
                           #& clinical_patient$residual_tumor=='R0'
                            #)
# con1_patientdata_of_not_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_clinical_not_nte$bcr_patient_barcode==TRUE
#                                & clinical_patient$history_neoadjuvant_treatment=='No'
#                                & clinical_patient$residual_tumor=='R0'
# )


#Condition2: Residual tumor=R0이고 neoadjuvant treatmenT=NO & stage1에 해당되는 두 환자 집단  
# con2_patientdata_of_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_nte$bcr_patient_barcode==TRUE 
#                                 & clinical_patient$history_neoadjuvant_treatment=='No'
#                                 & clinical_patient$residual_tumor=='R0'
#                                 & (clinical_patient$ajcc_pathologic_tumor_stage=='Stage IA'|clinical_patient$ajcc_pathologic_tumor_stage=='Stage IB')
# )
# con2_patientdata_of_not_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_clinical_not_nte$bcr_patient_barcode==TRUE
#                                     & clinical_patient$history_neoadjuvant_treatment=='No'
#                                     & clinical_patient$residual_tumor=='R0'
#                                     & (clinical_patient$ajcc_pathologic_tumor_stage=='Stage IA'|clinical_patient$ajcc_pathologic_tumor_stage=='Stage IB')
# )





