setwd("data")
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
  #writexl::write_xlsx(clinical$clinical_patient_luad,'patient data.xlsx')
  #writexl::write_xlsx(clinical$clinical_radiation_luad,'radiation data.xlsx')
  #writexl::write_xlsx(clinical$clinical_drug_luad,'drug data.xlsx')
  #writexl::write_xlsx(clinical$clinical_nte_luad,'nte data.xlsx')
  #writexl::write_xlsx(clinical$clinical_follow_up_v1.0_luad,'follow up data.xlsx')

##1. nte patient info list 생성
clinical_ntedata_of_followup<-subset(clinical_nte,clinical_nte$bcr_patient_barcode%in%clinical_followup$bcr_patient_barcode==TRUE)
followupdata_of_clinical_nte<-subset(clinical_followup,clinical_followup$bcr_patient_barcode%in%clinical_nte$bcr_patient_barcode==TRUE)
followupdata_of_nte<-subset(clinical_followup,clinical_followup$new_tumor_event_dx_days_to!='[Not Available]'&clinical_followup$new_tumor_event_dx_days_to!='[Not Applicable]')
followupdata_of_clinical_not_nte<-subset(clinical_followup,clinical_followup$new_tumor_event_dx_days_to=='[Not Applicable]')
patientdata_of_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%unique(followupdata_of_nte$bcr_patient_barcode)==TRUE)
patientdata_of_not_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%unique(followupdata_of_clinical_not_nte$bcr_patient_barcode)==TRUE)

luad_nte_info_list<-list(clinical_ntedata=clinical_nte,
                            clinical_ntedata_of_followup=clinical_ntedata_of_followup,
                            followupdata_of_clinical_nte=followupdata_of_clinical_nte,
                            followupdata_of_nte=followupdata_of_nte,
                            patientdata_of_nte=patientdata_of_nte
                            )
  writexl::write_xlsx(luad_nte_info_list,'luad_nte_info_list.xlsx')



  