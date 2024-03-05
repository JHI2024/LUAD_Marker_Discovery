source('[JHI] luad_nte_and_not_nte_clinical_data.R')
source('[JHI]download_luad_samples_data.R')
#-------
generate_rfs_dataframe<-function(patient_data_of_nte,patient_data_of_not_nte){
  
  
  # """
  # This function generates special form of data.frame that will be used for analyzing rfs of input-patients.
  #
  # Inputs
  #       patient_data_of_nte: clinical data of raw/controlled nte patients
  #       patient_data_of_not_nte: clinical data of raw/controlled nte patients
  # Return
  #       Data.frame that contains patient id, time for recurrence, censored state, and lots of gene exp (predictors). It can be analyzed by cox-hazard method.
  # """
  
  
  ##Extract the samples for recurrence and non-recurrence patients(who have whether last contact days or death days)
  con_followupdata_of_nte<-subset(followupdata_of_nte,followupdata_of_nte$bcr_patient_barcode%in%patient_data_of_nte$bcr_patient_barcode==TRUE)
  
  con_followupdata_of_not_nte<-subset(followupdata_of_clinical_not_nte,followupdata_of_clinical_not_nte$bcr_patient_barcode%in%patient_data_of_not_nte$bcr_patient_barcode==TRUE)
  con_followupdata_of_not_nte<-subset(con_followupdata_of_not_nte,con_followupdata_of_not_nte$last_contact_days_to!='[Not Available]' 
                                        |con_followupdata_of_not_nte$death_days_to!='[Not Applicable]')
  
  con_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con_followupdata_of_nte$bcr_patient_barcode==TRUE]
  con_not_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%con_followupdata_of_not_nte$bcr_patient_barcode==TRUE]
  
  #Sorting samples code with patient code
  nte_patient_order<-sort(colnames(con_nte_patient_samples))
  nte_patient_order_num=vector()
  for (i in 1:length(colnames(con_nte_patient_samples))){
    nte_patient_order_num[i]<-which(colnames(con_nte_patient_samples)%in%nte_patient_order[i]==TRUE)
  };con_nte_patient_samples<-con_nte_patient_samples[,nte_patient_order_num]
  not_nte_patient_order<-sort(colnames(con_not_nte_patient_samples))
  not_nte_patient_order_num=vector()
  for (i in 1:length(colnames(con_not_nte_patient_samples))){
    not_nte_patient_order_num[i]<-which(colnames(con_not_nte_patient_samples)%in%not_nte_patient_order[i]==TRUE)
  };con_not_nte_patient_samples<-con_not_nte_patient_samples[,not_nte_patient_order_num]
  
  #sample이 존재하고, nte patient만 extract
  con_followupdata_of_nte<-subset(con_followupdata_of_nte,
                                  con_followupdata_of_nte$bcr_patient_barcode%in%substr(colnames(con_nte_patient_samples),1,12)==TRUE)
    sum(substr(colnames(con_nte_patient_samples),1,12)==con_followupdata_of_nte$bcr_patient_barcode)==length(colnames(con_nte_patient_samples))
  con_followupdata_of_not_nte<-subset(con_followupdata_of_not_nte,con_followupdata_of_not_nte$bcr_patient_barcode%in%substr(colnames(con_not_nte_patient_samples),1,12)==TRUE)
    sum(substr(colnames(con_not_nte_patient_samples),1,12)==con_followupdata_of_not_nte$bcr_patient_barcode)==length(colnames(con_not_nte_patient_samples))
  
    
  ##Make the data.frame of rfs-analysis
    #1. for the nte patient
  con_nte_patient_samples<-t(con_nte_patient_samples)
  nte_rfs_data<-data.frame(patient=con_followupdata_of_nte$bcr_patient_barcode,
                           time=as.numeric(con_followupdata_of_nte$new_tumor_event_dx_days_to),
                           status=rep(2,length(con_followupdata_of_nte$bcr_patient_barcode)))
  nte_rfs_data<-cbind(nte_rfs_data,con_nte_patient_samples)
    #2. for the not nte patient
  centored_time<-vector()
  last_contact<-con_followupdata_of_not_nte$last_contact_days_to;death_date<-con_followupdata_of_not_nte$death_days_to
  for (i in 1:length(con_followupdata_of_not_nte$bcr_patient_barcode)){
    if (last_contact[i]!='[Not Available]'){
      centored_time[i]<-last_contact[i]
    } 
    else{
      centored_time[i]<-death_date[i]
    }
  }
  
  con_not_nte_patient_samples<-t(con_not_nte_patient_samples)
  not_nte_rfs_data<-data.frame(patient=con_followupdata_of_not_nte$bcr_patient_barcode,
                           time=as.numeric(centored_time),
                           status=rep(1,length(con_followupdata_of_not_nte$bcr_patient_barcode)))
  not_nte_rfs_data<-cbind(not_nte_rfs_data,con_not_nte_patient_samples)
    #3. combine each rfs data
  rfs_data<-rbind(nte_rfs_data,not_nte_rfs_data)
  return (rfs_data)
}

#Ex)
con2_patientdata_of_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_nte$bcr_patient_barcode==TRUE 
                                 & clinical_patient$history_neoadjuvant_treatment=='No'
                                 & clinical_patient$residual_tumor=='R0'
                                 & (clinical_patient$ajcc_pathologic_tumor_stage=='Stage IA'|clinical_patient$ajcc_pathologic_tumor_stage=='Stage IB')
 )
 con2_patientdata_of_not_nte<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%followupdata_of_clinical_not_nte$bcr_patient_barcode==TRUE
                                     & clinical_patient$history_neoadjuvant_treatment=='No'
                                     & clinical_patient$residual_tumor=='R0'
                                     & (clinical_patient$ajcc_pathologic_tumor_stage=='Stage IA'|clinical_patient$ajcc_pathologic_tumor_stage=='Stage IB')
 )

con2_rfs_data<-generate_rfs_dataframe(con2_patientdata_of_nte,con2_patientdata_of_not_nte)

