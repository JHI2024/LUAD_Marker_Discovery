source('[JHI] luad_nte_and_not_nte_clinical_data.R')
source('[JHI]download_luad_samples_data.R')
#-------
recurrence_patient<-subset(clinical_patient,clinical_patient$bcr_patient_barcode%in%unique(followupdata_of_nte$bcr_patient_barcode)==TRUE
)

#RFS and vital status of nte patients
nte_patient_rfs<-subset(followupdata_of_nte,duplicated(followupdata_of_nte$bcr_patient_barcode)==FALSE)
nte_patient_rfs<-data.frame(patient=nte_patient_rfs$bcr_patient_barcode,vital_status=nte_patient_rfs$vital_status,RFS=nte_patient_rfs$new_tumor_event_dx_days_to)
nte_patient_rfs_sample_list<-whole_tp_samples_list[whole_tp_samples_list%in%nte_patient_rfs$patient==TRUE]
nte_patient_rfs<-subset(nte_patient_rfs,nte_patient_rfs$patient%in%nte_patient_rfs_sample_list==TRUE & as.numeric(nte_patient_rfs$RFS)<=1825)

rfs_num<-as.numeric(nte_patient_rfs$RFS)

#Get gene expression data of these matched nte patients(The barcode of patients should correspond with sample's code)
nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%nte_patient_rfs$patient==TRUE]
patient_order<-sort(colnames(nte_patient_samples))
patient_order_num=vector()
for (i in 1:length(colnames(nte_patient_samples))){
  patient_order_num[i]<-which(colnames(nte_patient_samples)%in%patient_order[i]==TRUE)
};nte_patient_samples<-nte_patient_samples[,patient_order_num]
#sum(substr(colnames(nte_patient_samples),1,12)==nte_patient_rfs$patient)==length(colnames(nte_patient_samples))


##Calculate the lm p-value before further processing
lm_pval<-apply(nte_patient_samples,1,function(x){
  coefficients(summary(lm(x~rfs_num)))[,4]
  }
);lm_pval<-t(lm_pval);

lm_pval<-lm_pval[,2];lm_pval<-sort(lm_pval)

for (i in 1:length(lm_pval)){
  if (lm_pval[i]>=0.01){
    stop_point<<-i
    break
  }
}
rfs_marker_list<-sort(names(lm_pval[1:stop_point]))

rfs_marker_pos<-which(rownames(nte_patient_samples)%in%rfs_marker_list==TRUE)

pdf('rfs marker exp plot.pdf',width=12,height=10)
for (i in 1:length(rfs_marker_pos)){
  pos<-rfs_marker_pos[i]
  plot(nte_patient_samples[pos,]~rfs_num)
}
dev.off()
