library(TCGAbiolinks);library(dplyr);library(DT);library(SummarizedExperiment);library(EDASeq);library(edgeR)

## Download whole LUAD gene exp samples
query<-GDCquery(project = 'TCGA-LUAD',
                data.category ="Transcriptome Profiling",
                data.type="Gene Expression Quantification")
GDCdownload(query)
gene_exp_dataset<-GDCprepare(query)
gene_exp<-assay(gene_exp_dataset,'unstranded')
gene_exp<-TCGAanalyze_Normalization(tabDF = gene_exp,geneInfo = geneInfoHT)
gene_exp<- TCGAanalyze_Filtering(
  tabDF = gene_exp,
  method = "quantile", 
  qnt.cut =  0.25
)

#total sample 600개 중, tumor positive 경우만 select=>516개
tumor_positive<-TCGAquery_SampleTypes(
  barcode = colnames(gene_exp),
  typesample = c("TP")
);tumor_positive<-subset(tumor_positive,duplicated(substr(tumor_positive,1,12))==FALSE)

whole_tp_samples_list<-substr(tumor_positive,1,12)
tumor_gene_exp<-gene_exp[,tumor_positive]
ln_gene_exp<-log(tumor_gene_exp);ln_gene_exp[ln_gene_exp==-Inf]=0;ln_gene_exp<-ln_gene_exp+1


##nte patient samples and non-nte patient samples 
nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%patientdata_of_nte$bcr_patient_barcode==TRUE]
not_nte_patient_samples<-ln_gene_exp[,whole_tp_samples_list%in%patientdata_of_not_nte$bcr_patient_barcode==TRUE]
