library(pbapply)
library(dplyr)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(TCGAbiolinks)
library(GEOquery)

# Terms representing elements not available were defined in advance of the analysis
na.terms<-c(NA, "[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Not Evaluated]", "[Unknown]", "MX", "TX", "NX")

alive.days.thres=0*365
alive.month.thres=alive.days.thres/30.41667

death.days.thres=100000*365
death.month.thres=death.days.thres/30.41667

mad.exclusion<-0
qnorm=FALSE


# 1. TCGA data
target_cancer <- "TCGA-LUAD"

###### Clinical data
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

####### Mutation data
# query <- GDCquery(
#   project = "TCGA-LUAD", 
#   data.category = "Simple Nucleotide Variation", 
#   access = "open", 
#   legacy = FALSE, 
#   data.type = "Masked Somatic Mutation", 
# )
# GDCdownload(query)
# maf <- GDCprepare(query)
maf<-readRDS("tcga_luad_mutation.RDS")

######## Expression data
# query <- GDCquery(
#   project = target_cancer,
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification", 
#   workflow.type = "STAR - Counts"
# )
# GDCdownload(query)
# data <- GDCprepare(query)
data<-readRDS("tcga_luad_rnaseq.RDS")

clinical_df<-data.frame(data@colData@listData[c('barcode', 'sample_type', 'ajcc_pathologic_stage', 
                                                "days_to_last_follow_up", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", 
                                                "ajcc_staging_system_edition", 'race', 'gender', "vital_status", "days_to_death", "paper_expression_subtype")])

clinical_df$survival_time<-sapply(c(1:length(clinical_df$vital_status)), FUN=function(x){
  t=setdiff(c(clinical_df$days_to_last_follow_up[x], clinical_df$days_to_death[x]), na.terms)
  if(length(t)==0){
    return("[Not Available]")
  }
  return(max(as.numeric(t)))
})

clinical_df$relapse<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  if(is.na(match(barcode.cut, clinical[[2]]$bcr_patient_barcode))){
    return("NO")
  }else{
    return("YES")
  }
})

clinical_df$race<-sapply(clinical_df$race, FUN=function(x){
  if(x=="not reported")
    return("[Not Available]")
  else
    return(x)
})

clinical_df$adjuvant_chemo<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Pharmaceutical Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

clinical_df$adjuvant_radiation<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Radiation Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

tumor.ind<-which(clinical_df$sample_type=="Primary Tumor")

clinical_df<-clinical_df[tumor.ind,]

clinical_df$barcode.short<-sapply(clinical_df$barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

# smoking status
clinical_df$smoking_status<-sapply(clinical_df$barcode.short, FUN=function(barcode){
  ind<-match(barcode, clinical$clinical_patient_luad$bcr_patient_barcode)
  smoking.indicator<-clinical$clinical_patient_luad$tobacco_smoking_history_indicator[ind]
  if(smoking.indicator == "1"){
    res<-"Never smoker"
  }else if (smoking.indicator %in% c("2", "3", "4", "5")){
    res<-"Ever smoker"
  }else{
    res<-"[Not Available]"
  }
  return(res)
})

# Predominant subtype
clinical_df$predominant_subtype<-sapply(clinical_df$barcode.short, FUN=function(barcode){
  ind<-match(barcode, clinical$clinical_patient_luad$bcr_patient_barcode)
  x<-clinical$clinical_patient_luad$histologic_diagnosis...68[ind]
  if(x=="Lung Acinar Adenocarcinoma"){
    "Acinar"
  }else if(x=="Lung Micropapillary Adenocarcinoma"){
    "Micropapillary"
  }else if(x %in% c("Lung Bronchioloalveolar Carcinoma Mucinous", "Lung Mucinous Adenocarcinoma", "Mucinous (Colloid) Carcinoma")){
    "Mucinous"
  }else if(x=="Lung Papillary Adenocarcinoma"){
    "Papillary"
  }else{
    "[Not Available]"
  }
  
})

# EGFR mutation
mut.egfr<-maf %>% filter(Hugo_Symbol == "EGFR")
mut.egfr$barcode.short<-sapply(mut.egfr$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$egfr_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.egfr$barcode.short)
  if(is.na(ind))
    res<-"EGFR-WT"
  else
    res<-"EGFR-Mut"
})

# KRAS mutation
mut.kras<-maf %>% filter(Hugo_Symbol == "KRAS")
mut.kras$barcode.short<-sapply(mut.kras$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$kras_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.kras$barcode.short)
  if(is.na(ind))
    res<-"KRAS-WT"
  else
    res<-"KRAS-Mut"
})

# p53 mutation
mut.tp53<-maf %>% filter(Hugo_Symbol == "TP53")
mut.tp53$barcode.short<-sapply(mut.tp53$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$p53_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.tp53$barcode.short)
  if(is.na(ind))
    res<-"p53-WT"
  else
    res<-"p53-Mut"
})

# STK11 mutation
mut.stk11<-maf %>% filter(Hugo_Symbol == "STK11")
mut.stk11$barcode.short<-sapply(mut.stk11$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df$stk11_mut<-sapply(clinical_df$barcode.short, FUN=function(x){
  ind<-match(x, mut.stk11$barcode.short)
  if(is.na(ind))
    res<-"STK11-WT"
  else
    res<-"STK11-Mut"
})


data.luad <- assay(data, "fpkm_uq_unstrand")
data.luad.row.info<-rowRanges(data)
protein_coding_gene.ind<-which(data.luad.row.info@elementMetadata@listData$gene_type=="protein_coding")

data.luad.protein_coding<-data.luad[protein_coding_gene.ind, tumor.ind]
#zeros.length<-apply(data.luad.protein_coding, MARGIN=1, FUN=function(k){length(which(k==0))})
#data.luad.protein_coding<-data.luad.protein_coding[which(zeros.length<=180),]

expression_df<-t(data.luad.protein_coding)
colnames(expression_df)<-row.names(data.luad.protein_coding)
expression_df<-log(expression_df+1, base=2)

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(expression_df, MARGIN=2, mad)
  expression_df<-expression_df[,which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)])]
}
if(qnorm){
  expression_df<-t(normalize.quantiles(t(expression_df), keep.names=TRUE))
}

gene_info<-data.frame(gene_id=data.luad.row.info@elementMetadata@listData$gene_id, 
                      gene_symbol=data.luad.row.info@elementMetadata@listData$gene_name) %>% filter(gene_id %in% colnames(expression_df))

gene_info$ENTREZID<-pbsapply(gene_info$gene_symbol, FUN=function(x){
  tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                 keys=x, 
                                 columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")$ENTREZID[1], 
           error= function(e){ 
             return(NA)
           })
})


gene_info$ENTREZID<-pbsapply(c(1:length(gene_info$ENTREZID)), FUN=function(k){
  if(!is.na(gene_info$ENTREZID[k])){
    return(gene_info$ENTREZID[k])
  }else{
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                   keys=unlist(strsplit(gene_info$gene_id[k], ".", fixed=T))[1], 
                                   columns=c("ENSEMBL", "ENTREZID"), keytype = "ENSEMBL")$ENTREZID[1], 
             error= function(e){ 
               return(NA)
             })
  }
})


clinical_TCGA<-clinical_df
expression_TCGA<-expression_df

rm(clinical_df)
rm(expression_df)
rm(data.luad); rm(data.luad.protein_coding); rm(data.Platform); rm(data.luad.row.info)
rm(maf);rm(mut.egfr);rm(mut.kras);rm(mut.stk11);rm(mut.tp53);rm(query)

# 2. GSE3141
data <- getGEO(filename="GEO/GSE3141_series_matrix.txt.gz")
clindata <- data@phenoData@data

labels<-lapply(clindata$characteristics_ch1, FUN=function(x){
  strsplit(x, split = ";", fixed=T) %>% unlist
})

clinical_GSE3141<-data.frame(Accession=clindata$geo_accession,
                             histology=sapply(sapply(labels, FUN=function(x){
                               unlist(strsplit(x[1], split=": ", fixed=T))[3]
                             }), FUN=function(x){
                               if(x=="A"){
                                 "ADC"
                               }else if(x=="S"){
                                 "SCC"
                               }else{
                                 "[Not Available]"
                               }
                             }),
                             survival_time=sapply(labels, FUN=function(x){
                               as.numeric(unlist(strsplit(x[2], split=": ", fixed=T))[2])
                             }),
                             vital_status=sapply(sapply(labels, FUN=function(x){
                               unlist(strsplit(x[3], split=": ", fixed=T))[2]
                             }), FUN=function(x){
                               if(x=="0"){
                                 "Alive"
                               }else if(x=="1"){
                                 "Dead"
                               }else{
                                 "[Not Available]"
                               }
                             })
)

clinical_GSE3141$survival_time[106]<-72

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

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
data.MA<-log2(data.MA)

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE3141<-data.MA
adc.ind<-which(clinical_GSE3141$histology=="ADC")
clinical_GSE3141<-clinical_GSE3141[adc.ind,]
expression_GSE3141<-t(expression_GSE3141[, adc.ind])


# 3. GSE11969
data <- getGEO(filename="GEO/GSE11969_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE11969<-data.frame(Accession=clindata$geo_accession,
                              histology=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x=="Histology: AD"){
                                  "ADC"
                                }else if(x=="Histology: SQ"){
                                  "SCC"
                                }else if(x=="Histology: Normal lung tissue"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$characteristics_ch1.3, FUN=function(x){
                                if(x=="Sex: F"){
                                  "female"
                                }else if(x=="Sex: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(str_sub(clindata$`pTNM:ch1`,1,2), FUN=function(x){
                                if(x %in% c("T1","T2", "T3", "T4", "TX")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(str_sub(clindata$`pTNM:ch1`,3,4), FUN=function(x){
                                if(x %in% c("NX", "N0", "N1", "N2", "N3")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              # thie dataset has ajcc_pathologic_m label but only contains M0
                              ajcc_pathologic_stage=sapply(clindata$`pStage:ch1`, FUN=function(x){
                                if(x=="Stage: 1"){
                                  "Stage I"
                                }else if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`EGFR status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "EGFR-WT"
                                }else if(x=="Mut"){
                                  "EGFR-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`K-ras Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "KRAS-WT"
                                }else if (x=="Mut"){
                                  "KRAS-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`p53 Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "p53-WT"
                                }else if(x=="Mut"){
                                  "p53-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              grade=sapply(clindata$`Grade:ch1`, FUN=function(x){
                                if(x=="Moderately"){
                                  "Moderately"
                                }else if(x=="Poorly"){
                                  "Poorly"
                                }else if (x=="Well"){
                                  "Well"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`Survival (days):ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`Status:ch1`, FUN=function(x){
                                if(x=="Alive"){
                                  "Alive"
                                }else if(x=="Dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`Evidence of relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="N"){
                                  "NO"
                                }else if(x=="Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data


# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$`Gene symbol`[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[mapping.table$ID,]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}


expression_GSE11969<-t(data.MA)
adc.ind<-which(clinical_GSE11969$histology=="ADC")
clinical_GSE11969<-clinical_GSE11969[adc.ind,]
expression_GSE11969<-expression_GSE11969[adc.ind, ]


# 4. GSE13213
data <- getGEO(filename="GEO/GSE13213_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE13213<-data.frame(Accession=clindata$geo_accession,
                              histology=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x=="Histology: AD"){
                                  "ADC"
                                }else if(x=="Histology: SQ"){
                                  "SCC"
                                }else if(x=="Histology: Normal lung tissue"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$characteristics_ch1.3, FUN=function(x){
                                if(x=="Sex: F"){
                                  "female"
                                }else if(x=="Sex: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(str_sub(clindata$`TNM (Pathological):ch1`,1,2), FUN=function(x){
                                if(x %in% c("T1","T2", "T3", "T4", "TX")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(str_sub(clindata$`TNM (Pathological):ch1`,3,4), FUN=function(x){
                                if(x %in% c("NX", "N0", "N1", "N2", "N3")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`Stage (Pathological ):ch1`, FUN=function(x){
                                if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`EGFR status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "EGFR-WT"
                                }else if(x=="Mut"){
                                  "EGFR-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`K-ras Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "KRAS-WT"
                                }else if (x=="Mut"){
                                  "KRAS-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`p53 Status:ch1`, FUN=function(x){
                                if(x=="Wt"){
                                  "p53-WT"
                                }else if(x=="Mut"){
                                  "p53-Mut"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`Survival (days):ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`Status:ch1`, FUN=function(x){
                                if(x=="Alive"){
                                  "Alive"
                                }else if(x=="Dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`Evidence of relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="N"){
                                  "NO"
                                }else if(x=="Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Removing probes which contains NA values
na.ind<-which(apply(data.MA,MARGIN=1, FUN=function(x){
  length(which(is.na(x)))
})>0)
data.MA<-data.MA[-na.ind,]
data.Platform<-data.Platform[-na.ind,]

# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$GENE[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]


entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(na.omit(data.MA[entrez.to.id[[k]],])))
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE13213<-data.MA
adc.ind<-which(clinical_GSE13213$histology=="ADC")
clinical_GSE13213<-clinical_GSE13213[adc.ind,]
expression_GSE13213<-t(expression_GSE13213[, adc.ind])

# 5. GSE14814
data <- getGEO(filename="GEO/GSE14814_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE14814<-data.frame(Accession=clindata$geo_accession,
                              histology=sapply(clindata$`Histology type:ch1`, FUN=function(x){
                                if(x=="ADC"){
                                  "ADC"
                                }else if(x=="SQCC"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              gender=sapply(clindata$`Sex:ch1`, FUN=function(x){
                                if(x=="Female"){
                                  "female"
                                }else if(x=="Male"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`Stage:ch1`, FUN=function(x){
                                if(x=="I"){
                                  "Stage I"
                                }else if(x=="1B"){
                                  "Stage IB"
                                }else if(x=="2A"){
                                  "Stage IIA"
                                }else if(x=="2B"){
                                  "Stage IIB"
                                }else if(x=="II"){
                                  "Stage II"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              predominant_subtype=sapply(clindata$`predominant subtype:ch1`, FUN=function(x){
                                if(x=="Acinar"){
                                  "Acinar"
                                }else if(x=="Lepidic"){
                                  "Lepidic"
                                }else if(x=="Micropapillary"){
                                  "Micropapillary"
                                }else if(x=="invasive mucinous"){
                                  "Mucinous"
                                }else if(x=="Papillary"){
                                  "Papillary"
                                }else if(x %in% c("Solid", "solid")){
                                  "Solid"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              adjuvant_chemo=sapply(clindata$`Post Surgical Treatment:ch1`, FUN=function(x){
                                if(x=="OBS"){
                                  "no"
                                }else if(x=="ACT"){
                                  "yes"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`OS time:ch1`, FUN=function(x){
                                as.numeric(x)*12
                              }),
                              vital_status=sapply(clindata$`OS status:ch1`, FUN=function(x){
                                if(x=="Alive"){
                                  "Alive"
                                }else if(x=="Dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Removing probes which contains NA values
na.ind<-which(apply(data.MA,MARGIN=1, FUN=function(x){
  length(which(is.na(x)))
})>0)
data.MA<-data.MA[-na.ind,]
data.Platform<-data.Platform[-na.ind,]


# Mapping Probe ID to Entrez id
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
    return(colMeans(na.omit(data.MA[entrez.to.id[[k]],])))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE14814<-data.MA
adc.ind<-which(clinical_GSE14814$histology=="ADC")
clinical_GSE14814<-clinical_GSE14814[adc.ind,]
expression_GSE14814<-t(expression_GSE14814[, adc.ind])

### 6. GSE19188
# Lim (2017) data
clinical_Lim2017<-read.table("./LungCancer_validation_data/Clinicalinformation_processed.txt", header=T, sep="\t")

#### discovery data gene filtering by validation data genes
expression_Lim2017<-read.table("./LungCancer_validation_data/Combat_exprs.txt", header=T, sep="\t")
temp<-do.call(rbind, lapply(clinical_Lim2017$Accession, FUN=function(sample.id){
  return(expression_Lim2017[,which(sample.id==colnames(expression_Lim2017))])
}))

colnames(temp)<-rownames(expression_Lim2017)
rownames(temp)<-clinical_Lim2017$Accession
expression_Lim2017<-temp


### filter for ADC samples
adc.ind<-which(clinical_Lim2017$histology=="ADC")
clinical_Lim2017<-clinical_Lim2017[adc.ind, ]
for(i in 1:dim(clinical_Lim2017)[2]){
  clinical_Lim2017[which(is.na(clinical_Lim2017[,i])),i]<-"[Not Available]"
}
expression_Lim2017<-expression_Lim2017[adc.ind,]

## GSE19188
data <- getGEO(filename="GEO/GSE19188_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE19188<-clinical_Lim2017 %>% filter(Dataset=="GSE19188")
expression_GSE19188<-expression_Lim2017[clinical_GSE19188$Accession,]

# 7. GSE26939
data <- getGEO(filename="GEO/GSE26939_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE26939<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                x=setdiff(c(clindata$`sex:ch1`[k], clindata$`Sex:ch1`[k]), NA)
                                if(x=="F"){
                                  "female"
                                }else if(x=="M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking_status(0=nonsmoker,1=smoker):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="Smoker"){
                                  "Ever smoker"
                                }else if(x=="NeverSmoker"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(c(1:dim(clindata)[1]), FUN=function(k){
                                x=setdiff(c(clindata$`stage:ch1`[k], clindata$`Stage:ch1`[k]), NA)
                                if(length(x)==0){
                                  return("[Not Available]")
                                }
                                if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else if(x=="IIIB"){
                                  "Stage IIIB"
                                }else if(x=="IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`egfr (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "EGFR-Mut"
                                }else if(x=="0"){
                                  "EGFR-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "KRAS-Mut"
                                }else if(x=="0"){
                                  "KRAS-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              stk11_mut=sapply(clindata$`stk11 (0='wt',1='mutated'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="1"){
                                  "STK11-Mut"
                                }else if(x=="0"){
                                  "STK11-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              predimonant_subtype=sapply(clindata$`subtype:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }else{
                                  return(x)
                                }
                              }),
                              survival_time=sapply(clindata$`survival_months:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`survival_status(0='alive',1='dead'):ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="0"){
                                  "Alive"
                                }else if(x=="1"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)



data.MA=data@assayData$exprs
data.Platform<-data@featureData@data


# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$ORF[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[match(mapping.table$ID, rownames(data.MA)),]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)


# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE26939<-t(data.MA)

# 8. GSE29016
data <- getGEO(filename="GEO/GSE29016_series_matrix.txt.gz")
clindata <- data@phenoData@data
clinical_GSE29016<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$characteristics_ch1.5, FUN=function(x){
                                if(x=="Sex: F"){
                                  "female"
                                }else if(x=="Sex: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$`histology:ch1`, FUN=function(x){
                                if(x=="AC"){
                                  "ADC"
                                }else if(x=="histology: SqCC"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="smoker"){
                                  "Ever smoker"
                                }else if(x=="never_smoker"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`Stage:ch1`, FUN=function(x){
                                if(x=="IA"){
                                  "Stage IA"
                                }else if(x=="IB"){
                                  "Stage IB"
                                }else if(x=="IIA"){
                                  "Stage IIA"
                                }else if(x=="IIB"){
                                  "Stage IIB"
                                }else if(x=="IIIA"){
                                  "Stage IIIA"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(clindata$`tnm:ch1`, FUN=function(x){
                                t=substr(x, 1, 2)
                                if(t %in% c("T1", "T2", 'T3', "T4")){
                                  return(t)
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(clindata$`tnm:ch1`, FUN=function(x){
                                if(str_length(x)==7){
                                  n=substr(x, 4, 5)
                                }else if(str_length(x)==6){
                                  n=substr(x, 3, 4)
                                }else{
                                  n="[Not Available]"
                                }
                                
                                if(n %in% c("N0", "N1", "N2")){
                                  return(n)
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_m=sapply(clindata$`tnm:ch1`, FUN=function(x){
                                if(str_length(x)==7){
                                  m=substr(x, 6, 7)
                                }else if(str_length(x)==6){
                                  m=substr(x, 5, 6)
                                }else{
                                  m="[Not Available]"
                                }
                                
                                if(m %in% c("M0", "M1")){
                                  return(m)
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`egfr_mut:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x==1){
                                  "EGFR-Mut"
                                }else{
                                  "EGFR-WT"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras_mut:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x==1){
                                  "KRAS-Mut"
                                }else{
                                  "KRAS-WT"
                                }
                              }),
                              survival_time=sapply(clindata$`os:ch1`, FUN=function(x){
                                as.numeric(x)*12
                              }),
                              vital_status=sapply(clindata$`osbin:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="0"){
                                  "Alive"
                                }else if(x=="1"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

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
data.MA<-log2(data.MA)

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE29016<-t(data.MA)
adc.ind<-which(clinical_GSE29016$histology=="ADC")
clinical_GSE29016<-clinical_GSE29016[adc.ind,]
expression_GSE29016<-expression_GSE29016[adc.ind,]

# 9. GSE30219
data <- getGEO(filename="GEO/GSE30219_series_matrix.txt.gz")
clindata <- data@phenoData@data
clinical_GSE30219<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$characteristics_ch1.2, FUN=function(x){
                                if(x=="gender: F"){
                                  "female"
                                }else if(x=="gender: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$characteristics_ch1.3, FUN=function(x){
                                if(x=="histology: ADC"){
                                  "ADC"
                                }else if(x=="histology: SQC"){
                                  "SCC"
                                }else if(x=="histology: NTL"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_t=sapply(clindata$`pt stage:ch1`, FUN=function(x){
                                if(x %in% c("T1","T2", "T3", "T4", "TX")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_n=sapply(clindata$`pn stage:ch1`, FUN=function(x){
                                if(x %in% c("NX", "N0", "N1", "N2", "N3")){
                                  x
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`disease free survival in months:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  "[Not Available]"
                                }else if(x=="na"){
                                  "[Not Available]"
                                }else{
                                  x
                                }
                              }),
                              vital_status=sapply(clindata$`status:ch1`, FUN=function(x){
                                if(x=="ALIVE"){
                                  "Alive"
                                }else if(x=="DEAD"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              relapse=sapply(clindata$`relapse (event=1; no event=0):ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="0"){
                                  "NO"
                                }else if(x=="1"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              }),
                              disease_free_survival=sapply(clindata$`disease free survival in months:ch1`, FUN=function(x){
                                return(as.numeric(x))
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE30219<-data.MA
adc.ind<-which(clinical_GSE30219$histology=="ADC")
clinical_GSE30219<-clinical_GSE30219[adc.ind,]
expression_GSE30219<-t(expression_GSE30219[, adc.ind])

# 10. GSE31210
data <- getGEO(filename="GEO/GSE31210_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE31210<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$`gender:ch1`, FUN=function(x){
                                if(x=="female"){
                                  "female"
                                }else if(x=="male"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$`tissue:ch1`, FUN=function(x){
                                if(x=="primary lung tumor"){
                                  "ADC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`pathological stage:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x == "IA"){
                                  "Stage IA"
                                }else if(x == "IB"){
                                  "Stage IB"
                                }else if(x == "II"){
                                  "Stage II"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`gene alteration status:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="EGFR mutation +"){
                                  "EGFR-Mut"
                                }else{
                                  "EGFR-WT"
                                }
                              }),
                              kras_mut=sapply(clindata$`gene alteration status:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="KRAS mutation +"){
                                  "KRAS-Mut"
                                }else{
                                  "KRAS-WT"
                                }
                              }),
                              survival_time=sapply(clindata$`days before death/censor:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`death:ch1` , FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="alive"){
                                  "Alive"
                                }else if(x=="dead"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              disease_free_survival=sapply(clindata$`days before relapse/censor:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              relapse=sapply(clindata$`relapse:ch1`, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="not relapsed"){
                                  "NO"
                                }else if(x=="relapsed"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs %>% log(2)
data.Platform<-data@featureData@data
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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE31210<-t(data.MA)
adc.ind<-which(clinical_GSE31210$histology=="ADC")
clinical_GSE31210<-clinical_GSE31210[adc.ind,]
expression_GSE31210<-expression_GSE31210[adc.ind,]

# 11. GSE37745
data <- getGEO(filename="GEO/GSE37745_series_matrix.txt.gz")
clindata <- data@phenoData@data

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
                                  "NO"
                                }else if(x=="yes"){
                                  "YES"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              disease_free_survival=sapply(clindata$`days to recurrence / to last visit:ch1`, FUN=function(x){
                                as.numeric(x)
                              })
)



data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE37745<-data.MA
adc.ind<-which(clinical_GSE37745$histology=="ADC")
clinical_GSE37745<-clinical_GSE37745[adc.ind,]
expression_GSE37745<-t(expression_GSE37745[, adc.ind])


# 12. GSE41271
data <- getGEO(filename="GEO/GSE41271_series_matrix.txt.gz")
clindata <- data@phenoData@data
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
                              relapse=sapply(clindata$characteristics_ch1.9, FUN=function(x){
                                if(is.na(x))
                                  return("[Not Available]")
                                if(x=="recurrence: N"){
                                  "NO"
                                }else if(x=="recurrence: Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              }),
                              disease_free_survival=sapply(c(1:dim(clindata)[1]), FUN=function(k){
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
data.MA=data@assayData$exprs
data.Platform<-data@featureData@data
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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE41271<-t(data.MA)
adc.ind<-which(clinical_GSE41271$histology=="ADC")
clinical_GSE41271<-clinical_GSE41271[adc.ind,]
expression_GSE41271<-expression_GSE41271[adc.ind,]

# 13. GSE42127
data <- getGEO(filename="GEO/GSE42127_series_matrix.txt.gz")
clindata <- data@phenoData@data
clinical_GSE42127<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$characteristics_ch1.4, FUN=function(x){
                                if(x=="gender: F"){
                                  "female"
                                }else if(x=="gender: M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$characteristics_ch1.6, FUN=function(x){
                                if(x=="histology: Adenocarcinoma"){
                                  "ADC"
                                }else if(x=="histology: Squamous"){
                                  "SCC"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`final.pat.stage:ch1`, FUN=function(x){
                                if(x == "IA"){
                                  "Stage IA"
                                }else if(x == "IB"){
                                  "Stage IB"
                                }else if(x == "IIA"){
                                  "Stage IIA"
                                }else if(x == "IIB"){
                                  "Stage IIB"
                                }else if(x == "IIIA"){
                                  "Stage IIIA"
                                }else if(x == "IIIB"){
                                  "Stage IIIB"
                                }else if(x == "IV"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              adjuvant_chemo=sapply(clindata$`had_adjuvant_chemo:ch1`, FUN=function(x){
                                if(x=="FALSE"){
                                  "no"
                                }else if (x=="TRUE"){
                                  "yes"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`overall survival months:ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`survival status:ch1` , FUN=function(x){
                                if(x=="A"){
                                  "Alive"
                                }else if(x=="D"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data
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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE42127<-t(data.MA)
adc.ind<-which(clinical_GSE42127$histology=="ADC")
clinical_GSE42127<-clinical_GSE42127[adc.ind,]
expression_GSE42127<-expression_GSE42127[adc.ind,]

# 14. GSE50081
data <- getGEO(filename="GEO/GSE50081_series_matrix.txt.gz")
clindata <- data@phenoData@data

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
                                as.numeric(x)*12
                              }),
                              relapse=sapply(clindata$`recurrence:ch1`, FUN=function(x){
                                if(x=="N"){
                                  "NO"
                                }else if(x=="Y"){
                                  "YES"
                                }else
                                  "[Not Available]"
                              })
)

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data
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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}

if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE50081<-t(data.MA)
adc.ind<-which(clinical_GSE50081$histology=="ADC")
clinical_GSE50081<-clinical_GSE50081[adc.ind,]
expression_GSE50081<-expression_GSE50081[adc.ind,]

# 15.GSE63459
data <- getGEO(filename="GEO/GSE63459_series_matrix.txt.gz")
clindata <- data@phenoData@data

clinical_GSE63459<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$`Sex:ch1`, FUN=function(x){
                                if(x=="Female"){
                                  "female"
                                }else if(x=="Male"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              race=sapply(clindata$`race:ch1`, FUN=function(x){
                                if(x=="African American"){
                                  "black or african american"
                                }else if (x=="European American"){
                                  "white"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              histology=sapply(clindata$`histology:ch1`, FUN=function(x){
                                if(x=="adenocarcinoma"){
                                  "ADC"
                                }else if(x=="normal"){
                                  "Normal"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking status:ch1`, FUN=function(x){
                                if(x == "ever"){
                                  "Ever smoker"
                                }else if(x == "never"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              ajcc_pathologic_stage=sapply(clindata$`7th tnm stage:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x == "Ia"){
                                  "Stage IA"
                                }else if(x == "Ib"){
                                  "Stage IB"
                                }else if(x == "IIa"){
                                  "Stage IIA"
                                }else if(x == "IIb"){
                                  "Stage IIB"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras status:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="Mt"){
                                  "KRAS-Mut"
                                }else if(x=="WT"){
                                  "KRAS-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`tp53 mutation:ch1`, FUN=function(x){
                                if(is.na(x)){
                                  return("[Not Available]")
                                }
                                if(x=="Mt"){
                                  "p53-Mut"
                                }else if(x=="WT"){
                                  "p53-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$`survival after surgery (days):ch1`, FUN=function(x){
                                as.numeric(x)
                              }),
                              vital_status=sapply(clindata$`death due to cancer:ch1` , FUN=function(x){
                                if(x=="no"){
                                  "Alive"
                                }else if(x=="yes"){
                                  "Dead"
                                }else{
                                  "[Not Available]"
                                }
                              })
)
data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE63459<-t(data.MA)
adc.ind<-which(clinical_GSE63459$histology=="ADC")
clinical_GSE63459<-clinical_GSE63459[adc.ind,]
expression_GSE63459<-expression_GSE63459[adc.ind,]

# 16. GSE68465
series="GSE68465"
data <- getGEO(filename="GEO/GSE68465_series_matrix.txt.gz")
clinical_GSE68465<-read.csv("./LungCancer_validation_data/GSE68465_clinical.CSV")
clinical_GSE68465$survival_time<-sapply(c(1:length(clinical_GSE68465$vital_status)), FUN=function(x){
  if(clinical_GSE68465$last_contact_days_to[x] %in% na.terms){
    return("[Not Available]")
  }
  return(as.numeric(clinical_GSE68465$last_contact_days_to[x]))
})

clinical_GSE68465$race<-sapply(clinical_GSE68465$race, FUN=function(x){
  if(x=="Asian"){
    "asian"
  }else if(x=="Black or African American"){
    "black or african american"
  }else if(x=="Native Hawaiian or Other Pacific Islander"){
    "pacific islander"
  }else if (x=="White"){
    "white"
  }else{
    "[Not Available]"
  }
})
clinical_GSE68465$smoking_status<-sapply(clinical_GSE68465$tobacco_smoking_history_indicator, FUN=function(x){
  if(x %in% c("Smoked in the past", "Currently smoking")){
    "Ever smoker"
  }else if(x == "Never smoked"){
    "Never smoker"
  }else{
    "[Not Available]"
  }
})

clinical_GSE68465$grade<-sapply(clinical_GSE68465$histologic_grade, FUN=function(x){
  if(x=="Moderate Differentiation"){
    "Moderately"
  }else if(x=="POORLY DIFFERENTIATED"){
    "Poorly"
  }else if(x=="WELL DIFFERENTIATED"){
    "Well"
  }else{
    "[Not Available]"
  }
})

clinical_GSE68465$adjuvant_chemo<-sapply(clinical_GSE68465$pharmaceutical_tx_adjuvant, FUN=function(x){
  if(x=="NO"){
    "no"
  }else if (x=="YES"){
    "yes"
  }else{
    "[Not Available]"
  }
})

clinical_GSE68465$adjuvant_radiation<-sapply(clinical_GSE68465$radiation_treatment_adjuvant, FUN=function(x){
  if(x=="NO"){
    "no"
  }else if (x=="YES"){
    "yes"
  }else{
    "[Not Available]"
  }
})

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data
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
data.MA<-log2(data.MA)

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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE68465<-t(data.MA)

# 17. GSE72094
data <- getGEO(filename="GEO/GSE72094_series_matrix.txt.gz")
clindata <- data@phenoData@data
clinical_GSE72094<-data.frame(Accession=clindata$geo_accession,
                              gender=sapply(clindata$`gender:ch1`, FUN=function(x){
                                if(x=="F"){
                                  "female"
                                }else if(x=="M"){
                                  "male"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              race=sapply(clindata$`race:ch1`, FUN=function(x){
                                if(x %in% c("ASIAN INDIAN OR PAKISTANI", "VIETNAMESE", "THAI")){
                                  "asian"
                                }else if(x=="BLACK"){
                                  "black or african american"
                                }else if(x=="WHITE"){
                                  "white"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              smoking_status=sapply(clindata$`smoking_status:ch1`, FUN=function(x){
                                if(x== "Ever"){
                                  "Ever smoker"
                                }else if(x=="Never"){
                                  "Never smoker"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              
                              ajcc_pathologic_stage=sapply(clindata$characteristics_ch1.11, FUN=function(x){
                                if(x=="Stage: 1"){
                                  "Stage I"
                                }else if(x=="Stage: 1A"){
                                  "Stage IA"
                                }else if(x=="Stage: 1B"){
                                  "Stage IB"
                                }else if(x=="Stage: 2A"){
                                  "Stage IIA"
                                }else if(x=="Stage: 2B"){
                                  "Stage IIB"
                                }else if(x=="Stage: 3"){
                                  "Stage III"
                                }else if(x=="Stage: 3A"){
                                  "Stage IIIA"
                                }else if(x=="Stage: 3B"){
                                  "Stage IIIB"
                                }else if(x=="Stage: 4"){
                                  "Stage IV"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              egfr_mut=sapply(clindata$`egfr_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "EGFR-Mut"
                                }else if(x=="WT"){
                                  "EGFR-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              kras_mut=sapply(clindata$`kras_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "KRAS-Mut"
                                }else if(x=="WT"){
                                  "KRAS-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              p53_mut=sapply(clindata$`tp53_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "p53-Mut"
                                }else if(x=="WT"){
                                  "p53-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              stk11_mut=sapply(clindata$`stk11_status:ch1`, FUN=function(x){
                                if(x=="Mut"){
                                  "STK11-Mut"
                                }else if(x=="WT"){
                                  "STK11-WT"
                                }else{
                                  "[Not Available]"
                                }
                              }),
                              survival_time=sapply(clindata$characteristics_ch1.10, FUN=function(x){
                                res<-as.numeric(gsub("survival_time_in_days: ","", x))
                                if(is.na(res))
                                  res<-"[Not Available]"
                                return(res)
                              }),
                              vital_status=sapply(clindata$characteristics_ch1.9, FUN=function(x){
                                res<-gsub("vital_status: ","", x)
                                if(res=="NA")
                                  res<-"[Not Available]"
                                return(res)
                              }))


data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$EntrezGeneID[k], " ",fixed = T))
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

# MAD calculation
if(mad.exclusion!=0){
  mad.exp<-apply(data.MA, MARGIN=1, mad)
  data.MA<-data.MA[which(mad.exp>sort(mad.exp)[floor(length(mad.exp)*mad.exclusion)]),]
}
if(qnorm){
  data.MA<-normalize.quantiles(data.MA, keep.names=TRUE)
}

expression_GSE72094<-t(data.MA)

rm(data);rm(data.MA); rm(entrez.to.id);
rm(mapping.table); rm(adc.ind);
rm(alive.days.thres);rm(alive.month.thres);rm(death.days.thres);rm(death.month.thres)
