library(readxl)
setwd('C:\\Users\\조휘인\\Desktop\\2024 SP 개별연구\\GEO_Data_analysis')
result_cox_gse31210<-read.csv('gse31210 uni_cox_analysis_results.csv')
result_cox_gse37745<-read.csv('gse37745 uni_cox_analysis_results.csv')
result_cox_gse41271<-read.csv('gse41271 uni_cox_analysis_results.csv')
result_cox_gse500081<-read.csv('gse50081 uni_cox_analysis_results.csv')



#####
#Results of pvalue box plot으로 comparison
boxplot(sort(result_cox_gse31210$p.value)[1:50],sort(result_cox_gse37745$p.value)[1:50],sort(result_cox_gse41271$p.value)[1:50],sort(result_cox_gse500081$p.value)[1:50]
,xlab='GEO data',ylab='p-value',main='p-value of each group')
abline(h=4.425*10^-6,col='red')

#각 그룹 간 RFS 및 not rfs 환자 비율
table(rfs_data_GSE31210$status);table(rfs_data_GSE41271$status);table(rfs_data_GSE50081$status)

par(mfrow=c(2,2))
pie(table(rfs_data_GSE31210$status),
    labels = paste(round(table(rfs_data_GSE31210$status)/sum(table(rfs_data_GSE31210$status)),2)*100,'%'),main='Composition of GSE31210')
pie(table(rfs_data_GSE37745$status),
    labels = paste(round(table(rfs_data_GSE37745$status)/sum(table(rfs_data_GSE37745$status)),2)*100,'%'),main='Composition of GSE37745')
pie(table(rfs_data_GSE41271$status),
    labels = paste(round(table(rfs_data_GSE41271$status)/sum(table(rfs_data_GSE41271$status)),2)*100,'%'),main='Composition of GSE41271')
pie(table(rfs_data_GSE50081$status),
    labels = paste(round(table(rfs_data_GSE50081$status)/sum(table(rfs_data_GSE50081$status)),2)*100,'%'),main='Composition of GSE50081')

#각 그룹 간 rfs 및 last contact day 비교
non_rfs_31210<-subset(rfs_data_GSE31210,rfs_data_GSE31210$status==1);rfs_31210<-subset(rfs_data_GSE31210,rfs_data_GSE31210$status==2)
non_rfs_37745<-subset(rfs_data_GSE37745,rfs_data_GSE37745$status==1);rfs_37745<-subset(rfs_data_GSE37745,rfs_data_GSE37745$status==2)
non_rfs_41271<-subset(rfs_data_GSE41271,rfs_data_GSE41271$status==1);rfs_41271<-subset(rfs_data_GSE41271,rfs_data_GSE41271$status==2)
non_rfs_50081<-subset(rfs_data_GSE50081,rfs_data_GSE50081$status==1);rfs_50081<-subset(rfs_data_GSE50081,rfs_data_GSE50081$status==2)
par(mfrow=c(1,1))
boxplot(non_rfs_31210$time,non_rfs_37745$time,non_rfs_41271$time,non_rfs_50081$time,xlab='GSE data group',ylab='last contact days',main='Last contact days of non-nte patients')
boxplot(rfs_31210$time,rfs_37745$time,rfs_41271$time,rfs_50081$time,xlab='GSE data group',ylab='RFS days',main='RFS of nte patients')

#각 그룹간 nte vs not-nte days(time) 비교
par(mfrow=c(2,2))
boxplot(non_rfs_31210$time,rfs_31210$time,xlab='GSE31210',ylab='days(time)',main='Time comparision b/w non nte and nte')
boxplot(non_rfs_37745$time,rfs_37745$time,xlab='GSE37745',ylab='days(time)',main='Time comparision b/w non nte and nte')
boxplot(non_rfs_41271$time,rfs_41271$time,xlab='GSE41271',ylab='days(time)',main='Time comparision b/w non nte and nte')
boxplot(non_rfs_50081$time,rfs_50081$time,xlab='GSE50081',ylab='days(time)',main='Time comparision b/w non nte and nte')
boxplot()
boxplot(non_rfs_31210$G100616237,non_rfs_37745$G100616237,non_rfs_50081$G100616237)

####

####
#Outlier를 어떻게 제거할 것인가?
quantile(non_rfs_31210$time)



