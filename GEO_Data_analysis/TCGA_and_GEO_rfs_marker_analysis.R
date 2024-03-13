head(sort(as.numeric(result_cox_gse31210$p.value)))



rfs_marker_gse31210<-rownames(subset(result_cox_gse31210,as.numeric(result_cox_gse31210$p.value)<0.001))
rfs_marker_gse37745<-rownames(subset(result_cox_gse37745,as.numeric(result_cox_gse37745$p.value)<0.01))
rfs_marker_gse41271<-rownames(subset(result_cox_gse41271,as.numeric(result_cox_gse41271$p.value)<0.01))
rfs_marker_gse50081<-rownames(subset(result_cox_gse50081,as.numeric(result_cox_gse50081$p.value)<0.01))

sum(rfs_marker_gse37745%in%rfs_marker_gse41271==TRUE);rfs_marker_gse37745[rfs_marker_gse37745%in%rfs_marker_gse41271==TRUE]
sum(rfs_marker_gse37745%in%rfs_marker_gse50081==TRUE);rfs_marker_gse37745[rfs_marker_gse37745%in%rfs_marker_gse50081==TRUE]
sum(rfs_marker_gse41271%in%rfs_marker_gse50081==TRUE);rfs_marker_gse41271[rfs_marker_gse41271%in%rfs_marker_gse50081==TRUE]
sum(rfs_marker_gse31210%in%rfs_marker_gse37745==TRUE);rfs_marker_gse31210[rfs_marker_gse31210%in%rfs_marker_gse37745==TRUE]
sum(rfs_marker_gse31210%in%rfs_marker_gse41271==TRUE);rfs_marker_gse31210[rfs_marker_gse31210%in%rfs_marker_gse41271==TRUE]
sum(rfs_marker_gse31210%in%rfs_marker_gse50081==TRUE);rfs_marker_gse31210[rfs_marker_gse31210%in%rfs_marker_gse50081==TRUE]
org.Hs.eg()

keys(org.Hs.eg.db,keytype = 'ENSEMBL')

names(org.Hs.eg.db)
View(org.Hs.eg.db)
