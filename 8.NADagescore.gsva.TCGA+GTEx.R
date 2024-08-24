###gsea
source("D:/NAD/程序/0.GSVA.R")
library(GSVA)
CancerType<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

# 计算signature score -- 基于GSVA
load("D:/NAD/signature_score/1.gene_list_NAD_GTEx_TCGA.ensg.RData")
setwd("D:/NAD/signature_score/GSVA_score/TCGA+GTEx")
for(cancer in CancerType){
  load(file=paste0("D:/pan-cancer/TCGA-",cancer,"/TCGA-",cancer,"/data_exprs_tpm.Rdata"))
  indata <- as.matrix(log2(data_exprs_tpm+0.001))
  #开始进行ssGSEA
  gsva_NAD_list <-gsva(indata ,gene_list_NAD_age, method="ssgsea",kcdf='Gaussian',parallel.sz=24) 
  write.csv(gsva_NAD_list ,paste0("ssgsea/",cancer,".csv"))
}
for(cancer in CancerType){
  load(file=paste0("D:/pan-cancer/TCGA-",cancer,"/TCGA-",cancer,"/data_exprs_tpm.Rdata"))
  indata <- as.matrix(log2(data_exprs_tpm+0.001))
  #开始进行ssGSEA
  gsva_NAD_list <-gsva(indata ,gene_list_NAD_age, method="gsva",kcdf='Gaussian',parallel.sz=24) 
  write.csv(gsva_NAD_list ,paste0("gsva/",cancer,".csv"))
}
for(cancer in CancerType){
  load(file=paste0("D:/pan-cancer/TCGA-",cancer,"/TCGA-",cancer,"/data_exprs_tpm.Rdata"))
  indata <- as.matrix(log2(data_exprs_tpm+0.001))
  #开始进行ssGSEA
  gsva_NAD_list <-gsva(indata ,gene_list_NAD_age, method="zscore",kcdf='Gaussian',parallel.sz=24) 
  write.csv(gsva_NAD_list ,paste0("zscore/",cancer,".csv"))
}
for(cancer in CancerType){
  load(file=paste0("D:/pan-cancer/TCGA-",cancer,"/TCGA-",cancer,"/data_exprs_tpm.Rdata"))
  indata <- as.matrix(log2(data_exprs_tpm+0.001))
  #开始进行ssGSEA
  gsva_NAD_list <-gsva(indata ,gene_list_NAD_age, method="plage",kcdf='Gaussian',parallel.sz=24) 
  write.csv(gsva_NAD_list ,paste0("plage/",cancer,".csv"))
}
####以上都尝试输出  计算score
for(cancer in CancerType){
  gsva<-read.csv(paste0("gsva/",cancer,".csv"),row.names = 1)
  gsva[7,]<-log(exp((gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,]))/exp((gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,])))
  gsva[8,]<-log(exp((gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,]))/exp((gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,])))
  gsva[9,]<-(8*(gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,])-4*(gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,]))/12
  gsva[10,]<-(25*(gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,])-7*(gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,]))/32
  # gsva[9,]<-(8*gsva[2,]-4*gsva[3,])/12
  #  gsva[10,]<-(25*gsva[5,]-7*gsva[6,])/32
  row.names(gsva)[7:10]<-c("score_diff1","score_diff2","score1","score2")
  write.csv(gsva,paste0("gsva/",cancer,".csv"))
}
for(cancer in CancerType){
  gsva<-read.csv(paste0("ssgsea/",cancer,".csv"),row.names = 1)
  gsva[7,]<-log(exp((gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,]))/exp((gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,])))
  gsva[8,]<-log(exp((gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,]))/exp((gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,])))
  gsva[9,]<-(8*(gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,])-4*(gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,]))/12
  gsva[10,]<-(25*(gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,])-7*(gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,]))/32
# gsva[9,]<-(8*gsva[2,]-4*gsva[3,])/12
#  gsva[10,]<-(25*gsva[5,]-7*gsva[6,])/32
  #row.names(NAD_score)[7:8]<-c("score1","score2")
  row.names(gsva)[7:10]<-c("score_diff1","score_diff2","score1","score2")
  write.csv(gsva,paste0("ssgsea/",cancer,".csv"))
}
for(cancer in CancerType){
  gsva<-read.csv(paste0("plage/",cancer,".csv"),row.names = 1)
  gsva[7,]<-log(exp((gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,]))/exp((gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,])))
  gsva[8,]<-log(exp((gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,]))/exp((gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,])))
  gsva[9,]<-(8*(gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,])-4*(gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,]))/12
  gsva[10,]<-(25*(gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,])-7*(gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,]))/32
  # gsva[9,]<-(8*gsva[2,]-4*gsva[3,])/12
  #  gsva[10,]<-(25*gsva[5,]-7*gsva[6,])/32
  #row.names(NAD_score)[7:8]<-c("score1","score2")
  row.names(gsva)[7:10]<-c("score_diff1","score_diff2","score1","score2")
  write.csv(gsva,paste0("plage/",cancer,".csv"))
}
for(cancer in CancerType){
  gsva<-read.csv(paste0("zscore/",cancer,".csv"),row.names = 1)
  gsva[7,]<-log(exp((gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,]))/exp((gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,])))
  gsva[8,]<-log(exp((gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,]))/exp((gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,])))
  gsva[9,]<-(8*(gsva[2,]-mean(as.numeric(gsva[2,])))/sd(gsva[2,])-4*(gsva[3,]-mean(as.numeric(gsva[3,])))/sd(gsva[3,]))/12
  gsva[10,]<-(25*(gsva[5,]-mean(as.numeric(gsva[5,])))/sd(gsva[5,])-7*(gsva[6,]-mean(as.numeric(gsva[6,])))/sd(gsva[6,]))/32
  # gsva[9,]<-(8*gsva[2,]-4*gsva[3,])/12
  #  gsva[10,]<-(25*gsva[5,]-7*gsva[6,])/32
  #row.names(NAD_score)[7:8]<-c("score1","score2")
  row.names(gsva)[7:10]<-c("score_diff1","score_diff2","score1","score2")
  write.csv(gsva,paste0("zscore/",cancer,".csv"))
}
