library(TCGAbiolinks)
library(maftools)
library(stringr)
CancerType<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
for(cancer in CancerType){
  clinical <- GDCquery(project=paste("TCGA-",cancer,sep=""),
                       data.category = "Clinical",
                       file.type = "xml")
  GDCdownload(clinical)
  cliquery <- GDCprepare_clinic(clinical,clinical.info = "patient")
  colnames(cliquery)[1] <- "Tumor_Sample_Barcode"
  mut <- GDCquery_Maf(tumor = cancer,pipelines = "mutect2") #下载突变数据
  mut$Tumor_Sample_Barcode[1]
  mut$Tumor_Sample_Barcode <- str_sub(mut$Tumor_Sample_Barcode,1,12)
  young<-cliquery[which(abs(cliquery$days_to_birth)/365<=50),1]
  old<-cliquery[which(abs(cliquery$days_to_birth)/365>50),1]
  mut.young <- mut[(mut$Tumor_Sample_Barcode %in% young),]
  mut.old <- mut[(mut$Tumor_Sample_Barcode %in% old),]
  maf.young <- read.maf(maf=mut.young,isTCGA=T)
  maf.old <- read.maf(maf=mut.old,isTCGA=T)
  pt.vs.rt <- mafCompare(m1 =maf.young, m2 = maf.old, m1Name = 'young', m2Name = 'old', minMut = 5)
  print(pt.vs.rt)

  write.table(pt.vs.rt$results,paste(cancer,"Young0ldCompare.txt"),col.names=T,row.names = F,quote = F)
  write.table(pt.vs.rt$SampleSummary,paste(cancer,"Young0ldCompareSample.txt"),col.names=T,row.names = F,quote = F)
  print(cancer)
  }