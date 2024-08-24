###ctrp
CancerType<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
ctrp.candidate.all<-c()
prism.candidate.all<-c()
for(cancer in CancerType){
  ctrp.log2fc.p<- read.csv(paste0("D:/NAD/drug/NADctrp/",cancer,"ctrp.csv"))
  #ctrp.log2fc.p<-ctrp.log2fc.p[ctrp.log2fc.p$p.value<0.05,]
  ctrp.log2fc.p<-ctrp.log2fc.p[abs(ctrp.log2fc.p$low.high)>0.1,]
  ctrp.spearman<- read.csv(paste0("D:/NAD/drug/NADctrp/",cancer,"ctrp.Spearman.csv"))
  ctrp.spearman<-ctrp.spearman[ctrp.spearman$P<0.05,]
  ctrp.spearman<-ctrp.spearman[abs(ctrp.spearman$R)>0.35,]
  ctrp.candidate<-merge(ctrp.log2fc.p,ctrp.spearman,all=F)
  print(cancer)
  print(dim(ctrp.candidate)[1])
  if(dim(ctrp.candidate)[1]>0){
    ctrp.candidate$cancer<-cancer
    ctrp.candidate.all<-rbind(ctrp.candidate.all,ctrp.candidate)
  }
  
  prism.log2fc.p<- read.csv(paste0("D:/NAD/drug/NADprism/",cancer,"prism.csv"))
  #prism.log2fc.p<-prism.log2fc.p[prism.log2fc.p$p.value<0.05,]
  prism.log2fc.p<-prism.log2fc.p[abs(prism.log2fc.p$low.high)>0.1,]
  prism.spearman<- read.csv(paste0("D:/NAD/drug/NADprism/",cancer,"prism.Spearman.csv"))
  prism.spearman<-prism.spearman[prism.spearman$P<0.05,]
  prism.spearman<-prism.spearman[abs(prism.spearman$R)>0.35,]
  prism.candidate<-merge(prism.log2fc.p,prism.spearman,all=F)
  print(cancer)
  print(dim(prism.candidate)[1])
  if(dim(prism.candidate)[1]>0){
    prism.candidate$cancer<-cancer
    prism.candidate.all<-rbind(prism.candidate.all,prism.candidate)
  }
}
prism.candidate.all<-separate(prism.candidate.all,col=drug,into = c("drug","b"),sep="[ ]")[,-2]
write.csv(prism.candidate.all,"D:/NAD/drug/NAD_prism_original.csv",quote = F,row.names = F)
write.csv(ctrp.candidate.all,"D:/NAD/drug/NAD_ctrp_original.csv",quote = F,row.names = F)


#drugAnn1<-read.csv("F:/½ô³Ã/data/DrugData/DrugData/GeneDrugInformation.csv")
#drugAnn1<- drugAnn1[drugAnn1$associate_genes %in% unique(PSM_table$gene),]

#intersect(drugAnn1$generic_name,ctrp.candidate.all$drug)
#intersect(drugAnn1$generic_name,prism.candidate.all$drug)
