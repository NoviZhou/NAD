##GTEx pca主成分然后与年龄做相关#
###显著性module与age
#按照PMID 35000600的方法age的时间段取中位数
library(readr)
library(tidyr)
library(rtracklayer)
GTEx_MEGENA<-read_tsv("D:/NAD/MEGENA/GTEx_调整年龄前/GTEx_module.tsv/GTEx_module.tsv")
head(GTEx_MEGENA)
rm(res_enrich_all)
res_enrich_all<-c()
res_sig_all<-read.table("D:/NAD/coexpression/GTEx/MEGENA_sig.txt",header = T,sep="\t")
####gtf对应各自表达谱
#gtf <- rtracklayer::import('D:/GTEx/gencode.v23.annotation.gtf/gencode.v23.annotation.gtf')
#gtf <- as.data.frame(gtf)
#gtf <- gtf[gtf$type=="gene",c("gene_id","gene_name","gene_type")]
#write.table(gtf,"D:/GTEx/Gencodev23ENSG2SYMBOL.txt",row.names = F,col.names = T,sep="\t",quote=F)
ENSG_SYM<-read.table("D:/GTEx/Gencodev23ENSG2SYMBOL.txt",sep="\t",header = T)
names(ENSG_SYM)[1:2]<-c("Ensembl","Symbol")
rm(cor_matrix)
cor_matrix <- data.frame()
i=1
data_clinical<-read.table("D:/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",header=T,sep="\t")
data_clinical<-separate(data_clinical,col="AGE", into=c("min","max"),sep ="-" )
data_clinical$age<-(as.numeric(data_clinical$min)+as.numeric(data_clinical$max))/2
for(tissue in unique(res_sig_all$tissue)){
  tissue_sub<-GTEx_MEGENA[which(GTEx_MEGENA$Tissue==tissue),]
  enrich_module<-res_sig_all[which(res_sig_all$tissue==tissue),]
  TPM<-read.table(paste0("D:/GTEx/GTEx_SMTSD/",tissue,"_tpm_expr.txt"),row.names = 1,header = T,check.names = F)
  for(module in unique(enrich_module$Module)){
    modulegene_GTEx_MEGENA<-tissue_sub[which(tissue_sub$Module==module),]
    esmb<-merge(ENSG_SYM,modulegene_GTEx_MEGENA,all=F)
    ensg<-intersect(esmb$Ensembl,row.names(TPM))
    expr<-log2(TPM[ensg,]+1)
    pca.results <- prcomp(expr, center = TRUE, scale. = FALSE)
    pca.rotation <- pca.results$rotation
    pca.rotation
    low_dim_df <- as.data.frame(pca.rotation[,1])
    low_dim_df$sample<-row.names(low_dim_df)
    low_dim_df<-separate(low_dim_df,col="sample",into =c("a","b","c","d","e") ) 
    low_dim_df$SUBJID<-paste(low_dim_df$a,low_dim_df$b,sep="-")
    names(low_dim_df)[1]<-"PCA1"
    low_dim_df<-low_dim_df[,-c(2,3,4,5,6)]
    pca_age<-merge(data_clinical,low_dim_df,all=F)
    cor_test_res <- cor.test(pca_age$PCA1,pca_age$age,method="spearman" )
    cor_matrix[i,1] <- tissue
    cor_matrix[i,2] <- module
    cor_matrix[i,3] <- cor_test_res$estimate
    cor_matrix[i,4] <- cor_test_res$p.value
    i=i+1
  }
}
names(cor_matrix)<-c("Cancer","Module","R_value","P_value")
write.table(cor_matrix,"D:/NAD/coexpression/GTEx/cor_module_age.txt",row.names = F,col.names = T,sep="\t",quote=F)

sig_cor_matrix<-cor_matrix[cor_matrix$P_value<0.05,]
write.table(sig_cor_matrix,"D:/NAD/coexpression/GTEx/cor_module_age_sig.txt",row.names = F,col.names = T,sep="\t",quote=F)
