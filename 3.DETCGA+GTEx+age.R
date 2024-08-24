
CancerType<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
##NAD按照年龄的结果

TCGA_DEResults_CT <- readRDS("D:/NAD/expression/TCGA_DEResults_CT.rds")
NAD<-read.table("D:/NAD/genelist/NAD_GO+wiki+reactome+KEGGgene.txt",quote="",header=T)
NAD_wide<-NAD_DE_age[,1:2]#i=1先跑一下赋值
for(i in 1:length(names(TCGA_DEResults_CT))){
  DE_age<-TCGA_DEResults_CT[[i]]
  NAD_DE_age<-DE_age[intersect(row.names(DE_age),NAD$SYMBOL),]
  NAD_wide[,i]<-NAD_DE_age$P.Value
  write.table(NAD_DE_age,paste("D:/NAD/expression/age_sig/",names(TCGA_DEResults_CT)[i],"_TCGA_NAD_age.txt",sep=""),sep="\t",quote=F)
    }

names(NAD_wide)<-names(TCGA_DEResults_CT)
write.table(NAD_wide,"D:/NAD/expression/age_sig/NAD_age_pancancer.txt",sep="\t",quote=F)  

####不分癌症一起做
TCGA_DEResults_all <- readRDS("D:/NAD/expression/TCGADE_Symbols.rds")
NAD_TCGA_DEResults_all<-TCGA_DEResults_all[intersect(row.names(TCGA_DEResults_all),NAD$SYMBOL),]
write.table(NAD_TCGA_DEResults_all,"D:/NAD/expression/age_sig/NAD_TCGA_DEResult.txt",sep="\t",quote=F)  
###GTEx
GTEx_DEResults_all<-readRDS("D:/NAD/expression/GTExDE_Symbols.rds")
NAD_GTEx_DEResults_all<-GTEx_DEResults_all[intersect(row.names(GTEx_DEResults_all),NAD$SYMBOL),]
write.table(NAD_GTEx_DEResults_all,"D:/NAD/expression/age_sig/NAD_GTEx_DEResults_all.txt",sep="\t",quote=F)  

##GTEx计算每个组织的表达
setwd("D:/NAD/expression/GTEx/")
NAD_list<-list.files("D:/NAD/expression/GTEx/")
expr_GTEx<-expr_tissue#i=1先跑一下赋值
for(i in 1:length(NAD_list)){
NAD_expr<-read.table(NAD_list[i],header = T,sep="\t",row.names=1)
row.names(NAD_expr)<-NAD_expr$SYMBOL
NAD_expr<-NAD_expr[,-1]
expr_tissue<-as.data.frame(rowMeans(as.matrix(NAD_expr)))
names(expr_tissue)<-as.character(unlist(strsplit(NAD_list[i],split="_")))[1] 
expr_GTEx<-cbind(expr_GTEx,expr_tissue)  
}

write.table(expr_GTEx[,-1],"D:/NAD/expression/NAD_GTEx_pertissue.txt",sep="\t",quote = F,row.names = T)
expr_GTEx<-read.table("D:/NAD/expression/NAD_GTEx_pertissue.txt",sep="\t",quote = "")


###Deseq DE
#对所有进行差异表达

CancerType<-c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LIHC","LUAD",
              "LUSC","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","THCA","THYM","UCEC")
length(CancerType)
setwd("D:/RNA修饰/Counts(bugai)/Counts(bugai)/")
####无正常对照样本
#CancerType<-c("ACC","DLBC","LAML","LGG","MESO","OV","TGCT","UCS","UVM")
#length(CancerType)
library('DESeq2')
#NAD<-read.table("D:/NAD/genelist/NADwiki+GOgene.txt",quote="",header=T)
myMeanFun<-function(x){
  tapply(as.double(x),group,mean)  
}   
rm(allresultfile_all)
allresultfile_all<-c()

library(tidyr)
for(i in 1:length(CancerType)){
  expr_data<- read.table(paste("TCGA.",CancerType[i],".counts.ori.bugai.txt",sep=""),header=TRUE,sep='\t',check.names = FALSE,stringsAsFactors=FALSE) 
  expr_data<-separate(data =expr_data, col = gene_id, into = c("ENSEMBL", "banbenhao"), sep = "[.]")
  expr_regulator<-merge(NAD,expr_data,all=F)
  row.names(expr_regulator)<- expr_regulator$SYMBOL
  expr_regulator<-expr_regulator[,-c(1:4)]
  ##group分组
  group <- ifelse(as.numeric(substr(colnames(expr_regulator),14,15))<10,"Cancer","Normal");table(group)
  condition <- factor(group)
  coldata <- data.frame(row.names = colnames(expr_regulator), condition)
  dds <- DESeqDataSetFromMatrix(countData=expr_regulator, colData=coldata, design=~condition)
  dds <- DESeq(dds)
  res <- results(dds,pAdjustMethod = "fdr",contrast=c("condition", "Cancer","Normal"))
  ###calculate cancer and normal mean 
  meangroup <- t(apply(t(expr_regulator),2,myMeanFun))
  allresult <- data.frame(meangroup[rownames(res),],res)
  #allresult<-allresult[order(as.numeric(allresult[,"pvalue"])),]
  allresult<-data.frame(ID=row.names(allresult),allresult,stringsAsFactors = FALSE);head(allresult)
  allresultfile<-data.frame(as.character(allresult$ID),
                            CancerType[i],
                            signif(allresult$Cancer,4),
                            signif(allresult$Normal,4),
                            signif(allresult$log2FoldChange,4),
                            signif(allresult$pvalue,4),
                            signif(allresult$padj,4),stringsAsFactors = FALSE)
  colnames(allresultfile)<-c("Gene","Cancer Type","Cancer_MeanExpression","Normal_MeanExpression","log2FoldChange","PValue","FDR")
  allresultfile_all<-rbind(allresultfile_all,allresultfile)
  write.table(allresultfile,paste("D:/NAD/expression/TCGA_DEseq/",CancerType[i],"_Different_expression.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
  }
allresultfile_all<-allresultfile_all[order(allresultfile_all$Gene),]
write.table(allresultfile_all,"D:/NAD/expression/TCGA_DEseq_Different_expression.txt",quote=F,col.names = T,row.names = F,sep="\t")
