####做NAD相关lncRNA
library(stringr)
library(fgsea)
library(estimate)
Adjust_mrna <- function(mrna){
  mRow30 <- which(apply(mrna,1,function(v){return((sum(v==0)/length(v))>=0.3)}))
  mRemv <- mRow30
  if(length(mRemv)==0){
    mRNA_out0 <- mrna
  }else{
    mRNA_out0 <- mrna[-(mRemv),]
  }
  mRNA_exp_inter <- log2(mRNA_out0+0.001)
  return(mRNA_exp_inter)
}

Adjust_lncrna <- function(lncrna){
  lncRow50 <- which(apply(lncrna,1,quantile,probs=0.5)==0)
  lncRow90 <- which(apply(lncrna,1,quantile,probs=0.9)<=0.1)
  lncRemv <- union(lncRow50,lncRow90)
  if(length(lncRemv)==0){
    lncRNA_out0 <- lncrna
  }else{
    lncRNA_out0 <- lncrna[-(lncRemv),]
  }
  lncRNA_exp_inter <- log2(lncRNA_out0+0.001)
  return(lncRNA_exp_inter)
}
fun_mtx_pcr <- function(x,y,z){
  r12=cor(t(x),t(y))
  r13=cor(t(x),z)
  r23=cor(z,t(y))
  r123=r13%*%r23
  rup=r12-r123
  rd1=sqrt(1-r13*r13)
  rd2=sqrt(1-r23*r23)
  rd=rd1%*%rd2
  rm(r12, r13, r23, r123, rd1, rd2)
  gc()
  rrr=rup/rd
  return(rrr)
}
NAD_pathway<-list()
NAD_list<-read.table("D:/NAD/genelist/NAD_GO+wiki+reactome+KEGG.txt",header = T,sep = "\t")
names(NAD_list)[1]<-"SYMBOL"
NAD_gene<-read.table("D:/NAD/genelist/NAD_GO+wiki+reactome+KEGGgene.txt",header = T,sep = "\t")

for(pathway in unique(NAD_list$pathway)){
  temp<-NAD_list[which(NAD_list$pathway==pathway),]
  gene<-merge(NAD_gene,temp,all=F)
  NAD_pathway[[pathway]] <-gene$ENSEMBL
}
#save(NAD_pathway, file = "D:/NAD/genelist/NAD_pathway.Rdata")

gtf_mrna_v22 <- read.table("F:/小丫画图/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/gtfmRNA22.txt", header = T, sep = "\t")
head(gtf_mrna_v22)
gtf_lnc_v22 <- read.table("F:/小丫画图/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/gtflncRNA22.txt", header = T, sep = "\t")
head(gtf_lnc_v22)

CancerType<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
for(cancer in CancerType){
load(file=paste0("D:/pan-cancer/TCGA-",cancer,"/TCGA-",cancer,"/data_exprs_tpm.Rdata"))
# mRNA
# 用于immune pathway的富集
cancer_mrna <- data_exprs_tpm[rownames(data_exprs_tpm) %in% gtf_mrna_v22$gene_id, ]
cancer_mrna <- Adjust_mrna(cancer_mrna)
#write.table(cancer_mrna, "cancer_mrna.txt", sep = "\t")
cancer_mrna2 <- merge(gtf_mrna_v22[, c("gene_id", "gene_name")], data_exprs_tpm, by.x = 1, by.y = 0) %>%
  .[, -1]
cancer_mrna2 <- aggregate(.~gene_name, cancer_mrna2, max)
rownames(cancer_mrna2) <- cancer_mrna2$gene_name; cancer_mrna2 <- cancer_mrna2[, -1]
cancer_mrna2 <- Adjust_mrna(cancer_mrna2)
write.table(cancer_mrna2, paste0("E:/R-4.0.2/library/estimate/extdata/TCGA_",cancer,"_purityinput.txt"), sep = "\t", quote = F)

# lncRNA
cancer_lnc <- merge(gtf_lnc_v22[, c("gene_id", "gene_name")], data_exprs_tpm, by.x = 1, by.y = 0) %>% .[, -1]
cancer_lnc <- aggregate(.~gene_name, cancer_lnc, max)
rownames(cancer_lnc) <- cancer_lnc[, 1]; cancer_lnc <- cancer_lnc[, -1]
cancer_lnc <- Adjust_lncrna(cancer_lnc)
cancer_mrna[1:3, 1:3]
cancer_lnc[1:3, 1:3]
# estimate包计算浸润水平
#file.copy(from = paste0("F:/小丫画图/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/purityinput.txt"),to="E:/R-4.0.2/library/estimate/extdata/purityinput.txt")
setwd("E:/R-4.0.2/library/estimate/extdata/")
CancerExpr<-system.file("extdata", paste0("TCGA_",cancer,"_purityinput.txt"), package="estimate")
filterCommonGenes(input.f=CancerExpr, output.f="exp_deal.gct", id="GeneSymbol")
estimateScore("exp_deal.gct", "estimate_score_all.gct", platform="illumina")
est_score_all <- readLines("estimate_score_all.gct")
est_scores <- unlist(strsplit(est_score_all[grep("ESTIMATEScore",est_score_all)],"\t")) %>% 
  .[3:length(.)] %>% as.numeric()

#肿瘤纯度
Tumour_purity <- cos(0.6049872018+0.0001467884*est_scores)
samples_name <- unlist(strsplit(est_score_all[grep("NAME",est_score_all)],"\t")) %>% 
  .[3:length(.)]
names(Tumour_purity) <- samples_name %>% str_replace_all(., "\\.", "-")
#去除一些不要用的中间变量
rm(cancer_mrna2, samples_name, est_score_all, est_scores)
#这一步非常耗费内存:生成的相关矩阵是行为lncRNA，列为mRNA的矩阵 

all(colnames(cancer_lnc) == colnames(cancer_mrna))
all(colnames(cancer_mrna) == names(Tumour_purity))
n = ncol(cancer_lnc) #样本数目
gn = 1 #偏相关矫正的影响因素个数，这里只有肿瘤纯度一个。
pcor <- fun_mtx_pcr(cancer_lnc,cancer_mrna, Tumour_purity)
#计算统计量和P值
statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
p.value <- 2*pnorm(-abs(statistic))
rownames(pcor) <- rownames(cancer_lnc) ; rownames(p.value) <- rownames(cancer_lnc)
colnames(pcor) <- rownames(cancer_mrna) ; colnames(p.value) <- rownames(cancer_mrna)
#计算RS: 也就是按照pvalue排序，根据正负相关赋予pvalue一个方向
RS <- -log10(p.value)*sign(pcor)

#load("F:/小丫画图/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/FigureYa170ImmuLncRNA/pathways.rda")
k=0.995
fgseaRes_all <- c()

for(i in 1:nrow(RS)){
  ##remove the rows contain Infs
  if(sum( is.infinite(RS[i,])) != 0){
    next()
  }
  ranks <- RS[i,]
  fgseaRes <- fgsea(NAD_pathway, ranks, minSize=1, maxSize=5000, nperm=1000)
  sigValue <- c()
  for(j in 1:nrow(fgseaRes)){
    if(fgseaRes$ES[j]>0){
      sig_ij <- 1 - 2*fgseaRes$pval[j]
    }else{
      sig_ij <- 2*fgseaRes$pval[j] - 1
    }
    sigValue <- c(sigValue,sig_ij)
    
  }
  lncRNA <- rownames(RS)[i]
  fgseaRes_i <- cbind(lncRNA,fgseaRes,sigValue)
  fgseaRes_all <- rbind(fgseaRes_all,fgseaRes_i)
}
write.table(fgseaRes_all[,-9],paste0("D:/NAD/lncRNA/all_lnc_pathway/NAD_all_lncRNA_",cancer,".txt"),quote = F,sep = "\t",row.names = F)
#筛选符合标准的gsea结果
sig_ind <- which(abs(fgseaRes_all$sigValue) >= k & fgseaRes_all$padj < 0.05)
#即可求出lncRNA-pathway pairs
sig_pairs <-  as.data.frame(fgseaRes_all[sig_ind,])
write.table(sig_pairs[,-9] ,paste0("D:/NAD/lncRNA/NAD_lncRNA_",cancer,".txt"),quote = F,sep = "\t",row.names = F)
}

####画图
##统计
# 加载前面计算获得的NAD相关lncRNA数量，用于top y-axis
plot_input<-data.frame()
i=1
for(cancer in CancerType){
  
  sig<-read.table(paste0("D:/NAD/lncRNA/NAD_lncRNA_",cancer,".txt"),header=T,sep = "\t")
  all<-read.table(paste0("D:/NAD/lncRNA/all_lnc_pathway/NAD_all_lncRNA_",cancer,".txt"),header=T,sep = "\t")
  plot_input[i,1]<-cancer
  plot_input[i,2]<- length(unique(sig$lncRNA)) 
  plot_input[i,3]<- length(unique(all$lncRNA)) 
  plot_input[i,4]<- length(unique(sig$lncRNA))/length(unique(all$lncRNA)) 
  i=i+1
  }
colnames(plot_input) <- c("Cancer_type", "Sig", "All","Prop")

write.csv(plot_input, "D:/NAD/lncRNA/plot/plot_input.txt", quote = F, row.names = F)


pdat <- read.csv("D:/NAD/lncRNA/plot/plot_input.txt")
head(pdat)
pdat$Cancer_type <- factor(pdat$Cancer_type, levels = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
                                                        "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
                                                        "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
                                                        "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))
scaleFactor <- max(pdat$Sig)/ max(pdat$Prop)
library(ComplexHeatmap)
library(RColorBrewer) 
library(ggplot2)

#colName<-c("#847786","#878E9B","#8EB392","#D2B675","#F4B58B","#854A5E")
colName<-c("#EFD5D4","#FDD4D9","#A3B4CA","#BADDF1","#7EB5DD","#73ABD8","#5C92C2","#FAE0F5","#F885AE","#C16685","#C07B43","#A2DA9C","#9EBED3")
cancercolor<-colorRampPalette(colName)(33)
mycols <- c(colorRampPalette(brewer.pal(5, "YlOrRd"))(3),
            colorRampPalette(brewer.pal(9, "PuRd"))(9),
            colorRampPalette(brewer.pal(9, "Blues"))(9), 
            colorRampPalette(brewer.pal(9, "Purples"))(9),
            colorRampPalette(brewer.pal(5, "YlGn"))(3))


mycols<-cancercolor
p <- ggplot(pdat, aes(x=pdat$Cancer_type)) + 
  # 画柱形图
  geom_bar(aes(y = pdat$Sig / scaleFactor , fill = pdat$Cancer_type), stat = "identity", size = 0.8, width = 0.7) +
  scale_fill_manual(values = mycols) + #自定义bar的颜色
  
  # 画折线图
  geom_point(aes(y = pdat$Prop), shape = 21, fill = "white") +
  geom_line(aes(y = pdat$Prop, group = 1)) +
  
  scale_x_discrete(expand = expand_scale(mult = c(0.04,0.04)), #上下留空
                   limits = rev(levels(pdat$Cancer_type))) + 
  scale_y_continuous(name = "Proportions of NAD-related lncRNAs",  
                     sec.axis=sec_axis(trans = ~. * scaleFactor, 
                                       name = "Number of NAD-related lncRNAs",
                                       breaks = c(100,200,300,400)),
                     expand = c(0, 0)) + #画右侧Y轴
  labs(x="") + 
  theme_bw() + #去除背景色
  theme(panel.grid = element_blank()) + #去除网格线
  theme(legend.position = "none") + #不画图例
  coord_flip() #坐标轴互换，因此图中x轴对应纵坐标癌症名
p
setwd("D:/NAD/lncRNA/plot/")
ggsave("NADLncRNA.pdf", width = 5, height = 8)
####
