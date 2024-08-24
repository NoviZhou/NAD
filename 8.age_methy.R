####甲基化年龄 33种癌症
library(ChAMP)
library(data.table)
library(wateRmelon)"ACC","BLCA","BRCA","CESC","CHOL",
## 读取临床真实年龄
setwd("D:/NAD/age_calculate/甲基化年龄/")
CancerType<-c("COAD","DLBC","ESCA","GBM",
              "HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
              "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
              "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
setwd("F:/2022小丫/270 FigureYa270panMeth/")
for(cancer in CancerType){
  ## 快速加载DNA甲基化数据
  message("--",cancer,"...")
  load(file.path("tcga_methy450",paste0("TCGA-",cancer,"_methy450.RData"))) # 加载甲基化RData文件
  methy450 <- as.data.frame(methy450) # 转数据框
  methy450[1:3,1:3]
  rownames(methy450) <- methy450[,1]; methy450<- methy450[,-1]
  methy450[1:3,1:3] # 查看一下数据
  
  tum.sam <- colnames(methy450[,which(substr(colnames(methy450),14,15) == "01")]) # 提取原位癌样本
  nor.sam <- colnames(methy450[,which(substr(colnames(methy450),14,15) == "11")]) # 提取癌旁正常样本
  methy450 <- methy450[,c(tum.sam,nor.sam)] # 对数据做肿瘤、癌旁的排序
  methy450 <- as.data.frame(na.omit(methy450)) # 去除空值
  
  pd <- data.frame(Sample_Name = colnames(methy450),
                   Sample_Group = rep(c("T","N"),c(length(tum.sam),length(nor.sam))),
                   row.names = colnames(methy450),
                   stringsAsFactors = F)
  myFilter <- champ.filter(beta = as.matrix(methy450),
                           arraytype = "450K", # 平台使用的是甲基化450K
                           pd = pd,
                           filterNoCG = T, # 非CG位点被移除
                           filterSNPs = T, # 靠近SNP位点的CpG被移除
                           filterMultiHit = T, # 探针对应多个位点的被移除
                           filterXY = T, # 移除性染色体
                           fixOutlier = T, # 修正极端值
                           autoimpute = F) # 不填补缺失值
  # 速度较慢请耐心，我跑了30min左右
  myFilter <- champ.norm(as.matrix(myFilter$beta), 
                         method = "BMIQ", # 使用BMIQ法对数据进行标准化
                         plotBMIQ = F, # 不绘制相关标准化图像
                         cores = 1, # 单核运行
                         arraytype = "450K") # 平台采用甲基化450K
  save(myFilter,file = paste0("D:/NAD/age_calculate/甲基化年龄/data/",cancer,"_meth.filter.norm.rda")) # 保存数据，方便后续使用
   #计算甲基化年龄
  dnamage <- agep(myFilter[,tum.sam],
                       coeff=NULL)
  save(dnamage,file = paste0("D:/NAD/age_calculate/甲基化年龄/data/",cancer,"_dna_m_age.rda")) # 保存数据，方便后续使用
   # 计算生理年龄和实际临床年龄的差异
}
load(paste0("D:/NAD/age_calculate/甲基化年龄/data/ACC_dna_m_age.rda"))
