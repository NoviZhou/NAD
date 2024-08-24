# 设置目录
data_home <- "D:/NAD/"
res_home <- "D:/NAD/"
source(paste0("D:/NAD/程序/0.checkDir.R"))
checkDir(paste0(res_home, "PSM.omics/RData"))
checkDir(paste0(res_home, "PSM.omics/pdf"))
checkDir(paste0(res_home, "PSM.omics/table"))
# 导入特需函数
source(paste0("D:/NAD/程序/0.cal.R"))
source(paste0("D:/NAD/程序/0.GIPW_function_omega.R"))
source(paste0("D:/NAD/程序/0.check_balance.R"))


# 导入数据
load(file=paste0(res_home, "PSM.omics/RData/dum_form_list.RData"))
load(file=paste0("D:/pan-cancer/psm.data/data_exprs_list.RData"))

# 设置并行计算
library(dummies)
library(doMC)
library(foreach)
registerDoMC(24)

# library(dummies)
folder <- "1.mRNA"
scripts.dir <- paste0(res_home, "PSM.omics/",folder)
checkDir(scripts.dir)

# setwd("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/")
# if (!file.exists(folder)) { dir.create(folder) }

setwd(scripts.dir)
# analysis没有什么用，只是用来占位和目录命名的
sum.mRNAAll <- data.frame()
analysis<-"NADscore"
for(cancer in names(dum_form_list)){
  print(paste("Running for", cancer, "..."))
  # 提取数据
  form <- dum_form_list[[cancer]]$form
  data.dum <- dum_form_list[[cancer]]$data.dum
  mRNAseq.pri <- log2(t(data_exprs_list[[cancer]])+1)
  mRNAseq.pri <- rm.zero.col(mRNAseq.pri)

  outdir_tmp <- paste(scripts.dir, "/",cancer,"_",analysis,sep="")
  checkDir(outdir_tmp)
  mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE, weight="MW", mirror.plot=FALSE, cancer, data.type= "mRNAseq", outdir=outdir_tmp, perm=FALSE)
  sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE)
  summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE, cutoff=0.05)
  write.summary(sum.mRNA, cancer, analysis,"mRNA")
  write.result(mRNAseq.result, cancer, analysis,"mRNA")
  save(mRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(mRNAseq.result$fdr < 0.05)) > 0){
    sum.mRNA <- data.frame(sum.mRNA)
    sum.mRNA$class <- rep(cancer,times=nrow(sum.mRNA))
    if(nrow(sum.mRNAAll) == 0){
      sum.mRNAAll <- sum.mRNA
    }else{
      sum.mRNAAll <- rbind(sum.mRNAAll,sum.mRNA)
    }
    perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
    {
      ## mRNAseq, only for KIRC        
      perm.mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight="MW",mirror.plot=FALSE, cancer, "mRNA", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
      perm.sum.mRNA <- summarize.fdr(mRNAseq.pri, perm.mRNAseq.result)
      
      write(c(seed, perm.sum.mRNA$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
      save(seed,perm.mRNAseq.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
      
    }
    
    
    cutoff <- 0.05
    seedV <- 1:100
    perm.cal(cancer, analysis, "mRNAseq", mRNAseq.pri, cutoff=cutoff, seedV=seedV)
    
    if(FALSE)
    {
      mRNA.ttest <- myttest(data.dum, mRNAseq.pri, cancer,"mRNA")
      sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNA.ttest)
      save(mRNAseq.ttest, file=paste(cancer,"_mRNA_ttest.RData", sep=""))
    }
  }
}

write.table(sum.mRNAAll,file="mRNAseqlog2.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)

