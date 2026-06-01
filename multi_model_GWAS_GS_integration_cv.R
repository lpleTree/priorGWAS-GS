source("gapit_functionsw.txt")
library(openxlsx)
library(tidyverse)
library(data.table)
library(ggplot2)
library(caret)
library(dplyr)
library(caret)
library(Matrix)
#source("http://zzlab.net/GAPIT/gapit_functions.txt")

GD = fread(file = "sample849_012tab.txt", header = TRUE)
myY = read.table(file = "杨树849个样本表型数据.txt", header = TRUE)[ ,c(1,2)]   #每次只输入一个性状
GM = read.table(file = "GM.txt", header = TRUE)
#sum(is.na(GD));sum(is.na(myY));sum(is.na(GM))
# #########GWAS 分析不同的模型######
myGAPIT=GAPIT(
  Y=myY[,c(1,2)], #fist column is ID
  GD=GD,
  GM=GM,
  PCA.total=5,
  model=c("FarmCPU", "BLINK", "MLM", "GLM"),
  Multiple_analysis=TRUE)


files <- list.files(pattern = "^GAPIT\\.Association\\.Filter_GWAS_results\\.csv$")# 查找以 GAPIT.Association.Filter_GWAS_results.csv 开头的文件，在确定有显著位点情况下才能进行

if (length(files) > 0) {    #检查是否找到文件
  file_to_read <- files[1]  #如果有多个匹配文件，选择第1个
  p.sig_sum <- read.csv(file_to_read)
}

summary_counts <- p.sig_sum %>%  #统计第7列的类别数量
  group_by(across(7)) %>%
  summarise(count = n())

most_common_category <- summary_counts %>%
  arrange(desc(count)) %>%
  slice(1) %>%
  pull(1)

best_model <- sub("\\..*$", "", most_common_category[[1]])       #从第一列中提取.号之前的字符

dir.create("differ-stra")  #建立新的文件夹，设置不同策略的工作路径，查看路径，执行不同策略
setwd("differ-stra")
getwd()

############################################################
###############编写基于GWAS辅助的GS预测函数#################
############################################################
GWAS_assisted_GS <- function(GD, GM, myY, models, p.levels, n_repeats, n_folds) {
  results <- data.frame(models = character(),repeats = integer(), folds = integer(), PandSNPnumber = numeric(), strategy = character(), r = numeric())  #创建一个数据框来存储结果
  set.seed(123)  # 保证结果可重复
  folds_indices <- createMultiFolds(myY[, 2], k = n_folds, times = n_repeats)
  names(folds_indices) <- gsub("Fold0+([1-9])", "Fold\\1", names(folds_indices)) #正则表达式操作来格式化结果名称。例如，始终移除Fold后面多余的零。
  
  for (m in 1:length(models)) {  #外部循环：重复次数
    #m=1
    model <- models[m] 
    for (repeats in 1:n_repeats) {
      #repeats=1
      for (fold in 1:n_folds) {
        #fold=1
        Fold.Rep <- as.name(paste0("Fold", fold, ".Rep", repeats))  #获取当前折的索引作为测试集
        train_indices <- folds_indices[[Fold.Rep]]
        test_indices <- setdiff(1:nrow(myY),  train_indices)        #除去测试集的索引，剩余的样本作为训练集
        source("gapit_functions王佳博.txt")
        myGAPIT = GAPIT(                                            #在训练群体中进行GWAS分析，7个模型均进行，其中有一个是最佳模型
          Y = myY[train_indices, c(1,2)], # 第一列是 ID
          GD = GD,
          GM = GM,
          PCA.total = 5,
          model = c("FarmCPU", "BLINK", "MLM", "GLM"),
          Multiple_analysis = F)
        
       source("./gapit_functions.txt")
        #source("./gapit_functions20250824.txt")
        for (p.level in p.levels) {    #针对不同的 p.level 进行策略选择
          #p.level=0.01
          all_files <- list.files()    #从当前文件夹读取最佳模型的结果
          pattern <- paste0("GAPIT.Association.GWAS_Results.", best_model, ".*")
          target_files <- grep(pattern, all_files, value = TRUE)
          
          for (file in target_files) {  #读取每个符合条件的文件
            file_path <- file.path(getwd(), file)
            best_model_GWAS <- read.csv(file_path)
          }
          
          ##策略2，根据GWAS的结果选择不同P水平的SNP子集
          p02 = best_model_GWAS[, c(1, 4)]
          selected_indices02 <- which(p02$P.value < p.level)    #根据不同P值筛选标记的索引
          selected_snps_p02 <- p02[selected_indices02, 1]         #根据索引提取标记的名称
          subset_matrix_p_GD02 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p02))
          subset_matrix_p_GM02 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD02), ]   #提取GM矩阵
          
          stra2 <- GAPIT(
            Y=myY[train_indices,c(1,2)],
            GD=subset_matrix_p_GD02, 
            GM=subset_matrix_p_GM02,
            PCA.total=8,
            model="model", 
            file.out=F)
          
          m2=merge(myY[test_indices,c(1,2)],stra2$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
          R2=cor(m2[,5],m2[,2])^2
          r2=sqrt(R2)
          results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str2", r = r2)) 

          
          ##策略4，以上一步最佳模型作为训练群体的GWAS分析模型，将P最小的作为固定效应，根据P值构建不同P.level的SNP子集
          best_sig_snp_1 <- best_model_GWAS[which.min(best_model_GWAS$P.value), ] #筛选最显著的SNP--筛选最小的 P.value 所对应的第一列的SNP位点编号
          snp_value <- best_sig_snp_1$SNP
          best_sig_matrix_GD_1 <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_value)))   #根据筛选的 SNP 构建子集GD矩阵作为固定效应

          p04 = best_model_GWAS[, c(1, 4)]
          selected_indices04 <- which(p04$P.value < p.level)    #根据不同P值筛选标记的索引
          selected_snps_p04 <- p04[selected_indices04, 1]         #根据索引提取标记的名称
          subset_matrix_p_GD04 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p04))
          subset_matrix_p_GM04 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD04), ]   #提取GM矩阵


          pattern <- paste0("GAPIT.Association.PVE.", best_model, ".*")
          best.model.sigSNP <- list.files("../", pattern, full.names = TRUE)
          all.gwas <- read.csv(best.model.sigSNP)
          best_sig_snp_name <- unique(all.gwas$SNP)        #根据索引提取标记的名称
          best_sig_snp_GD <- data.frame(GD %>% dplyr::select(Taxa, all_of( best_sig_snp_name))) #提取GD矩阵
          stra4 <- GAPIT(
            Y=myY[train_indices,c(1,2)],
            GD=subset_matrix_p_GD04,
            GM=subset_matrix_p_GM04,
            PCA.total=5,
            CV=best_sig_snp_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
            model="model",
            SNP.test=FALSE,
            memo="MAS+model",file.out=F)

          m4=merge(myY[test_indices,c(1,2)],stra4$Pred[,c(1,3,5,11)],by.x="Taxa",by.y="Taxa")
          R4=cor(m4[,5],m4[,2])^2
          r4=sqrt(R4)
          results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str4", r = r4))


          ##策略6,在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的位点的并集作为固定效应，根据P值构建不同P.level的SNP子集
          p06 = best_model_GWAS[, c(1, 4)]
          selected_indices06 <- which(p06$P.value < p.level)    #根据不同P值筛选标记的索引
          selected_snps_p06 <- p06[selected_indices06, 1]         #根据索引提取标记的名称
          subset_matrix_p_GD06 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p06)) #提取GD矩阵
          subset_matrix_p_GM06 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD06), ]   #提取GM矩阵

          all.gwas <- read.csv("../GAPIT.Association.Filter_GWAS_results.csv")
          u.selected_snps <- unique(all.gwas$SNP)        #根据索引提取标记的名称
          U_subset_matrix_GD <- data.frame(GD %>% dplyr::select(Taxa, all_of(u.selected_snps))) #提取GD矩阵
          stra6 <- GAPIT(
            Y=myY[train_indices,c(1,2)],
            GD=subset_matrix_p_GD06,
            GM=subset_matrix_p_GM06,
            PCA.total=5,
            CV=U_subset_matrix_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
            model="model",
            SNP.test=FALSE,
            memo="MAS+model",file.out=F)

          m6=merge(myY[test_indices,c(1,2)],stra6$Pred[,c(1,3,5,11)],by.x="Taxa",by.y="Taxa")
          R6=cor(m6[,5],m6[,2])^2
          r6=sqrt(R6)

          # 将策略6的结果添加到数据框中
          results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str6", r = r6))

        }
      }
    }
    
  }
  return(results)
}

#执行函数
result <- GWAS_assisted_GS(GD, GM, myY,
                           models =c("sBLUP", "gBLUP", "cBLUP"), #
                           p.levels <- c(0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5,1.0),
                           n_repeats = 1, 
                           n_folds =10)
write.csv(result, "traitq_stra2和4和6_P值-20251023补充性状1.csv")



