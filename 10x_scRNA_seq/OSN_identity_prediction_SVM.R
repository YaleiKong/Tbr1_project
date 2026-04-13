rm(list=ls())

# install.packages("LiblineaR")
# install.packages("e1071")
# install.packages("caret")   # 混淆矩阵、训练辅助

library(LiblineaR)
library(e1071)
# library(caret)
library(pROC)
library(PRROC)

## Load required packages
.libPaths(new = "/lustre/home/acct-medlqian/medlqian-loop3/R/x86_64-redhat-linux-gnu-library/4.2.2")

library(Seurat) # 文章中用的v4.1.0
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(purrr) 
library(tidyverse) 
library(stringr)
library(RColorBrewer)
library(patchwork)
library(colorspace) # darken
library(harmony)
library(ggvenn)
library(Matrix)
library(clustree)
library(ComplexUpset)
library(openxlsx)

set.seed(42)

setwd("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/TAAR_OR_identity_model/")

#### 1. 读取对象并定义 TAAR / OR 标签 ####
# 提取WT单细胞数据中的表达单个受体的mOSN
load("../20250518_WT_10X_scRNA_Seq_seurat_analysis/integrated_analysis_seuratV5/rdata/mOSN_expressing_single_receptor.rdata")

# 不区分 class I OR 和 class II OR
mOSN_scRNA_single_receptor$Taar_OR_class <- ifelse(
  mOSN_scRNA_single_receptor$receptor_class == "Taar",
  "Taar",
  "OR"
)

Idents(mOSN_scRNA_single_receptor) <- "Taar_OR_class"

#### 2. 先划分训练/测试 ####
meta <- mOSN_scRNA_single_receptor@meta.data
meta$cell_id <- rownames(meta)
meta$label <- mOSN_scRNA_single_receptor$Taar_OR_class

## 联合分层：orig.ident + label
meta$strata <- paste(meta$orig.ident, meta$label, sep = "_")

## 分层抽取 80% 进入训练集
train_cells <- unlist(lapply(split(meta$cell_id, meta$strata), function(cells) {
  n <- length(cells)
  
  n_train <- max(1, min(n - 1, floor(0.8 * n)))
  
  return(sample(cells, size = n_train))
}))
# split(meta$cell_id, meta$strata) 这一步是把所有细胞按 meta$strata 分组。
# lapply(..., function(cells) { ... }) 对每一个分组，分别执行后面的抽样函数。

test_cells <- setdiff(meta$cell_id, train_cells)

## 检查训练/测试标签分布
table(meta[train_cells, "label"])
# OR  Taar 
# 46992   576 

table(meta[train_cells, "receptor_class"])
# class_I_OR class_II_OR        Taar 
# 5704       41288         576 

table(meta[test_cells, "label"])
# OR  Taar 
# 11757   148 

#### 3. 仅在训练集上计算 DEGs (Taar vs OR) ####
# 导入差异基因

#### 4. 提取表达矩阵 ####
train_markers_use <- read.xlsx("../20250518_WT_10X_scRNA_Seq_seurat_analysis/integrated_analysis_seuratV5/files/Taar_OR_train_only_DE_results_mOSN_expressing_single_receptor.xlsx",
                               sheet = "train_markers_filtered")

features <- train_markers_use$gene

length(features)
# 191

mat <- GetAssayData(mOSN_scRNA_single_receptor, assay = "SCT", layer = "data")

features_exclude_Tbr1 <- features[-which(features == "Tbr1")]

#### 5. 构造 x ####
# 把表达矩阵变成 cells × genes，符合 glmnet 输入格式
x1_train <- t(as.matrix(mat[features_exclude_Tbr1, train_cells]))
x1_test  <- t(as.matrix(mat[features_exclude_Tbr1, test_cells]))

## 可选：去掉训练集中零方差基因
keep <- apply(x1_train, 2, function(z) sd(z) > 0)
x1_train <- x1_train[, keep, drop = FALSE]
x1_test  <- x1_test[, keep, drop = FALSE]

dim(x1_train)
# 47568   190

#### 6. 构造 y 构造二分类标签 ####
y_train <- ifelse(meta[train_cells, "label"] == "Taar", 1, -1)

y_test  <- ifelse(meta[test_cells, "label"] == "Taar", 1, -1)

#### 7. 用训练集均值和标准差做标准化 ####
# SVM 对特征尺度敏感。e1071::svm() 默认会做缩放；如果你用 LiblineaR，通常自己先标准化更稳妥。e1071 文档也明确说明其默认会对变量做 scaling。
# 因为 SVM（以及很多基于距离/内积的模型）对特征尺度很敏感，先标准化是为了让每个基因在模型里处于可比较的量纲，避免某些数值范围大的基因不成比例地支配分类边界。
# 标准化后，每个基因都被拉到类似尺度，模型才更公平地比较它们的贡献。

# 本质上是在做每个基因的 z-score 标准化：标准化后均值变成 0 标准化后标准差变成 1
train_center <- colMeans(x1_train)
train_scale  <- apply(x1_train, 2, sd)
train_scale[train_scale == 0] <- 1

x_train_sc <- scale(x1_train, center = train_center, scale = train_scale)
x_test_sc  <- scale(x1_test,  center = train_center, scale = train_scale)

#### 8. 先训练一个 Linear SVM ####
fit_svm <- LiblineaR(
  data = x_train_sc,
  target = y_train,
  type = 1,      # linear SVM的一种常用分类设置
  cost = 1,
  bias = TRUE
)

#### 9. 预测测试集 ####
pred_cls <- predict(fit_svm, x_test_sc, decisionValues = TRUE)

# 预测类别
pred_label <- pred_cls$predictions

y_test <- ifelse(meta[test_cells, "label"] == "Taar", 1, -1)

table(pred_label, y_test)
# y_test
# pred_label    -1     1
# -1 11744    17
# 1     13   131

#### 10. 交叉验证选参数 ####
## 对一组不同的 cost 参数，进行 5 折交叉验证，比较每个 cost 的平均分类准确率，最后选出表现最好的 cost。

## 1. 构建分层 5 折 fold_id
set.seed(42)

nfolds <- 5
fold_id <- rep(NA, length(y_train))

# 正类（例如 TAAR = 1）
idx_pos <- which(y_train == 1)
fold_id[idx_pos] <- sample(rep(1:nfolds, length.out = length(idx_pos)))

# 负类（例如 OR = -1）
idx_neg <- which(y_train == -1)
fold_id[idx_neg] <- sample(rep(1:nfolds, length.out = length(idx_neg)))

# 检查每折类别分布
print(table(fold_id, y_train))
print(prop.table(table(fold_id, y_train), margin = 1))

## 2. 对不同 cost 做分层 5 折 CV
cost_grid <- c(0.01, 0.1, 1, 10, 100)
cv_acc <- numeric(length(cost_grid))
cv_f1 <- numeric(length(cost_grid))

for (i in seq_along(cost_grid)) {
  cost_i <- cost_grid[i]
  
  fold_acc <- c()
  
  for (k in 1:nfolds) {
    tr_idx <- which(fold_id != k)
    va_idx <- which(fold_id == k)
    
    center_k <- colMeans(x1_train[tr_idx, , drop = FALSE])
    scale_k  <- apply(x1_train[tr_idx, , drop = FALSE], 2, sd)
    scale_k[scale_k == 0] <- 1
    
    x_train_sc <- scale(
      x1_train[tr_idx, , drop = FALSE],
      center = center_k,
      scale  = scale_k
    )
    
    fit_k <- LiblineaR(
      data   = x_train_sc,
      target = y_train[tr_idx],
      type   = 1,
      cost   = cost_i,
      bias   = TRUE
    )
    
    x_valid_sc <- scale(
      x1_train[va_idx, , drop = FALSE],
      center = center_k,
      scale  = scale_k
    )
    
    pred_k <- predict(fit_k, x_valid_sc)
    
    acc_k <- mean(pred_k$predictions == y_train[va_idx])
    fold_acc <- c(fold_acc, acc_k)
    
    y_pred <- pred_k$predictions
    y_true <- y_train[va_idx]
    
    TP <- sum(y_pred == 1  & y_true == 1)
    FP <- sum(y_pred == 1  & y_true == -1)
    FN <- sum(y_pred == -1 & y_true == 1)
    
    precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
    recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
    fold_f1[k] <- ifelse(precision + recall == 0, 0,
                         2 * precision * recall / (precision + recall))
  }
  
  cv_acc[i] <- mean(fold_acc)
  cv_f1[i] <- mean(fold_f1)
}

## 3. 查看结果并选 best_cost
cv_result <- data.frame(
  cost = cost_grid,
  cv_acc = cv_acc,
  cv_f1 = cv_f1
)

print(cv_result)
# cost    cv_acc     cv_f1
# 1 1e-02 0.9971199 0.8701635
# 2 1e-01 0.9972461 0.8793262
# 3 1e+00 0.9968887 0.8658168
# 4 1e+01 0.9957955 0.8248795
# 5 1e+02 0.9944291 0.7909579

best_cost <- cost_grid[which.max(cv_acc)]
print(best_cost)
# 0.1

#### 5. 构造 x ####
# 把表达矩阵变成 cells × genes，符合 glmnet 输入格式
x1_train <- t(as.matrix(mat[features_exclude_Tbr1, train_cells]))
x1_test  <- t(as.matrix(mat[features_exclude_Tbr1, test_cells]))

## 可选：去掉训练集中零方差基因
keep <- apply(x1_train, 2, function(z) sd(z) > 0)
x1_train <- x1_train[, keep, drop = FALSE]
x1_test  <- x1_test[, keep, drop = FALSE]

dim(x1_train)
# 47568   190

#### 6. 构造 y 构造二分类标签 ####
y_train <- ifelse(meta[train_cells, "label"] == "Taar", 1, -1)

y_test  <- ifelse(meta[test_cells, "label"] == "Taar", 1, -1)

#### 7. 用训练集均值和标准差做标准化 ####
# SVM 对特征尺度敏感。e1071::svm() 默认会做缩放；如果你用 LiblineaR，通常自己先标准化更稳妥。e1071 文档也明确说明其默认会对变量做 scaling。
# 因为 SVM（以及很多基于距离/内积的模型）对特征尺度很敏感，先标准化是为了让每个基因在模型里处于可比较的量纲，避免某些数值范围大的基因不成比例地支配分类边界。
# 标准化后，每个基因都被拉到类似尺度，模型才更公平地比较它们的贡献。

# 本质上是在做每个基因的 z-score 标准化：标准化后均值变成 0 标准化后标准差变成 1
train_center <- colMeans(x1_train)
train_scale  <- apply(x1_train, 2, sd)
train_scale[train_scale == 0] <- 1

x_train_sc <- scale(x1_train, center = train_center, scale = train_scale)
x_test_sc  <- scale(x1_test,  center = train_center, scale = train_scale)

#### 11. 重新训练最终模型 (best_cost) ####
final_svm <- LiblineaR(
  data   = x_train_sc,
  target = y_train,
  type   = 1,
  cost   = best_cost,
  bias   = TRUE
)
saveRDS(final_svm, file = "rdata/final_svm.rds")

w <- as.vector(final_svm$W)
names(w) <- colnames(x_train_sc)
w
write.csv(as.data.frame(w), file = "files/linear_SVM_genes_weight.csv")

# Psmd14       Gucy1b2        Tmbim6           Ckb         Hint1        Stoml3        Tbc1d4         Acsm4       Sult1d1         Tafa4          Nsg2         Stbd1 
# 0.3697319468  0.0027607671  0.1461062825  0.1941029484  0.3102004547 -0.4595040116  0.0117881687 -0.0245313959  0.1425821168  0.0583682821  0.0574735727  0.0775598598 
# Ndfip1         Mslnl          Ppa1         Acss2         Tenm2         Rgs16        Stmnd1         Myo1b          Snca           Ptn        Actr3b           Cd9 
# 0.1426182878  0.0614108744 -0.0823766696  0.0530317569  0.0322913734  0.0259332265  0.0125095440  0.0729662101  0.0580023536  0.0336553126 -0.0099614584 -0.0685597172 
# A430035B10Rik         Psmc5          Nfix        Cartpt       Slc27a2        Adam28 2810001G20Rik          Ebf1          Rtn1         Lsamp          Ece1         Fetub 
# 0.0347306995 -0.0143939452 -0.0055686821  0.0076382515  0.0347423990  0.1088523169 -0.0217962659 -0.1091339981 -0.0697041013  0.0919203334  0.0574300730 -0.0387959816 
# Ccn3         Rexo2         Tenm4        Mrpl27         Unc5b          Peg3        Homer2       Fam129a          Faim      H2-M10.3          Ech1          Cd55 
# 0.0077537785 -0.0065584184 -0.0100838766  0.0440463533 -0.0223521243  0.0044671066  0.0209594764  0.0098480430 -0.0731853812 -0.0224788298 -0.0491118111 -0.0319349464 
# Pdlim1         Cebpb         Olig1         Acss3        Ssx2ip        Sorbs2         Prrg4       Smarcd2       Tspan12        Papss2         Lactb         Cidea 
# -0.0314202497 -0.0546191673 -0.0715755201  0.0144248023  0.0331932180 -0.0369400117  0.0401052553  0.0098097085  0.0065376432 -0.0385119866  0.0144780481 -0.1001384464 
# Ucp2          Glul          Pth2         Plcb4          Dmkn          Syt7          Bmp6        Kcnmb3       Ankrd10        Tnni3k         Pcdh7          Gldc 
# -0.0595552724 -0.0495445764  0.0759758493 -0.0028547281 -0.0566983645 -0.0405921093  0.0873407313  0.0359866997 -0.0043373919 -0.0076207815 -0.1810860825 -0.0598456407 
# Celf4          Tpi1         Dock8     Hist2h2be       Htatsf1          Nfib        Dbndd2         Auts2         Rpl11        Sema7a       Mkrn2os            F8 
# 0.0080992513  0.0018731482  0.0154625467  0.0821397085  0.0139403558  0.0742758529 -0.0109123205 -0.0539610968 -0.0548803840 -0.0689459208  0.0321322221 -0.0204740393 
# Lgmn         Tshz1         Hjurp          Mitf        Plaat3        Trim45          Pfn2 6430548M08Rik        Rnf185       Gm13589         Sept3          Aff3 
# 0.0099963743 -0.0686383527 -0.0613598527  0.0078453868  0.0499227335 -0.0738336338  0.0202115519 -0.0252151715 -0.0099671862 -0.1008931708 -0.0244050764 -0.0105960931 
# Cflar          Map4         Tppp3         Rims3         Cldn3         Ifi27       Gm31557 4930515B02Rik        Arxes2       Gm13648    AC149090.1         Pcdh9 
# 0.0440077580  0.0673808957 -0.0565783608 -0.0529291613 -0.0218809921 -0.0482906829 -0.1111834123  0.0135523288  0.0389523543 -0.0198951538 -0.0558048011 -0.0702868433 
# B230118H07Rik          Phf1        Elavl2         Olig2          Rbm4         Tubb3         Med27        Gm6878         Cotl1         Gdf11        Mboat1        Bicdl2 
# -0.0011650345  0.0080798935  0.0284021584 -0.0782198565  0.0006005545  0.0090989333  0.0438586111 -0.0173814844 -0.0295451986  0.0064236854 -0.0542578899  0.0279337156 
# Bcl11b      Slc25a33       Osbpl1a        Cox7a1       Pcdh11x      Hectd2os        Fam89a      Tmem132a           Cpe      Ppp1r15a          Caly       Gm16579 
# -0.0518218085 -0.0132152239 -0.0293784004  0.0275978102  0.0254284272 -0.0212553577  0.0140812074 -0.1105548401 -0.0447674679 -0.0228141020 -0.0541797339 -0.0684492049 
# Ssbp2         Wdr54        Cyb5r1         Trib3       Gm16196           Mt2        Pde10a      Cacna2d3       Phactr1        Cab39l          Btg1         Mef2b 
# -0.0648270769  0.0006598536 -0.0244747349  0.0085609886 -0.0054593708  0.0199924482  0.0287930425 -0.0366656285 -0.0570358120 -0.0565381207  0.0394042302 -0.0570400752 
# Cd36          Nqo1         Sgms1      Pafah1b3        Erlin1         Prpf4         Kcnh7          Rhoq        Plxna1         Rab39       Ppp2r3a         Clip4 
# -0.0550955028 -0.3775338705 -0.0050119243 -0.0553943682 -0.0187102555 -0.0712702861  0.0536858396  0.0218735515  0.0426459641  0.0281934993  0.0339106712 -0.0691240424 
# Bace1       Ankrd54        Pcdh17      BC051077 9530034E10Rik         Calb2         Tmtc2          Fmn1          Mtx2       Tmem159        Ctnna1         Mlxip 
# 0.0518345723 -0.0263834202  0.0316978262  0.0519272158 -0.0238234506 -0.1659305078 -0.0584975323  0.1072732283 -0.0063064823 -0.0729286360 -0.0096570437  0.0062025770 
# Bri3bp       Plekha1          Nsmf       Ciapin1         Dhrs3        Lpgat1         Bcar3       Map3k19          Mgmt          Sv2b         Ddit3       Kirrel3 
# 0.0301189189  0.0278814016  0.0219730267 -0.0164405570 -0.0231118049  0.0149633527 -0.0128249861  0.0274110583  0.0002865444 -0.0178295504  0.0026755173 -0.0150409783 
# Erh        Pcdh10          Tep1       Cyp26b1        Arxes1          Dio2         Susd4          Nfia        Necab3      C1qtnf12          <NA> 
#   0.0072379125 -0.0104767449  0.0159896610 -0.0170622017 -0.0113037675  0.0461949195 -0.0005205602 -0.0870759875 -0.0662154990 -0.0521596487 -2.7392383604 

# <NA> Linear SVM 的截距/bias 项

#### 12. 预测测试集 ####
final_svm <- readRDS("rdata/final_svm.rds")
pred_test <- predict(final_svm, x_test_sc, decisionValues = TRUE)

pred_label <- pred_test$predictions

conf_mat <- table(Pred = pred_label, True = y_test)
conf_mat
# True
# Pred    -1     1
# -1 11746    17
# 1     11   131

##### plot #####
conf_df <- as.data.frame(conf_mat)

conf_df_prop <- conf_df %>%
  group_by(True) %>%
  mutate(prop = Freq / sum(Freq),
         label = paste0(Freq, "\n", round(prop * 100, 1), "%"))

pdf("plot/SVM_prediction_WT_test_mOSN.pdf")
ggplot(conf_df_prop, aes(x = True, y = Pred, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_classic(base_size = 14) +
  labs(
    x = "True label",
    y = "Predicted label",
    fill = "Fraction",
    title = "Confusion matrix on test set"
  )
dev.off()

##### ROC ####
## 用训练集训练好的 final_svm 预测测试集
pred_test <- predict(final_svm, x_test_sc, decisionValues = TRUE)

## 取连续分数
svm_score <- if (is.matrix(pred_test$decisionValues)) {
  pred_test$decisionValues[, 1]
} else {
  pred_test$decisionValues
}

## 把标签改成 0/1
## 假设 TAAR = 1, OR = -1
y_test01 <- ifelse(y_test == 1, 1, 0)

## 画 ROC 并计算 AUC
roc_obj <- roc(response = y_test01, predictor = svm_score)

tail(data.frame(
  threshold = roc_obj$thresholds,
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities,
  fpr = 1 - roc_obj$specificities
), 20)

co <- coords(
  roc_obj,
  x = "all",
  ret = c("specificity", "sensitivity"),
  transpose = FALSE
)

fpr <- 1 - co$specificity
tpr <- co$sensitivity

pdf("plot/SVM_ROC_curve.pdf")
plot(
  fpr, tpr,
  type = "l",
  xlim = c(-0.1, 1.1),
  ylim = c(-0.1, 1.1),
  xaxs = "i",
  yaxs = "i",
  xlab = "False positive rate",
  ylab = "True positive rate",
  main = "Linear SVM ROC curve",
  bty = "l"
)

abline(a = 0, b = 1, lty = 2, col = "grey60")

auc_value <- roc_obj$auc
auc_value

text(
  x = 0.75, y = 0.1,
  labels = paste0("AUC = ", round(as.numeric(auc_value), 3)),
  cex = 1.2
)

dev.off()

##### metrics #####
TP <- conf_mat["1", "1"]
TN <- conf_mat["-1", "-1"]
FP <- conf_mat["1", "-1"]
FN <- conf_mat["-1", "1"]

recall <- TP / (TP + FN)

precision <- TP / (TP + FP)

accuracy <- (TP + TN) / sum(conf_mat)

f1 <- 2 * precision * recall / (precision + recall)

c(accuracy = accuracy,
  recall = recall,
  precision = precision,
  f1 = f1)
# accuracy    recall precision        f1 
# 0.9976480 0.8851351 0.9225352 0.9034483 

metric_df <- data.frame(
  metric = c("Accuracy", "Precision", "Recall", "F1"),
  value  = c(accuracy, precision, recall, f1)
)

ggplot(metric_df, aes(x = metric, y = value)) +
  geom_col() +
  ylim(0, 1) +
  theme_bw() +
  labs(title = "SVM Performance Metrics", y = "Value")

##### 取 decision values ####
# Linear SVM 默认最自然的输出不是概率，而是到分离超平面的距离/决策值。
# e1071 文档说明分类是基于决策函数的符号；对线性核，还可以提取超平面系数。

# decision value > 0 → 模型判为正类，也就是更像 TAAR
# decision value < 0 → 模型判为负类，也就是更像 OR
# 绝对值越大 → 离边界越远，分类越“明确”
# 接近 0 → 靠近边界，说明这个样本更模糊、更不确定

svm_score <- as.vector(pred_test$decisionValues[,1])
head(svm_score)
length(svm_score)

score_df <- data.frame(
  score = svm_score,
  true_class = factor(y_test, levels = c(-1, 1), labels = c("OR", "TAAR"))
)

# 统计每组样本量
n_df <- score_df %>%
  count(true_class) %>%
  mutate(label = paste0(true_class, "\n(n=", n, ")"))

score_df <- score_df %>%
  left_join(n_df[, c("true_class", "label")], by = "true_class")

pdf("plot/SVM_prediction_WT_test_mOSN_decisioin_value.pdf")
ggplot(score_df, aes(x = label, y = score, fill = true_class)) +
  geom_violin(trim = FALSE, width = 0.9, alpha = 0.6, colour = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, colour = "black") +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  geom_hline(yintercept = 0, linetype = 2, colour = "red", linewidth = 0.7) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Decision Value Distribution by True Class",
    y = "SVM Decision Value"
  )

dev.off()

#### 13. 评估 EGFP/Cre+ mOSN ####
load("../20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/seurat_V5_analysis/rdata/EGFP_Cre_mOSN.rdata")

gene_info <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T)
all_functional_olfactory_receptors <- gene_info[which(grepl("Taar|Olfr", gene_info$gene_name) & gene_info$gene_type == "protein_coding" & gene_info$annotation_source == "HAVANA"),]
dim(all_functional_olfactory_receptors)
# 1152   14
TAAR_info <- all_functional_olfactory_receptors[grep("Taar", all_functional_olfactory_receptors$gene_name),]
Taar_name <- TAAR_info$gene_name
# 15

ClassI_ClassII_OR_info <- read.csv("/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
rownames(ClassI_ClassII_OR_info) <- gsub("\\.ps", "-ps",  rownames(ClassI_ClassII_OR_info))

ClassI_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_I" & ClassI_ClassII_OR_info$gene_type == "protein_coding")]
length(ClassI_OR_name)
# 131

ClassII_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_II" & ClassI_ClassII_OR_info$gene_type == "protein_coding")]
length(ClassII_OR_name)
# 1006

get_single_receptor_mOSN <- function(object, group_name) {
  receptor_genes <- intersect(all_functional_olfactory_receptors$gene_name, rownames(object))
  receptor_counts <- FetchData(object, vars = receptor_genes, assay = "SCT", layer = "counts")
  
  object$receptor_number <- rowSums(receptor_counts > 1)
  single_receptor_obj <- subset(object, subset = receptor_number == 1)
  
  Idents(single_receptor_obj) <- "orig.ident"
  obj_subset <- subset(single_receptor_obj, idents = group_name)
  receptor_counts_subset <- FetchData(obj_subset, vars = receptor_genes, assay = "SCT", layer = "counts")
  
  obj_subset$receptor_name <- apply(receptor_counts_subset, 1, function(x) {
    colnames(receptor_counts_subset)[x > 1]
  })
  
  obj_subset$receptor_class <- case_when(
    obj_subset$receptor_name %in% Taar_name ~ "Taar",
    obj_subset$receptor_name %in% ClassI_OR_name ~ "Class_I_OR",
    obj_subset$receptor_name %in% ClassII_OR_name ~ "Class_II_OR",
    TRUE ~ "None"
  )
  
  return(obj_subset)
}

single_mOSN <- get_single_receptor_mOSN(EGFP_Cre_mOSN, c("Het-1", "Homo"))

final_svm <- readRDS(file = "rdata/final_svm.rds")

predict_svm_seurat <- function(seu, feature_genes, train_center, train_scale, model, assay = "SCT", layer = "data") {
  
  # 取表达矩阵：gene x cell
  expr <- GetAssayData(seu, assay = assay, layer = layer)
  
  # 只保留训练用基因
  common_genes <- intersect(feature_genes, rownames(expr))
  missing_genes <- setdiff(feature_genes, rownames(expr))
  
  # 先建一个完整矩阵，缺失基因补0
  mat <- matrix(
    0,
    nrow = ncol(expr),
    ncol = length(feature_genes),
    dimnames = list(colnames(expr), feature_genes)
  )
  
  # 填入已有基因
  mat[, common_genes] <- t(as.matrix(expr[common_genes, , drop = FALSE]))
  
  # 按训练集参数标准化
  mat_sc <- scale(mat, center = train_center[feature_genes], scale = train_scale[feature_genes])
  
  # 防止某些scale后出现NA
  mat_sc[is.na(mat_sc)] <- 0
  
  # 预测
  pred <- predict(model, mat_sc, decisionValues = TRUE)
  
  out <- data.frame(
    cell = rownames(mat_sc),
    svm_pred = pred$predictions,
    svm_score = as.vector(pred$decisionValues),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

seu_list <- SplitObject(single_mOSN, split.by = "orig.ident")
names(seu_list)

pred_list <- lapply(names(seu_list), function(id) {
  seu_sub <- seu_list[[id]]
  
  pred_df <- predict_svm_seurat(
    seu = seu_sub,
    feature_genes = features_exclude_Tbr1,
    train_center = train_center,
    train_scale = train_scale,
    model = final_svm,
    assay = "SCT",
    layer = "data"
  )
  
  pred_df$orig.ident <- id
  pred_df
})

pred_all <- do.call(rbind, pred_list)
head(pred_all)

single_mOSN$svm_pred <- pred_all$svm_pred[match(colnames(single_mOSN), pred_all$cell)]
single_mOSN$svm_score <- pred_all$svm_score[match(colnames(single_mOSN), pred_all$cell)]

single_mOSN$svm_identity <- ifelse(single_mOSN$svm_pred == 1, "TAAR_identity", "OR_identity")
table(single_mOSN$svm_identity, useNA = "ifany")

table(single_mOSN$orig.ident, single_mOSN$svm_identity)
# OR_identity TAAR_identity
# Het-1          23            78
# Homo           62             2
write.csv(table(single_mOSN$orig.ident, single_mOSN$svm_identity), file = "files/linear_SVM_EGFP_Cre_mOSN_identity.csv")
table(single_mOSN$svm_identity, single_mOSN$receptor_class, single_mOSN$orig.ident)
# , ,  = Het-1
# 
# 
# Class_I_OR Class_II_OR Taar
# OR_identity            7           9    7
# TAAR_identity          4           2   72
# 
# , ,  = Homo
# 
# 
# Class_I_OR Class_II_OR Taar
# OR_identity           44          18    0
# TAAR_identity          2           0    0

df <- data.frame(
  cell = colnames(single_mOSN),
  group = single_mOSN$orig.ident,           # 改成你的实际分组列名
  Taar_OR_class = single_mOSN$receptor_class,
  TAAR_OR_pred = single_mOSN$svm_identity
)

pdf("plot/SVM_prediction_EGFP_Cre_mOSN.pdf")
ggplot(df, aes(x = group, fill = TAAR_OR_pred)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Fraction of cells",
    fill = "Predicted identity",
    title = "Predicted identities in EGFP and Cre mOSNs"
  )
dev.off()

