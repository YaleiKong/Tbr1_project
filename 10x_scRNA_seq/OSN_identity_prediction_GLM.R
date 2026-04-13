rm(list=ls())

# install.packages("glmnet")
library(glmnet)
# 对于每个细胞，模型把一堆基因表达值组合起来，先算一个线性分数

# install.packages("pROC")
library(pROC)

# install.packages("PRROC")
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
set.seed(42)

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
table(meta[test_cells, "label"])
# OR  Taar 
# 11757   148 

#### 3. 仅在训练集上计算 DEGs (Taar vs OR) ####
train_obj <- subset(mOSN_scRNA_single_receptor, cells = train_cells)
Idents(train_obj) <- "Taar_OR_class"
train_obj <- PrepSCTFindMarkers(train_obj)

train_markers <- FindMarkers(
  train_obj,
  ident.1 = "Taar",
  ident.2 = "OR",
  logfc.threshold = 0.25,
  test.use = "wilcox",
  only.pos = FALSE
)

## 加 gene 列
train_markers$gene <- rownames(train_markers)

## 筛选：保留显著且 |log2FC| >= 0.585
train_markers_use <- train_markers[
  train_markers$p_val_adj < 0.05 & abs(train_markers$avg_log2FC) >= 0.585,
]

## 去掉受体基因本身，避免模型只靠 Taar / Olfr 直接分类
train_markers_use <- train_markers_use[
  !grepl("Taar|Olfr", train_markers_use$gene),
]

## 保存训练集 DEGs
write.xlsx(
  list(
    train_markers_all = cbind(gene = rownames(train_markers), train_markers),
    train_markers_filtered = train_markers_use
  ),
  file = "files/Taar_OR_train_only_DE_results_mOSN_expressing_single_receptor.xlsx"
)

#### 4-1. 提取表达矩阵 ####
mat <- GetAssayData(mOSN_scRNA_single_receptor, assay = "SCT", layer = "data")

train_markers_use <- read.xlsx("files/Taar_OR_train_only_DE_results_mOSN_expressing_single_receptor.xlsx",
                               sheet = "train_markers_filtered")

features <- train_markers_use$gene

features <- intersect(features, rownames(mat))

features_exclude_Tbr1 <- features[-which(features == "Tbr1")]

#### 5-1. 构造 x ####
# 把表达矩阵变成 cells × genes，符合 glmnet 输入格式
x1_train <- t(as.matrix(mat[features_exclude_Tbr1, train_cells]))
x1_test  <- t(as.matrix(mat[features_exclude_Tbr1, test_cells]))

## 可选：去掉训练集中零方差基因
keep <- apply(x1_train, 2, function(z) sd(z) > 0)
x1_train <- x1_train[, keep, drop = FALSE]
x1_test  <- x1_test[, keep, drop = FALSE]

dim(x1_train)
# 47568   190

#### 6-1. 构造 y 构造二分类标签 ####
y_train <- ifelse(meta[train_cells, "label"] == "Taar", 1, 0)
y_test  <- ifelse(meta[test_cells, "label"] == "Taar", 1, 0)

table(y_train)
# y_train
# 0     1 
# 46992   576 

table(y_test)
# y_test
# 0     1 
# 11757   148 

#### 7-1. 训练二分类模型 ####
cvfit1 <- cv.glmnet(
  x = x1_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)
# x：训练特征矩阵，cells × genes
# y：0/1 标签
# family = "binomial"：二分类逻辑回归
# alpha = 1：LASSO（L1 惩罚）
# nfolds = 5：5 折交叉验证

## 保存模型
saveRDS(cvfit1, "../TAAR_OR_identity_model/rdata/TAAR_OR_glmnet_model_without_Tbr1.rds")
saveRDS(features_exclude_Tbr1, "../TAAR_OR_identity_model/rdata/TAAR_OR_glmnet_features_without_Tbr1.rds")

model_info1 <- list(
  model = cvfit1,
  features = colnames(x1_train),   # 最好直接保存训练矩阵列名
  lambda_use = "lambda.1se",
  assay = "SCT",
  layer = "data"
)

saveRDS(model_info1, "../TAAR_OR_identity_model/rdata/TAAR_OR_glmnet_model_info_without_Tbr1.rds")

## coef()：返回截距和每个基因的系数。
coef_mat1 <- coef(cvfit1, s = "lambda.1se")
# lambda.min：交叉验证误差最小的那个 lambda
# lambda.1se：在“误差不比最小值明显差”的前提下，选择的更大、更保守的 lambd
# s = "lambda.1se"：使用更保守的 lambda。

coef_df1 <- data.frame(
  gene = rownames(coef_mat1),
  coef = as.numeric(coef_mat1)
)

coef_df1_only_coef_not_0 <- coef_df1 %>% 
  filter(coef != 0)
dim(coef_df1_only_coef_not_0)
# 98  2
coef_df1_only_coef_not_0
# gene         coef
# 1    (Intercept) -4.284409111
# 2         Psmd14  2.739501797
# 3        Gucy1b2  0.062163991
# 4         Tmbim6  0.412685820
# 5            Ckb  0.343837699
# 6          Hint1  1.015616888
# 7         Stoml3 -1.575369552
# 8         Tbc1d4  1.101802126
# 9        Sult1d1  0.460623965
# 10         Tafa4  0.384405216
# 11          Nsg2  0.245656089
# 12         Stbd1  0.320004646
# 13        Ndfip1  0.192892475
# 14         Mslnl  0.076076627
# 15         Acss2  0.163186880
# 16         Tenm2  0.075818928
# 17         Rgs16  1.409588981
# 18        Stmnd1  0.015455949
# 19         Myo1b  0.885606525
# 20          Snca  0.027579392
# 21           Ptn  0.201500814
# 22           Cd9 -0.223735332
# 23        Adam28  0.379666573
# 24          Ebf1 -0.372408176
# 25          Rtn1 -0.215412095
# 26         Lsamp  1.090755777
# 27          Ece1  0.156671756
# 28         Fetub -0.078400877
# 29          Ccn3  0.046242095
# 30       Fam129a  0.135732051
# 31          Faim -0.340587153
# 32        Pdlim1 -0.096981458
# 33         Olig1 -0.214937512
# 34        Sorbs2 -0.097647799
# 35         Prrg4  0.499014779
# 36         Cidea -0.319965215
# 37          Glul -0.118808080
# 38          Pth2  0.088862706
# 39          Dmkn -0.007979118
# 40          Bmp6  0.568678132
# 41        Kcnmb3  0.065948617
# 42         Pcdh7 -0.615244787
# 43          Gldc -0.092636229
# 44         Celf4 -0.016489764
# 45         Auts2 -0.065183654
# 46         Rpl11 -0.092846964
# 47        Sema7a -0.109630726
# 48       Mkrn2os  0.625536847
# 49         Tshz1 -0.201893653
# 50         Hjurp -0.205250235
# 51        Plaat3  0.504753373
# 52        Trim45 -0.221158454
# 53       Gm13589 -0.258304923
# 54         Sept3 -0.024578115
# 55          Map4  0.487312718
# 56         Tppp3 -0.097499903
# 57         Rims3 -0.335770861
# 58         Cldn3 -0.005784063
# 59         Ifi27 -0.008377397
# 60       Gm31557 -0.012016555
# 61 4930515B02Rik  0.002151506
# 62       Gm13648 -0.086077306
# 63    AC149090.1 -0.216314546
# 64         Pcdh9 -0.010239611
# 65        Elavl2  0.005847070
# 66         Olig2 -0.267236802
# 67          Rbm4 -0.029496803
# 68         Med27  0.274131186
# 69         Cotl1 -0.174295976
# 70         Gdf11  0.201548532
# 71        Bicdl2  0.117163936
# 72        Bcl11b -0.145992391
# 73        Cox7a1  0.132342937
# 74      Tmem132a -0.340466709
# 75           Cpe -0.106362315
# 76          Caly -0.206779119
# 77       Gm16579 -0.254818588
# 78         Ssbp2 -0.130926345
# 79        Pde10a  0.054538090
# 80      Cacna2d3 -0.172836038
# 81       Phactr1 -0.051275563
# 82        Cab39l -0.150676346
# 83         Mef2b -0.024679792
# 84          Nqo1 -0.780607003
# 85      Pafah1b3 -0.143654035
# 86         Prpf4 -0.168917509
# 87         Kcnh7  0.570656636
# 88       Ppp2r3a  0.007329422
# 89         Clip4 -0.129200138
# 90       Ankrd54 -0.085102012
# 91         Calb2 -0.318984830
# 92         Tmtc2 -0.206408291
# 93          Fmn1  1.043695517
# 94       Tmem159 -0.235631209
# 95        Lpgat1  0.071686989
# 96          Nfia -0.178217029
# 97        Necab3 -0.056927643
# 98      C1qtnf12 -0.063866411
write.csv(coef_df1_only_coef_not_0, file = "files/GLM_model_genes_weight.csv")

#### 8-1. 预测测试集 ####
cvfit1 <- readRDS("rdata/TAAR_OR_glmnet_model_without_Tbr1.rds")
features_exclude_Tbr1 <- readRDS("rdata/TAAR_OR_glmnet_features_without_Tbr1.rds")

prob_test1 <- as.numeric(
  predict(cvfit1, newx = x1_test, s = "lambda.1se", type = "response")
)
# newx = x_test：测试集矩阵
# s = "lambda.1se"：使用保守 lambda
# type = "response"：返回概率而不是线性预测值

# lambda.1se
# GSE120199_P0_rep1_AACTCCCAGATCCCGC-1 2.671058e-04
# GSE120199_P0_rep1_AAGTCTGGTTTCGCTC-1 1.750924e-04

#### 9-1. 评估 测试组 ####

## 测试不同threshold的
ths <- seq(0.1, 0.9, by = 0.05)

metrics <- lapply(ths, function(th) {
  pred <- ifelse(prob_test1 > th, 1, 0)
  tp <- sum(pred == 1 & y_test == 1)
  tn <- sum(pred == 0 & y_test == 0)
  fp <- sum(pred == 1 & y_test == 0)
  fn <- sum(pred == 0 & y_test == 1)
  
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  specificity <- ifelse(tn + fp == 0, NA, tn / (tn + fp))
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall == 0),
               NA,
               2 * precision * recall / (precision + recall))
  
  data.frame(
    threshold = th,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1 = f1
  )
})

metrics_df <- do.call(rbind, metrics)
# Precision “我捞出来的这些细胞，有多少是真的？”
# Recall “所有真的 TAAR 细胞里，我抓到了多少？”
# Specificity “所有 OR 细胞里，我成功排除了多少？”
# F1 “precision 和 recall 总体平衡得怎么样？”
# Accuracy “全部细胞里，总共判对了多少？”

metrics_df
# threshold precision    recall specificity        f1
# 1       0.10 0.7329843 0.9459459   0.9956622 0.8259587
# 2       0.15 0.8242424 0.9189189   0.9975334 0.8690096
# 3       0.20 0.8645161 0.9054054   0.9982138 0.8844884
# 4       0.25 0.8926174 0.8986486   0.9986391 0.8956229
# 5       0.30 0.9084507 0.8716216   0.9988943 0.8896552
# 6       0.35 0.9275362 0.8648649   0.9991494 0.8951049
# 7       0.40 0.9338235 0.8581081   0.9992345 0.8943662
# 8       0.45 0.9407407 0.8581081   0.9993196 0.8975265
# 9       0.50 0.9398496 0.8445946   0.9993196 0.8896797
# 10      0.55 0.9520000 0.8040541   0.9994897 0.8717949
# 11      0.60 0.9500000 0.7702703   0.9994897 0.8507463
# 12      0.65 0.9565217 0.7432432   0.9995747 0.8365019
# 13      0.70 0.9553571 0.7229730   0.9995747 0.8230769
# 14      0.75 0.9626168 0.6959459   0.9996598 0.8078431
# 15      0.80 0.9705882 0.6689189   0.9997448 0.7920000
# 16      0.85 0.9791667 0.6351351   0.9998299 0.7704918
# 17      0.90 0.9756098 0.5405405   0.9998299 0.6956522

# 阈值是基于内部测试集调出来的；
# 再固定应用到外部数据。
# 根据f1“整体最平衡”，选择0.45
metrics_df$threshold[which(metrics_df$f1 == max(metrics_df$f1))]
# 0.45
acc_df <- do.call(rbind, lapply(ths, function(th) {
  pred <- ifelse(prob_test1 > th, 1, 0)
  acc <- mean(pred == y_test)
  
  data.frame(
    threshold = th,
    accuracy = acc
  )
}))

acc_df
# threshold  accuracy
# 1       0.10 0.9950441
# 2       0.15 0.9965561
# 3       0.20 0.9970601
# 4       0.25 0.9973961
# 5       0.30 0.9973121
# 6       0.35 0.9974801
# 7       0.40 0.9974801
# 8       0.45 0.9975640
# 9       0.50 0.9973961
# 10      0.55 0.9970601
# 11      0.60 0.9966401
# 12      0.65 0.9963881
# 13      0.70 0.9961361
# 14      0.75 0.9958841
# 15      0.80 0.9956321
# 16      0.85 0.9952961
# 17      0.90 0.9941201

# 阈值设置为0.45，准确率也较高
pred_test1 <- ifelse(prob_test1 > 0.45, 1, 0)

conf_mat1 <- table(Pred = pred_test1, True = y_test)
print(conf_mat1)
# True
# Pred     0     1
# 0 11749    21
# 1     8   127

acc1 <- mean(pred_test1 == y_test)

print(acc1)
# 0.997564

#### 9-1 评估 测试组 plot ####

##### threshold–metrics 曲线 #####

plot_df <- metrics_df %>%
  select(threshold, precision, recall, specificity, f1) %>%
  pivot_longer(
    cols = -threshold,
    names_to = "metric",
    values_to = "value"
  )

pdf("plot/threshold_metrics.pdf")
ggplot(plot_df, aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_classic(base_size = 14) +
  labs(
    x = "Threshold",
    y = "Metric value",
    title = "Model performance across thresholds"
  )
dev.off()

##### ROC 曲线 #####

roc_obj <- roc(y_test, prob_test1)

# AUC 看的是：如果随机抽一个真实 Taar 细胞和一个真实 OR 细胞，模型把 Taar 的分数排在 OR 前面的概率。
auc_val <- auc(roc_obj)
auc_num <- as.numeric(auc_val)
auc_num

roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

pdf("plot/roc_curve.pdf")
ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    title = paste0("ROC curve (AUC = ", round(auc_val, 3), ")")
  )
dev.off()

##### PR 曲线 #####

pr <- pr.curve(
  scores.class0 = prob_test1[y_test == 1],
  scores.class1 = prob_test1[y_test == 0],
  curve = TRUE
)
round(pr$auc.integral, 3)

pr_df <- data.frame(
  recall = pr$curve[, 1],
  precision = pr$curve[, 2]
)

pdf("plot/pr_curve.pdf")
ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 14) +
  labs(
    x = "Recall",
    y = "Precision",
    title = paste0("PR curve (AUPRC = ", round(pr$auc.integral, 3), ")")
  )
# 横轴 = Recall
# 纵轴 = Precision
dev.off()

# Taar 是少数类
# OR 是多数类
# 
# 这种场景下，AUPRC 比 ROC-AUC 往往更有信息量。
# 
# 因为 PR 曲线直接关注的是：
# 
# 我找回了多少 Taar（recall）
# 我预测出来的 Taar 有多真（precision）
# AUPRC 高 说明模型好：  找回更多正类时，precision 还能保持高 预测成 Taar 的细胞里，假阳性少
# AUPRC 低 说明模型差： 想提高 recall，就得牺牲很多 precision 或者模型抓出来的“正类”里混了很多假阳性

# 模型在识别少数类（Taar）时，能否在较高召回的同时保持较高精度。

##### heatmap (WT 测试组) #####
conf_df1 <- as.data.frame(conf_mat1)

conf_df1_prop <- conf_df1 %>%
  group_by(True) %>%
  mutate(prop = Freq / sum(Freq),
         label = paste0(Freq, "\n", round(prop * 100, 1), "%"))

pdf("plot/prediction_WT_test_mOSN.pdf")
ggplot(conf_df1_prop, aes(x = True, y = Pred, fill = prop)) +
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

#### 10-1 评估 EGFP/Cre+ mOSN ####
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

model_info1 <- readRDS("../TAAR_OR_identity_model/rdata/TAAR_OR_glmnet_model_info_without_Tbr1.rds")

predict_taar_or_glmnet <- function(seu, model_info, threshold = 0.5) {
  cvfit <- model_info$model
  features_train <- model_info$features
  lambda_use <- model_info$lambda_use
  assay_use <- model_info$assay
  layer_use <- model_info$layer
  
  ## 提取表达矩阵
  mat <- GetAssayData(seu, assay = assay_use, layer = layer_use)
  
  ## 共同基因
  common_genes <- intersect(features_train, rownames(mat))
  missing_genes <- setdiff(features_train, rownames(mat))
  
  message("Matched genes: ", length(common_genes), "/", length(features_train))
  message("Missing genes: ", length(missing_genes))
  
  ## 构造 newx，缺失基因补 0
  x_new <- matrix(
    0,
    nrow = ncol(mat),
    ncol = length(features_train),
    dimnames = list(colnames(mat), features_train)
  )
  
  x_new[, common_genes] <- t(as.matrix(mat[common_genes, ]))
  
  x_new <- as.matrix(x_new)
  storage.mode(x_new) <- "double"
  
  ## 预测 TAAR 概率
  prob <- as.numeric(
    predict(cvfit, newx = x_new, s = lambda_use, type = "response")
  )
  
  pred <- ifelse(prob > threshold, "Taar", "OR")
  
  ## 写回对象
  seu$TAAR_score <- prob
  seu$TAAR_OR_pred <- pred
  
  return(seu)
}
# seu：新对象
# model_info：保存的模型与元信息
# threshold：分类阈值
# common_genes：训练和新数据共有基因
# missing_genes：新数据缺失的训练特征
# x_new：按训练特征顺序构造的新数据矩阵，缺失基因补 0

# threshold = 0.45
obj <- predict_taar_or_glmnet(single_mOSN, model_info1, threshold = 0.45)

df <- data.frame(
  cell = colnames(obj),
  group = obj$orig.ident,           # 改成你的实际分组列名
  Taar_OR_class = obj$receptor_class,
  TAAR_OR_pred = obj$TAAR_OR_pred,
  TAAR_score = obj$TAAR_score
)

table(df$Taar_OR_class, df$group)
# Het-1 Homo
# Class_I_OR     11   46
# Class_II_OR    11   18
# Taar           79    0

write.csv(table(df$TAAR_OR_pred, df$group), file = "files/GLM_EGFP_Cre_mOSN_identity.csv")

table(df$TAAR_OR_pred, df$Taar_OR_class, df$group)
# , ,  = Het-1
# 
# 
# Class_I_OR Class_II_OR Taar
# OR            6           9    5
# Taar          5           2   74
# 
# , ,  = Homo
# 
# 
# Class_I_OR Class_II_OR Taar
# OR           43          18    0
# Taar          3           0    0

pdf("plot/prediction_EGFP_Cre_mOSN.pdf")
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



