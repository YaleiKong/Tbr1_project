################################################################################
#### Joint integration and transcriptional-distance analysis of WT and Tbr1 mOSNs ####
####
#### This script combines WT single-receptor mOSNs with Tbr1 Het and Homo mOSNs,
#### performs SCT/CCA integration, transfers receptor annotations, and compares
#### the resulting cell populations in PCA space. Run sections sequentially in
#### the analysis directory because several inputs and outputs use relative paths.
################################################################################

#### 1. Initialize the R environment and load packages ####
library(Seurat)
library(ggplot2)
library(grid)
library(ggpubr)
library(stringr)
library(dplyr)
library(tidyr)
library(corpcor)


rm(list=ls())

#### 2. Reconstruct the WT single-receptor mOSN count object ####
## Re-read the six WT count matrices, then retain cells previously identified
## as mature OSNs expressing exactly one receptor.
GSE173947_samples <- paste("GSE173947_homecage",1:6,sep="_")
GSE173947_scRNA_dirs <- str_c("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_10X_scRNA_Seq_seurat_analysis/GSE173947_SandeepRobertDatta_Cell_2021/raw_data/", GSE173947_samples)

WT_scRNA_seq <- lapply(1:length(GSE173947_scRNA_dirs),function(i){
  CreateSeuratObject(counts = Read10X(GSE173947_scRNA_dirs[i]),
                     project = GSE173947_samples[i],assay = "RNA",
                     min.cells = 0,
                     min.features = 0) 
})
for ( i in 1:length(WT_scRNA_seq)){
  WT_scRNA_seq[[i]][["percent.mt"]] = PercentageFeatureSet(WT_scRNA_seq[[i]], pattern="^mt-")
}

merged_WT_scRNA_seq <- merge(
  x = WT_scRNA_seq[[1]],
  y = WT_scRNA_seq[2:6],
  add.cell.ids = GSE173947_samples
)

load("integrated_analysis_seuratV5_Datta_data/rdata/mOSN_scRNA_single_receptor.rdata")
mOSN_scRNA_single_receptor_names <- colnames(mOSN_scRNA_single_receptor)
WT_mOSN_scRNA_seq <- subset(merged_WT_scRNA_seq, cells = mOSN_scRNA_single_receptor_names)
# An object of class Seurat 
# 32285 features across 44745 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 6 layers present: counts.GSE173947_homecage_1, counts.GSE173947_homecage_2, counts.GSE173947_homecage_3, counts.GSE173947_homecage_4, counts.GSE173947_homecage_5, counts.GSE173947_homecage_6

#### 3. Reconstruct the Tbr1 Het and Homo mOSN count object ####
## Re-read the two Cell Ranger matrices and retain cells annotated as mOSNs in
## the previously saved Tbr1 Seurat object.
Tbr1_samples <- c("Het-1","Homo")
Tbr1_scRNA_dirs <- str_c('/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/cellranger_V7.0.1_output_renew_ref_exclude_intron/', Tbr1_samples,  "/outs/filtered_feature_bc_matrix/")
Tbr1_scRNA_seq <- lapply(1:length(Tbr1_scRNA_dirs),function(i){
  CreateSeuratObject(counts = Read10X(Tbr1_scRNA_dirs[i]),
                     project = Tbr1_samples[i], assay = "RNA",
                     min.cells = 0,
                     min.features = 0) 
})

for ( i in 1:length(Tbr1_scRNA_seq)){
  Tbr1_scRNA_seq[[i]][["percent.mt"]] = PercentageFeatureSet(Tbr1_scRNA_seq[[i]], pattern="^mt-")
}

load("../20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/seurat_V5_analysis/rdata/CCAIntegration_obj_annotated_cluster.rdata")
Idents(CCAIntegration_obj)
Tbr1_mOSN <- subset(
  CCAIntegration_obj,
  idents = "mOSN"
)
Tbr1_mOSN_names <- colnames(
  subset(
    CCAIntegration_obj,
    idents = "mOSN"
  )
)

Tbr1_scRNA_seq_merge_obj <- merge(Tbr1_scRNA_seq[[1]], y = Tbr1_scRNA_seq[[2]])

Tbr1_mOSN_scRNA_seq <- subset(
  Tbr1_scRNA_seq_merge_obj,
  cells = Tbr1_mOSN_names
)
Tbr1_mOSN_scRNA_seq
# An object of class Seurat 
# 32287 features across 1544 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 2 layers present: counts.Het-1, counts.Homo

#### 4. Integrate WT and Tbr1 mOSNs ####
## SCTransform is applied to the merged raw-count object before CCA integration.
all_merged_obj <- merge(WT_mOSN_scRNA_seq, Tbr1_mOSN_scRNA_seq)
# An object of class Seurat 
# 32287 features across 46289 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 8 layers present: counts.GSE173947_homecage_1.1, counts.GSE173947_homecage_2.1, counts.GSE173947_homecage_3.1, counts.GSE173947_homecage_4.1, counts.GSE173947_homecage_5.1, counts.GSE173947_homecage_6.1, counts.Het-1.2, counts.Homo.2

all_merged_obj <- SCTransform(all_merged_obj, vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50)

DefaultAssay(all_merged_obj)

CCAIntegration_obj <- IntegrateLayers(
  object = all_merged_obj,
  method = CCAIntegration,
  assay = "SCT",
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated",
  verbose = TRUE
)

CCAIntegration_obj <- RunUMAP(CCAIntegration_obj, dims = 1:30, 
                              reduction = "integrated", reduction.name = "umap")
save(CCAIntegration_obj, file = "integrated_analysis_seuratV5_Datta_data/Tbr1_integration/CCAIntegration_WT_mOSN_single_receptor_with_Tbr1_mOSN.rdata")

#### 5. Plot the integrated WT and Tbr1-KO mOSNs ####
## Assign a simplified WT/Tbr1-KO label and calculate the small UMAP-axis
## coordinates used in the final integrated UMAP figure.
CCAIntegration_obj$integration_label <- ifelse(CCAIntegration_obj$orig.ident %in% c("Het-1", "Homo"),
                                               "Tbr1-KO", "WT")
Idents(CCAIntegration_obj) <- "integration_label"
## 提取 UMAP 坐标范围
umap_coord <- Embeddings(
  CCAIntegration_obj,
  reduction = "umap"
)

x_range <- range(umap_coord[, 1])
y_range <- range(umap_coord[, 2])

## 左下角坐标轴起点
x0 <- x_range[1] + 0.05 * diff(x_range)
y0 <- y_range[1] + 0.05 * diff(y_range)

## 坐标轴长度
x_len <- 0.18 * diff(x_range)
y_len <- 0.18 * diff(y_range)
pdf("integrated_analysis_seuratV5_Datta_data/Tbr1_integration/umap_plot_WT_Tbr1_KO_integration_red.pdf", width = 6, height = 5)
DimPlot(
  CCAIntegration_obj,
  group.by = "integration_label",
  reduction = "umap",
  cols = c("#fd8d3c", '#fff5eb'),
  pt.size = 0.3
) +
  ## 去掉原始坐标轴
  theme_void(base_size = 12) +
  
  ## UMAP 1 横轴箭头
  annotate(
    "segment",
    x = x0,
    xend = x0 + x_len,
    y = y0,
    yend = y0,
    linewidth = 0.6,
    arrow = arrow(
      length = unit(0.12, "cm"),
      type = "closed"
    )
  ) +
  
  ## UMAP 2 纵轴箭头
  annotate(
    "segment",
    x = x0,
    xend = x0,
    y = y0,
    yend = y0 + y_len,
    linewidth = 0.6,
    arrow = arrow(
      length = unit(0.12, "cm"),
      type = "closed"
    )
  ) +
  
  ## 横轴标签
  annotate(
    "text",
    x = x0 + x_len / 2,
    y = y0 - 0.04 * diff(y_range),
    label = "UMAP 1",
    size = 3.5
  ) +
  
  ## 纵轴标签
  annotate(
    "text",
    x = x0 - 0.04 * diff(x_range),
    y = y0 + y_len / 2,
    label = "UMAP 2",
    angle = 90,
    size = 3.5
  ) +
  theme(
    aspect.ratio = 1,
    legend.position = "right"
  )
#### 6. Transfer WT receptor annotations to the integrated object ####
## Cell names are matched explicitly so annotations remain aligned after merging.
CCAIntegration_obj$receptor <- 'none'
idx <- match(
  colnames(CCAIntegration_obj),
  colnames(mOSN_scRNA_single_receptor)
)

CCAIntegration_obj$receptor[
  !is.na(idx)
] <- mOSN_scRNA_single_receptor$receptor[
  idx[!is.na(idx)]
]

table(CCAIntegration_obj$receptor)

## Recreate the receptor-class reference lists using the same annotation files
## and filtering criteria as in the original Tbr1 scRNA-seq analysis.
gene_info <- read.delim(
  "/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt",
  header = TRUE
)

all_functional_olfactory_receptors <- gene_info[
  grepl("Taar|Olfr", gene_info$gene_name) &
    gene_info$gene_type == "protein_coding" &
    gene_info$annotation_source == "HAVANA",
]

TAAR_info <- all_functional_olfactory_receptors[
  grep("Taar", all_functional_olfactory_receptors$gene_name),
]
Taar_name <- TAAR_info$gene_name

ClassI_ClassII_OR_info <- read.csv(
  "/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv",
  header = TRUE,
  row.names = 1
)
rownames(ClassI_ClassII_OR_info) <- gsub(
  "\\.ps",
  "-ps",
  rownames(ClassI_ClassII_OR_info)
)

ClassI_OR_name <- rownames(ClassI_ClassII_OR_info)[
  ClassI_ClassII_OR_info$OR_Class == "Class_I" &
    ClassI_ClassII_OR_info$gene_type == "protein_coding"
]

ClassII_OR_name <- rownames(ClassI_ClassII_OR_info)[
  ClassI_ClassII_OR_info$OR_Class == "Class_II" &
    ClassI_ClassII_OR_info$gene_type == "protein_coding"
]

CCAIntegration_obj$receptor_class <- ifelse(CCAIntegration_obj$receptor %in% ClassI_OR_name, "Class_I_OR", 
                                            ifelse(CCAIntegration_obj$receptor %in% ClassII_OR_name, "Class_II_OR",
                                                   ifelse(CCAIntegration_obj$receptor %in% Taar_name, "Taar",
                                                          "none")))

EGFP_Cre_mOSN <- subset(CCAIntegration_obj, subset = EGFP > 0 | Cre > 0)
#### 7. Prepare the cell populations for transcriptional-distance analysis ####
## Extract the first 30 PCs and select the three cell populations using the
## same definitions as in the original analysis.
pca_embedding <- Embeddings(
  CCAIntegration_obj,
  reduction = "pca"
)[, 1:30]

Taar_cells <- WhichCells(
  CCAIntegration_obj,
  expression = receptor_class == "Taar"
)
Taar_pca <- pca_embedding[Taar_cells, , drop = FALSE]

Homo_EGFP_Cre_cells <- WhichCells(
  EGFP_Cre_mOSN,
  expression = orig.ident == "Homo"
)
Homo_pca <- pca_embedding[Homo_EGFP_Cre_cells, , drop = FALSE]

Het_EGFP_Cre_cells <- WhichCells(
  EGFP_Cre_mOSN,
  expression = orig.ident == "Het-1"
)
Het_pca <- pca_embedding[Het_EGFP_Cre_cells, , drop = FALSE]
#### 8. Visualize population positions and transcriptional distances ####
##### 8.1 Integrated UMAP ####
plot_obj <- subset(
  CCAIntegration_obj,
  subset = receptor_class != "none" |
    EGFP > 0 |
    Cre > 0
)
plot_obj$label <- ifelse(plot_obj$receptor_class=="Taar", "WT TAAR OSN", 
                         ifelse(plot_obj$receptor_class %in% c("Class_I_OR", "Class_II_OR"), "WT OR OSN",
                                ifelse(plot_obj$orig.ident == "Het-1", "Tbr1-KO Het EGFP/Cre mOSN",
                                       ifelse(plot_obj$orig.ident == "Homo", "Tbr1-KO Homo EGFP/Cre mOSN",
                                              NA))))
plot_obj$label <- factor(plot_obj$label, levels = c("WT OR OSN", "WT TAAR OSN", "Tbr1-KO Het EGFP/Cre mOSN", "Tbr1-KO Homo EGFP/Cre mOSN"))
## 固定分组及图例顺序
plot_obj$label <- factor(
  plot_obj$label,
  levels = c(
    "WT OR OSN",
    "WT TAAR OSN",
    "Tbr1-KO Het EGFP/Cre mOSN",
    "Tbr1-KO Homo EGFP/Cre mOSN"
  )
)

## 提取 UMAP 坐标范围
umap_coord <- Embeddings(
  plot_obj,
  reduction = "umap"
)

x_range <- range(umap_coord[, 1], na.rm = TRUE)
y_range <- range(umap_coord[, 2], na.rm = TRUE)

x_span <- diff(x_range)
y_span <- diff(y_range)

## 左下角小坐标轴的位置
x0 <- x_range[1] + 0.06 * x_span
y0 <- y_range[1] + 0.06 * y_span

## 小坐标轴长度
x_len <- 0.16 * x_span
y_len <- 0.16 * y_span

p <- DimPlot(
  plot_obj,
  group.by = "label",
  reduction = "umap",
  cols = c(
    "WT TAAR OSN" = "#4292C6",
    "WT OR OSN" = "#C6DBEF",
    "Tbr1-KO Het EGFP/Cre mOSN" = "#D9F0A3",
    "Tbr1-KO Homo EGFP/Cre mOSN" = "#238443"
  ),
  pt.size = 0.5,
  order = TRUE
) +
  ## 去掉完整坐标轴、刻度和背景
  theme_void(base_size = 10) +
  
  ## UMAP 1 横向箭头
  annotate(
    geom = "segment",
    x = x0,
    xend = x0 + x_len,
    y = y0,
    yend = y0,
    linewidth = 0.55,
    colour = "black",
    arrow = grid::arrow(
      length = unit(0.10, "cm"),
      type = "closed"
    )
  ) +
  
  ## UMAP 2 纵向箭头
  annotate(
    geom = "segment",
    x = x0,
    xend = x0,
    y = y0,
    yend = y0 + y_len,
    linewidth = 0.55,
    colour = "black",
    arrow = grid::arrow(
      length = unit(0.10, "cm"),
      type = "closed"
    )
  ) +
  
  ## UMAP 1 标签
  annotate(
    geom = "text",
    x = x0 + x_len / 2,
    y = y0 - 0.045 * y_span,
    label = "UMAP 1",
    size = 3.2
  ) +
  
  ## UMAP 2 标签
  annotate(
    geom = "text",
    x = x0 - 0.04 * x_span,
    y = y0 + y_len / 2,
    label = "UMAP 2",
    angle = 90,
    size = 3.2
  ) +
  theme(
    ## 主绘图区为正方形
    aspect.ratio = 1,
    
    legend.position = "right",
    legend.title = element_blank(),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

pdf("integrated_analysis_seuratV5_Datta_data/Tbr1_integration/umap_TAAR_OR_Tbr1_Het_Homo_EGFP.pdf", width = 7,height = 5)
p
dev.off()
#### 9. Compare Euclidean and shrinkage Mahalanobis distances ####

##### 9.1 Extract the first 30 PCs ####


n_pc <- 30

Taar_mat <- as.matrix(
  Taar_pca[, seq_len(n_pc), drop = FALSE]
)

Het_mat <- as.matrix(
  Het_pca[, seq_len(n_pc), drop = FALSE]
)

Homo_mat <- as.matrix(
  Homo_pca[, seq_len(n_pc), drop = FALSE]
)



##### 9.2 Define the WT TAAR reference state ####


## WT TAAR 中心
Taar_center <- colMeans(
  Taar_mat,
  na.rm = TRUE
)

## 原始协方差矩阵，仅用于比较条件数
Taar_cov <- cov(
  Taar_mat,
  use = "pairwise.complete.obs"
)

raw_kappa <- kappa(Taar_cov)

message(
  "Condition number of raw covariance matrix: ",
  round(raw_kappa, 3)
)

##### 9.3 Estimate a shrinkage covariance matrix ####


Taar_cov_shrink <- corpcor::cov.shrink(
  Taar_mat,
  verbose = FALSE
)

shrink_kappa <- kappa(Taar_cov_shrink)

message(
  "Condition number of shrinkage covariance matrix: ",
  round(shrink_kappa, 3)
)

## 后续 Mahalanobis distance 使用收缩协方差
Taar_cov_inv_shrink <- solve(
  Taar_cov_shrink
)



##### 9.4 Define distance functions ####


## Euclidean distance
get_euclidean <- function(x, center) {
  
  x <- as.matrix(x)
  
  delta <- sweep(
    x,
    MARGIN = 2,
    STATS = center,
    FUN = "-"
  )
  
  sqrt(
    rowSums(delta^2)
  )
}


## Shrinkage Mahalanobis distance
get_mahal <- function(
    x,
    center,
    cov_inv
) {
  
  x <- as.matrix(x)
  
  delta <- sweep(
    x,
    MARGIN = 2,
    STATS = center,
    FUN = "-"
  )
  
  mahal_sq <- rowSums(
    (delta %*% cov_inv) * delta
  )
  
  ## 避免浮点误差产生极小负数
  mahal_sq[mahal_sq < 0] <- 0
  
  sqrt(mahal_sq)
}

##### 9.5 Calculate Het and Homo distances from the WT TAAR state ####


distance_df <- dplyr::bind_rows(
  
  data.frame(
    group = "Het",
    Euclidean = get_euclidean(
      Het_mat,
      Taar_center
    ),
    Mahalanobis = get_mahal(
      Het_mat,
      Taar_center,
      Taar_cov_inv_shrink
    )
  ),
  
  data.frame(
    group = "Homo",
    Euclidean = get_euclidean(
      Homo_mat,
      Taar_center
    ),
    Mahalanobis = get_mahal(
      Homo_mat,
      Taar_center,
      Taar_cov_inv_shrink
    )
  )
)

distance_df$group <- factor(
  distance_df$group,
  levels = c("Het", "Homo")
)



##### 9.6 Reshape the distance table for plotting ####


plot_df <- distance_df %>%
  pivot_longer(
    cols = c(
      Euclidean,
      Mahalanobis
    ),
    names_to = "metric",
    values_to = "distance"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Euclidean",
        "Mahalanobis"
      ),
      labels = c(
        "Euclidean distance",
        "Shrinkage Mahalanobis distance"
      )
    )
  )


## 检查距离分布
plot_df %>%
  group_by(metric, group) %>%
  summarise(
    n = n(),
    mean = mean(distance),
    median = median(distance),
    sd = sd(distance),
    .groups = "drop"
  )



##### 9.7 Plot the distance distributions ####



p_distance <- ggplot(
  plot_df,
  aes(
    x = group,
    y = distance,
    fill = group,
    colour = group
  )
) +
  geom_violin(
    alpha = 0.35,
    trim = FALSE,
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.6,
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("Het", "Homo")
    ),
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_fill_manual(
    values = c(
      "Het" =
        scales::alpha("#D9F0A3", 0.5),
      "Homo" =
        scales::alpha("#238443", 0.5)
    )
  ) +
  scale_colour_manual(
    values = c(
      "Het" = "#8FCB50",
      "Homo" = "#238443"
    )
  ) +
  labs(
    x = NULL,
    y = "Distance from WT TAAR state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 4 / 3,
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


##### 9.8 Save the comparison plot ####

pdf(
  paste0(
    "integrated_analysis_seuratV5_Datta_data/",
    "Tbr1_integration/",
    "euclidean_shrinkage_mahalanobis_distance_vlnplot.pdf"
  ),
  width = 8,
  height = 4
)

print(p_distance)

dev.off()
