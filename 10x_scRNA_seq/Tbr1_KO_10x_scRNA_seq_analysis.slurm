################################################################################################################
# cellranger_count_Het_1_exclude_intron.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input_dir=/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/raw_data/
output_dir=/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/cellranger_V7.0.1_output_renew_ref_exclude_intron/
reference_dir=/lustre/home/acct-medlqian/medlqian-loop3/database/Mus_musculus/GRCm38_mm10/cellranger/refdata_with_EGFP_Cre_build_cellrangerV7.0.1_20260221/Mouse_EGFP_Cre/

## generate single cell feature counts by cellranger count
cd ${output_dir}

TENX_IGNORE_DEPRECATED_OS=1  cellranger count --id=Het-1 \
--fastqs=${input_dir} \
--sample=Het-1 \
--transcriptome=${reference_dir} \
--include-introns false \
--localcores=10 --localmem=40 

################################################################################################################
# cellranger_count_Homo_exclude_intron.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

input_dir=/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/raw_data/
output_dir=/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/cellranger_V7.0.1_output_renew_ref_exclude_intron/
reference_dir=/lustre/home/acct-medlqian/medlqian-loop3/database/Mus_musculus/GRCm38_mm10/cellranger/refdata_with_EGFP_Cre_build_cellrangerV7.0.1_20260221/Mouse_EGFP_Cre/

## generate single cell feature counts by cellranger count
cd ${output_dir}

TENX_IGNORE_DEPRECATED_OS=1  cellranger count --id=Homo \
--fastqs=${input_dir} \
--sample=Homo \
--transcriptome=${reference_dir} \
--include-introns false \
--localcores=10 --localmem=40 

################################################################################################################
# seurat_v5_analysis.R
################################################################################################################

################################################################################################################
# seurat_v5_analysis.R
################################################################################################################
.libPaths(new = "/lustre/home/acct-medlqian/medlqian-loop3/R/x86_64-redhat-linux-gnu-library/4.2.2")
rm(list = ls()) 
set.seed(42)

#### load packages ####
suppressPackageStartupMessages(library(Seurat)) # seurat v5
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot)) # ggplot的辅助插件
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(gridExtra)) # grid.arrange() merge multiple plot in one pdf
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grDevices)) 
suppressPackageStartupMessages(library(gridExtra))
library(DropletUtils)
# remotes::install_version("matrixStats", version="1.1.0") # restart your session and run previous scripts
library(sctransform)
library(glmGamPoi)
# devtools::install_github("immunogenomics/lisi")
library(lisi)
library(harmony)
# library(SoupX) # An R package for the estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data.
library(purrr) # purrr 是 R 语言里 tidyverse 系列的一个函数式编程工具包
library(RColorBrewer)
library(VennDiagram)
library(ggplot2)
library(reshape2)
library(ComplexUpset)
library(patchwork)
library(RColorBrewer)

setwd("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20260221_Tbr1_KO_Homo_EGFP_P_N_10X_scRNA_Seq_200G_Tbr1_KO_Het_EGFP_P_N_150G/")

#### Creat output directory ####
out_dir <- "./seurat_V5_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}
samples <- c("Het-1","Homo")
scRNA_dirs <- str_c(getwd(), "/cellranger_V7.0.1_output_renew_ref_exclude_intron/", samples,  "/outs/filtered_feature_bc_matrix/")

#### CreateSeuratObject ####
objList <- lapply(1:length(scRNA_dirs),function(i){
  CreateSeuratObject(counts = Read10X(scRNA_dirs[i]),
                     project = samples[i], assay = "RNA",
                     min.cells = 0,
                     min.features = 0) 
})
objList
# [[1]]
# An object of class Seurat 
# 32287 features across 5016 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 1 layer present: counts
# 
# [[2]]
# An object of class Seurat 
# 32287 features across 5027 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 1 layer present: counts

#### filter ####
for (i in 1:length(objList)){
  objList[[i]][["percent.mt"]] = PercentageFeatureSet(objList[[i]], pattern="^mt-")
  objList[[i]][["percent.hb"]] = PercentageFeatureSet(objList[[i]], pattern="^Hb[ab]")
}

for (i in 1:length(objList)) {
  cat("==== Sample", samples[i], "====\n")
  cat("nFeature_RNA summary:\n")
  print(summary(objList[[i]]@meta.data$nFeature_RNA))
  
  cat("nCount_RNA summary:\n")
  print(summary(objList[[i]]@meta.data$nCount_RNA))
  
  if ("percent.mt" %in% colnames(objList[[i]]@meta.data)) {
    cat("percent.mt summary:\n")
    print(summary(objList[[i]]@meta.data$percent.mt))
  }
  
  if ("percent.hb" %in% colnames(objList[[i]]@meta.data)) {
    cat("percent.hb summary:\n")
    print(summary(objList[[i]]@meta.data$percent.hb))
  }
  
  cat("\n")  # 空行分隔不同样本
}
# ==== Sample Het-1 ====
#   nFeature_RNA summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 77    2898    4360    3806    5140   10554 
# nCount_RNA summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 500   10500   17216   17056   22674  223645 
# percent.mt summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.858   2.358   3.964   3.044  86.965 
# percent.hb summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1771  0.2189  7.5435  0.2925 96.0049 
# 
# ==== Sample Homo ====
#   nFeature_RNA summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41     646    3419    2977    4589   10522 
# nCount_RNA summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 500    3666   14338   15638   21852  152443 
# percent.mt summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.6348  1.5946  3.3657  2.6925 97.4745 
# percent.hb summary:
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.03509  0.05545 12.98010  0.31512 98.28300

cell_counts_by_condition <- lapply(objList, function(obj) {
  data.frame(
    Sample = obj@project.name,
    Total = ncol(obj),
    nFeature_200_plus     = sum(obj$nFeature_RNA >= 200),
    nFeature_300_plus     = sum(obj$nFeature_RNA >= 300),
    nFeature_400_plus     = sum(obj$nFeature_RNA >= 400),
    nFeature_500_plus     = sum(obj$nFeature_RNA >= 500),
    nFeature_600_plus     = sum(obj$nFeature_RNA >= 600),
    nFeature_700_plus     = sum(obj$nFeature_RNA >= 700),
    nFeature_5000_minus   = sum(obj$nFeature_RNA <= 5000),
    nFeature_6000_minus   = sum(obj$nFeature_RNA <= 6000),
    nFeature_7000_minus   = sum(obj$nFeature_RNA <= 7000),
    nFeature_8000_minus   = sum(obj$nFeature_RNA <= 8000),
    
    nCount_700_plus       = sum(obj$nCount_RNA >= 700),
    nCount_1500_plus       = sum(obj$nCount_RNA >= 1500),
    nCount_20000_minus    = sum(obj$nCount_RNA <= 20000),
    nCount_30000_minus    = sum(obj$nCount_RNA <= 30000),
    nCount_40000_minus    = sum(obj$nCount_RNA <= 40000),
    
    percent_mt_under_10   = sum(obj$percent.mt < 10),
    
    percent_hb_under_1 = sum(obj$percent.hb < 1),
    percent_hb_under_5 = sum(obj$percent.hb < 5)
  )
})
# Warning: Default search for "data" layer in "RNA" assay yielded no results; utilizing "counts" layer instead.

cell_counts_by_condition_df <- do.call(rbind, cell_counts_by_condition)
print(cell_counts_by_condition_df)
write.csv(cell_counts_by_condition_df, file = "seurat_V5_analysis/files/differenct_filter_condition_cell_counts.csv")

for ( i in 1:length(objList)){
  pdf(str_c(out_dir,"plots/",samples[i],"_qc_plot.pdf"),width=9)
  p1=VlnPlot(objList[[i]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=ggplot(data=objList[[i]][["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
  p3=ggplot(data=objList[[i]][["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
  p4=ggplot(data=objList[[i]][["percent.mt"]],aes(x=percent.mt))+geom_density()
  p5=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}

for (i in 1:length(objList)){
  objList[[i]] <- subset(objList[[i]],subset = nFeature_RNA >= 700 & nFeature_RNA <= 8000 &
                           nCount_RNA >= 1500 & nCount_RNA <= 40000 &
                           percent.mt < 10  & percent.hb < 1)
}
objList
# [[1]]
# An object of class Seurat 
# 32287 features across 3964 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 1 layer present: counts
# 
# [[2]]
# An object of class Seurat 
# 32287 features across 3205 samples within 1 assay 
# Active assay: RNA (32287 features, 0 variable features)
# 1 layer present: counts

for ( i in 1:length(objList)){
  pdf(str_c(out_dir,"plots/",samples[i],"_filtered_qc_plot.pdf"),width=9)
  p1=VlnPlot(objList[[i]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=ggplot(data=objList[[i]][["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
  p3=ggplot(data=objList[[i]][["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
  p4=ggplot(data=objList[[i]][["percent.mt"]],aes(x=percent.mt))+geom_density()
  p5=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}

sample_summary_tbl <- map2_dfr(
  objList,
  seq_along(objList),
  ~{
    tibble(
      sample = paste0("objList[[", .y, "]]"),
      n_cells = ncol(.x),
      median_genes = median(.x$nFeature_RNA, na.rm = TRUE),
      median_umi   = median(.x$nCount_RNA,   na.rm = TRUE)
    )
  }
)

sample_summary_tbl
# # A tibble: 2 × 4
# sample       n_cells median_genes median_umi
# <chr>          <int>        <dbl>      <dbl>
#   1 objList[[1]]    3964         4645      18628
# 2 objList[[2]]    3205         4065      16900

for (i in 1:length(objList)){
  objList[[i]] <-SCTransform(objList[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
  objList[[i]] <- RunPCA(objList[[i]])
}

#### integration ####
# merge object
merge_obj <- merge(objList[[1]], y = c(objList[[2]]))
merge_obj[["RNA"]]
# Assay (v5) data with 32287 features for 7169 cells
# First 10 features:
#   Xkr4, Gm1992, Gm19938, Gm37381, Rp1, Sox17, Gm37587, Gm37323, Mrpl15, Lypla1 
# Layers:
#   counts.Het-1, counts.Homo 
merge_obj[["SCT"]]
# SCTAssay data with 20733 features for 7169 cells, and 2 SCTModel(s) 
# First 10 features:
#   Xkr4, Gm19938, Rp1, Sox17, Mrpl15, Lypla1, Tcea1, Rgs20, Gm16041, Atp6v1h 
DefaultAssay(merge_obj)

merge_obj <- SCTransform(merge_obj, vars.to.regress = c("percent.mt")) %>%
  RunPCA(verbose = FALSE)
# 本来已经在单样本时进行了SCTransform，但合并后 VariableFeatures 没被写入或被清空，直接RunPCA报错，重新运行SCTransform
# 是分两层分别进行的SCTransform

merge_obj[["RNA"]]
# Assay (v5) data with 32287 features for 7169 cells
# First 10 features:
#   Xkr4, Gm1992, Gm19938, Gm37381, Rp1, Sox17, Gm37587, Gm37323, Mrpl15, Lypla1 
# Layers:
#   counts.Het-1, counts.Homo 

merge_obj[["SCT"]]
# SCTAssay data with 20733 features for 7169 cells, and 2 SCTModel(s) 
# Top 10 variable features:
#   Mgp, Col1a1, S100a5, Col1a2, Apoe, Tyrobp, Gstm1, Cbr2, Reg3g, Krt15

# CCAIntegration
CCAIntegration_obj <- IntegrateLayers(object = merge_obj, method = CCAIntegration, 
                                      assay = "SCT", normalization.method = "SCT",
                                      orig.reduction = "pca",
                                      new.reduction = "integrated.CCA",
                                      verbose = FALSE)

#### tune the PCs ####
# tune the number of PCs
markers <- c( # OE
  "Nqo1","Ncam2", # dorsal/ventral OSNs
  "Omp", "Adcy3", "Gnal", "Cnga4", "Cnga2",  "Gng13", "Stoml3", "Rtp1", # Mature OSNs
  "Gap43", "Gnas", "Gng8", "Ablim1", "Dbn1", "Dpysl3", "Dpysl5", "Crmp1", "Marcksl1",  "Stmn1", "Stmn2", "Ebf2", "Atf5",  "Hdac2", "Cbx8",  #Immature OSNs
  "Neurog1", "Neurod1","Neurod2","Tex15", "Top2a", "Mki67", #INPs, immediate neuronal precursors 
  "Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a",#GBCs
  "Krt5","Krt14", "Trp63", #HBCs
  'Gucy2d', "Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3", # Ms4a-expressing chemosensory receptor cell (Ms4)
  "Ermn","Cyp2g1","Cyp1a2", "Sult1c1", "Gpx6", "Sox2", "Notch2", "Notch3", "Hes1", "Hes5",  "Hey1",  # Sustentacular cells
  "Coch", "Ascl3", "Cftr", "Foxi1","Hepacam2", # olfactory microvillar cells (Trpm5-)
  "Sox9", "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Ltc4s", # mv (Trpm5+)
  
  # 固有层 (lamina propria)
  'Lum', 'Dcn', 'Clec11a', # fibroblasts/stromal cells
  "Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3", # Ensheathing glia
  "Sox9","Sox10", "Krt18", "Muc5ac", "Gpx3", # Bowman's gland
  "Eng","Sox17",# Pericytes # Vascular wall cells
  "Tagln", # vascular smooth muscle cell # Vascular wall cells
  "Pecam1", "Kdr", # Endothelial cells # Vascular wall cells
  
  "Col1a1","Bglap", #Osteogenic cells
  
  # blood and immune cells
  "Cd3d", "Cd3e", "Cd8a", # cd8+ t cell
  "Cd3d", "Cd3e", "Cd4", "Il7r", # cd4+ t cell
  "Cx3cr1", # natural killer cell
  "Cd37",   "Ms4a4a", # b cell
  "Mzb1", "Sdc1", "Cd79a", # plasma cell
  "Hbb-bs","Hbb-bt",# Erythrocytes
  "S100a9","S100a8",# Neutrophils
  "Cd14",  "Lyz2", "S100a4", # monocyte
  "C1qa", "C1qb", "C1qc", "Ms4a7" # macrophage
)
markers <- markers[markers %in% rownames(CCAIntegration_obj)]
all(markers %in% rownames(CCAIntegration_obj))
unique_markers <- unique(markers)

marker_sets <- list(
  mOSN = c("Omp", "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1"),
  imOSN = c("Gap43","Gng8", "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", 
            "Dpysl3", "Stmn1", "Stmn2", "Ebf2", "Cbx8", "Trib3"),
  INP = c("Sox11","Neurog1", "Neurod1","Neurod2", "Tex15","Top2a", "Mki67", "Lhx2", "Ebf1"),
  GBC = c("Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a"),
  HBC = c("Krt5","Krt14", "Trp63", "Cxcl14", "Meg3"),
  Ms4a = c("Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3"),
  SUS = c("Sox2","Ermn","Cyp2g1","Cyp1a2", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Hey1", "Sult1c1"),
  Ensheathing = c("Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3"),
  mv1 = c("Coch", "Ascl3", "Cftr", "Hepacam2"),
  mv2 = c("Sox9", "Trpm5"),
  blood_immune = c(
    "Cd3d", "Cd3e", "Cd8a", "Cd4", "Il7r", "Cx3cr1",
    "Cd37",  "Cd79a", "Ms4a4a", "Mzb1", "Sdc1",
    "Hbb-bs","Hbb-bt", "S100a9","S100a8", "Cd14", "Clec10a",
    "Lyz2", "S100a4", "C1qa", "C1qb", "C1qc", "Ms4a7", "Hba-a1"
  ),
  EGFP_or_Cre = c("EGFP", "Cre")
)

mat <- GetAssayData(CCAIntegration_obj, assay = "SCT", layer = "data")  # 通常用 "data" slot

# 加入 meta.data
for (name in names(marker_sets)) {
  genes <- intersect(marker_sets[[name]], rownames(mat)) 
  CCAIntegration_obj[[paste0(name, "_marker_sum")]] <- Matrix::colSums(mat[genes, , drop = FALSE])
}

for (nPCs in seq(20,50,5)){
  CCAIntegration_obj <- FindNeighbors(CCAIntegration_obj, dims = 1:nPCs, reduction = "integrated.CCA")
  CCAIntegration_obj <- FindClusters(CCAIntegration_obj, resolution = 0.2, reduction = "integrated.CCA")
  CCAIntegration_obj <- RunUMAP(CCAIntegration_obj, dims = 1:nPCs, reduction = "integrated.CCA", reduction.name = "umap")
  
  DefaultAssay(CCAIntegration_obj) <- "SCT"
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(CCAIntegration_obj, reduction = "umap", label=TRUE, raster=FALSE) + theme(legend.position = "bottom")
  p2 <- DimPlot(CCAIntegration_obj, reduction = "umap", group.by="orig.ident", raster=FALSE) + theme(legend.position = "bottom")
  print(p1)
  print(p2)
  dev.off()
  
  DefaultAssay(CCAIntegration_obj) <- "SCT"
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_markers_sum_UMAP.pdf"))
  p <- FeaturePlot(CCAIntegration_obj, reduction = "umap", features=paste0(names(marker_sets),"_marker_sum"), combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
  print(p)
  dev.off()
  
  DefaultAssay(CCAIntegration_obj) <- "SCT"
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(CCAIntegration_obj, reduction = "umap", features=unique_markers, combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
  print(p)
  dev.off()
  
  DefaultAssay(CCAIntegration_obj) <- "SCT"
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_CD36_UMAP.pdf"))
  p <- FeaturePlot(CCAIntegration_obj, reduction = "umap", features='Cd36', combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
  print(p)
  dev.off()
}

#### tune the resolution (35 PCs) ####
CCAIntegration_obj <- FindNeighbors(CCAIntegration_obj, reduction = "integrated.CCA", dims = 1:35)
CCAIntegration_obj <- RunUMAP(CCAIntegration_obj, dims = 1:35, reduction = "integrated.CCA", reduction.name = "umap.CCA")

pdf("seurat_V5_analysis/plots/integrated_OE_scRNA_PC35_markers_UMAP_without_order.pdf")
p <- FeaturePlot(CCAIntegration_obj, reduction = "umap.CCA", features=unique_markers, combine = F,
                 cols=c("lightgrey", "red"), order=F, raster=FALSE)
print(p)
dev.off()

pdf("seurat_V5_analysis/plots/integrated_OE_scRNA_PC35_markers_splitByGenotype_UMAP_without_order.pdf", width = 12, height = 6)
for (markers in unique_markers) {
  p <-FeaturePlot(CCAIntegration_obj, reduction = "umap.CCA", features=markers, 
                  cols=c("lightgrey", "red"), order=F, split.by = "orig.ident")
  print(p)
}
dev.off()

pdf("seurat_V5_analysis/plots/integrated_OE_scRNA_PC35_marker_sum_splitByGenotype_UMAP.pdf", width = 12, height = 6)
for (name in names(marker_sets)) {
  plot_feature <- paste0(name, "_marker_sum")  # 构建 meta.data 中的列名
  p <- FeaturePlot(CCAIntegration_obj, reduction = "umap.CCA", features = plot_feature, order = F,
                   split.by = "orig.ident")
  print(p)
}
dev.off()

# 运行多个 resolution
resolutions <- seq(0.1, 1.5, by = 0.1)
for (res in resolutions) {
  CCAIntegration_obj <- FindClusters(CCAIntegration_obj, resolution = res)
  
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC35_resolution",res,"_UMAP.pdf"))
  p1 <- DimPlot(CCAIntegration_obj, reduction = "umap.CCA", label=TRUE, group.by=str_c("SCT_snn_res.", res))
  p2 <- DimPlot(CCAIntegration_obj, reduction = "umap.CCA", group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  
  pdf(width = 15, height = 80, str_c(out_dir,"plots/integrated_OE_scRNA_PC35_resolution", res,"_vlnplot_marker.pdf"))
  p1 <- VlnPlot(CCAIntegration_obj, group.by=str_c("SCT_snn_res.", res), features = unique_markers, pt.size = 0, ncol = 3) &
    xlab(NULL) & ylab(NULL)
  print(p1)
  dev.off()
  
  pdf(width = 24, height = 20, str_c(out_dir,"plots/integrated_OE_scRNA_PC35_resolution",res,"_dotplot_marker.pdf"))
  p2 <- DotPlot(CCAIntegration_obj, group.by=str_c("SCT_snn_res.", res), features = rev(unique_markers), cols = c("lightgrey", "dodgerblue3")) +
    RotatedAxis()+
    coord_flip ()
  print(p2)
  dev.off()
}

save(CCAIntegration_obj, file = paste0(out_dir, "rdata/CCAIntegration_obj_35PCs_different_rsolutions.rdata"))

#### Findmarker (35 PCs resolution 0.9) #### 
load("seurat_V5_analysis/rdata/CCAIntegration_obj_35PCs_different_rsolutions.rdata")
Idents(CCAIntegration_obj) <- "SCT_snn_res.0.9"

CCAIntegration_obj <- PrepSCTFindMarkers(CCAIntegration_obj)
# 在你用 SCTransform（尤其是多样本/多 layer → 多个 SCTModel）之后，帮你把 FindMarkers/FindAllMarkers 需要的 SCT 统计量“对齐并校准”，避免因为不同 SCT 模型/不同测序深度校正基线不一致而产生偏差或直接报错
# SCT 的表达不是简单的 log-normalized counts，而是基于 NB 回归得到的残差/校正值，并且每个样本（或每个 layer）可能有自己的一套 SCTModel。
# PrepSCTFindMarkers 统一/校准用于 DE 的计数基线（library size / corrected counts）它会根据每个 SCTModel 的回归结果，生成/更新用于差异分析的一些内部量（常见叫法是 corrected counts 或与 median_umi 相关的校正），让不同模型之间的 DE 有一个一致的参照。

clusters <- levels(Idents(CCAIntegration_obj))

markers_list <- map(clusters, ~ FindMarkers(
  CCAIntegration_obj,
  ident.1 = .x,
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
))
# ROC 非常合适找高 AUC 的分离基因，而不是严格的差异表达显著性。
# roc 在 Seurat 里主要基于表达值做 AUC 排序（更像分类打分），它不走需要“归一化方法选择”的那条逻辑分支，所以 norm.method 不参与计算。

names(markers_list) <- clusters
save(markers_list, file = "seurat_V5_analysis/rdata/CCAIntegration_obj_35PCs_res.0.9_clusters_markers_list.rdata")

load("seurat_V5_analysis/rdata/CCAIntegration_obj_35PCs_res.0.9_clusters_markers_list.rdata")
markers_df <- bind_rows(
  lapply(clusters, function(cl) {
    df <- markers_list[[cl]]
    df$cluster <- cl
    df$gene <- rownames(df)
    df
  }),
  .id = "cluster_id"
)
write.csv(markers_df, file = "seurat_V5_analysis/files/CCAIntegration_obj_35PCs_res.0.9_clusters_markers_df.csv")

head(markers_df[which(markers_df$cluster == 12),], n = 100)

#### annotate the cluster (35 PCs 0.9 resolution) ####
# mOSN: 0, 2, 
# iOSN: 1, 3, 
# late INP: 9, 
# early INP: 15, 
# GBC: 10, 16, 17
# HBC: 5, 23, 
# SUS: 8, 11
# MV: 18

# OEC: 19
# vasulature: 4, 6, 22
# immune cell: 13, 14

# respiratory epithelial cell： 12, 21
# Immature excitatory neuron: 7
# Immature interneuron: 20

Idents(CCAIntegration_obj) <- "SCT_snn_res.0.9"
CCAIntegration_obj$celltype <- case_when(
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(0, 2) ~ "mOSN",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(1, 3) ~ "iOSN",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 9 ~ "late INP",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 15 ~ "early INP",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(10, 16, 17) ~ "GBC",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(5, 23) ~ "HBC",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(8, 11) ~ "SUS",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 18 ~ "MV",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 19 ~ "OEC",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(4, 6, 22) ~ "vasculature cell",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(13, 14) ~ "immune cell",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% c(12, 21) ~ "respiratory epithelial cell",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 7 ~ "immature excitatory neuron",
  CCAIntegration_obj$SCT_snn_res.0.9 %in% 20 ~ "immature interneuron",
  TRUE ~ NA_character_
)

CCAIntegration_obj$celltype <- factor(CCAIntegration_obj$celltype, levels = c("mOSN", "iOSN", "late INP","early INP", "GBC", "HBC", 
                                                                              "SUS","MV", 
                                                                              "OEC", "vasculature cell", "immune cell", 
                                                                              "respiratory epithelial cell", 
                                                                              "immature excitatory neuron", "immature interneuron"))

CCAIntegration_obj$celltype_number <- paste0(CCAIntegration_obj$celltype, "_", CCAIntegration_obj$SCT_snn_res.0.9)

CCAIntegration_obj$celltype_number <- factor(CCAIntegration_obj$celltype_number, 
                                             levels = c(paste0("mOSN", "_", c(0, 2)),
                                                        paste0("iOSN", "_", c(1, 3)), 
                                                        paste0("late INP", "_", 9), 
                                                        paste0("early INP", "_", 15), 
                                                        paste0("GBC","_",  c(10, 16, 17)),
                                                        paste0("HBC","_",  c(5, 23)),
                                                        paste0("SUS","_",  c(8, 11)),
                                                        paste0("MV", "_", 18),
                                                        paste0("OEC","_",  19),
                                                        paste0("vasculature cell","_", c(4, 6, 22)),
                                                        paste0("immune cell", "_", c(13, 14)),
                                                        paste0("respiratory epithelial cell", "_", c(12, 21)),
                                                        paste0("immature excitatory neuron","_", 7),
                                                        paste0("immature interneuron","_", 20)))
save(CCAIntegration_obj, file = "seurat_V5_analysis/rdata/CCAIntegration_obj_annotated_cluster.rdata")

#### umap ####
pdf("seurat_V5_analysis/plots/umap_CCAIntegration_obj_cluster_or_celltype.pdf", height = 6, width = 8)
DimPlot(CCAIntegration_obj, reduction = "umap.CCA", group.by = "celltype", label = F)
DimPlot(CCAIntegration_obj, reduction = "umap.CCA", group.by = "SCT_snn_res.0.9", label = T)
dev.off()

pdf("seurat_V5_analysis/plots/umap_CCAIntegration_obj_cluster_or_celltype_recolor.pdf", height = 6, width = 8)
DimPlot(CCAIntegration_obj, reduction = "umap.CCA", group.by = "celltype", 
        label = F, cols = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#ff7f00',
                            '#cab2d6','#6a3d9a',
                            '#cccc00','#b15928',
                            '#FF66B2','#66c2a5','#8da0cb','#b266ff'))

DimPlot(CCAIntegration_obj, reduction = "umap.CCA", 
        group.by = "SCT_snn_res.0.9", label = T)
dev.off()

#### markers Dotplot ####
markers_to_display <- c("Omp", "Gap43", "Tex15", "Neurod1", "Top2a", "Ascl1","Krt5",
                        "Cyp2g1", "Coch", "Atp1a2", "Dcn", "Eng", "Cd14", "Krt19", "Eomes", "Dlx5")

Idents(CCAIntegration_obj) <- "celltype_number"

CCAIntegration_obj$celltype_number <- factor(CCAIntegration_obj$celltype_number, 
                                             levels = rev(levels(CCAIntegration_obj$celltype_number)))


pdf("seurat_V5_analysis/plots/dotplot_selected_markers.pdf", width = 12, height = 8)
DotPlot(CCAIntegration_obj, features = markers_to_display) + 
  RotatedAxis() + 
  scale_color_gradient(low = "grey90", high = '#3a78d7')
dev.off()

CCAIntegration_obj$celltype_number <- factor(CCAIntegration_obj$celltype_number, 
                                             levels = rev(levels(CCAIntegration_obj$celltype_number)))

#### cell count barplot ####
celltype_cell_number_df <- table(
  celltype = CCAIntegration_obj$celltype,
  genotype   = CCAIntegration_obj$orig.ident
) %>%
  as.data.frame() %>%
  dplyr::rename(cellnumber = Freq)
# celltype genotype cellnumber
# 1                         mOSN  Het-1       1024
# 2                         iOSN  Het-1        833
# 3                     late INP  Het-1        212
# 4                    early INP  Het-1        130
# 5                          GBC  Het-1        333
# 6                          HBC  Het-1        202
# 7                          SUS  Het-1        184
# 8                           MV  Het-1         54
# 9                          OEC  Het-1         67
# 10            vasculature cell  Het-1        339
# 11                 immune cell  Het-1        117
# 12 respiratory epithelial cell  Het-1        144
# 13  immature excitatory neuron  Het-1        268
# 14        immature interneuron  Het-1         57
# 15                        mOSN   Homo        520
# 16                        iOSN   Homo        485
# 17                    late INP   Homo        106
# 18                   early INP   Homo         15
# 19                         GBC   Homo        234
# 20                         HBC   Homo        249
# 21                         SUS   Homo        428
# 22                          MV   Homo         50
# 23                         OEC   Homo         23
# 24            vasculature cell   Homo        490
# 25                 immune cell   Homo        292
# 26 respiratory epithelial cell   Homo        181
# 27  immature excitatory neuron   Homo        100
# 28        immature interneuron   Homo         32

pdf('seurat_V5_analysis/plots/barplot_celltype_percent_per_genotype.pdf', width = 8, height = 6)
ggplot(celltype_cell_number_df, aes(x = genotype, y = cellnumber, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(celltype_cell_number_df$celltype)))) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Group", y = "Cell type proportion", fill = "Cell type",
       title = "Cell type composition in Het vs Homo") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(celltype_cell_number_df, aes(x = genotype, y = cellnumber, fill = celltype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(celltype_cell_number_df$celltype)))) +
  labs(
    x = "Cell type",
    y = "Cell number",
    fill = "Group",
    title = "Cell type distribution by group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()

#### QC plot ####
load("seurat_V5_analysis/rdata/CCAIntegration_obj_annotated_cluster.rdata")

meta_Het <- CCAIntegration_obj@meta.data %>%
  filter(orig.ident == "Het-1") 

meta_Homo <- CCAIntegration_obj@meta.data %>%
  filter(orig.ident == "Homo") 

meta_combined <- rbind(meta_Het, meta_Homo)

# 比较 nCount_RNA 和 nFeature_RNA
summary_table <- data.frame(
  Metric = c("nCount_RNA", "nFeature_RNA"),
  Het_Median = c(median(meta_Het$nCount_RNA), median(meta_Het$nFeature_RNA)),
  Homo_Median = c(median(meta_Homo$nCount_RNA), median(meta_Homo$nFeature_RNA)),
  p_value_Wilcox = c(
    wilcox.test(meta_Het$nCount_RNA, meta_Homo$nCount_RNA)$p.value,
    wilcox.test(meta_Het$nFeature_RNA, meta_Homo$nFeature_RNA)$p.value
  )
)
summary_table
# Metric Het_Median Homo_Median p_value_Wilcox
# 1   nCount_RNA      18628       16900   4.892133e-13
# 2 nFeature_RNA       4645        4065   2.558032e-86

median_df <- meta_combined %>%
  group_by(orig.ident) %>%
  summarise(nCount_RNA_med = median(nCount_RNA),
            nFeature_RNA_med = median(nFeature_RNA))

pdf("seurat_V5_analysis/plots/QC_plots.pdf", width = 10, height = 6)
p1 <- ggplot(meta_combined, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  geom_text(data = median_df, 
            aes(x = orig.ident, y = nCount_RNA_med, label = paste0("Median: ", round(nCount_RNA_med, 0))),
            vjust = -0.8, size = 4, color = "black") +
  scale_y_log10() +
  # 将 Y 轴按 log10 缩放，非常常用于展示 RNA count，因为其分布极不均匀。
  labs(y = "nCount_RNA (log10)", x = "", title = "Sequencing Depth per Cell") +
  theme_classic(base_size = 10) +
  # base_size = 10：默认字体大小为 10。
  theme(legend.position = "none")

p2 <- ggplot(meta_combined, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_y_log10() +
  geom_text(data = median_df, aes(x = orig.ident, y = nFeature_RNA_med, label = paste0("Median: ", round(nFeature_RNA_med, 0))),
            vjust = -0.8, size = 4, color = "black") +
  labs(y = "nFeature_RNA (log10)", x = "", title = "Detected Genes per Cell") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

p1 + p2

p3 <- ggplot(meta_combined, aes(x = nCount_RNA, color = orig.ident)) +
  stat_ecdf(linewidth = 1.2) +
  # ECDF（empirical cumulative distribution function, 经验累积分布曲线）
  # linewidth = 1.2：设置曲线线条宽度为 1.2，使曲线更醒目。
  scale_x_log10() +
  labs(x = "nCount_RNA (log10)", y = "Cumulative Fraction", title = "Cumulative Distribution of nCount_RNA") +
  theme_classic(base_size = 10)

p4 <- ggplot(meta_combined, aes(x = nFeature_RNA, color = orig.ident)) +
  stat_ecdf(linewidth = 1.2) +
  scale_x_log10() +
  labs(x = "nFeature_RNA (log10)", y = "Cumulative Fraction", title = "Cumulative Distribution of nFeature_RNA") +
  theme_classic(base_size = 10)

p3 + p4
dev.off()

#### tune the UMI threshold of receptor expression in Het_mOSN ####
Idents(CCAIntegration_obj) <- "celltype"
mOSN <- subset(CCAIntegration_obj, idents = "mOSN")
Idents(mOSN) <- "orig.ident"
Het_mOSN <- subset(mOSN, idents = "Het-1")
Het_mOSN_scRNA <- DietSeurat(
  object=Het_mOSN,
  layers = c("counts", "data"),
  assays="SCT"
)
# DietSeurat
# Save memory
# Reduce file size for sharing/saving
# Speed up processing
# Focus only on specific components of the data

# identiﬁcation of the OR expressed in each OSN
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

setdiff(grep("Olfr", all_functional_olfactory_receptors$gene_name, value = TRUE),
        c(ClassI_OR_name, ClassII_OR_name))
# character(0)

Het_mOSN_receptor_counts <- FetchData(
  object = Het_mOSN_scRNA,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(Het_mOSN_scRNA)),
  layer = "counts"
)

dim(Het_mOSN_receptor_counts)
# 1024  642

Het_mOSN_calculate_OSNs_percentage <- function(cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(Het_mOSN_receptor_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(Het_mOSN_receptor_counts>=cutoff)
  }
  
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  
  Het_mOSN_receptor_threshold_df <- data.frame(cutoff=cutoff,
                                               Olfr_N=c("0","1",">1"),
                                               percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(Het_mOSN_receptor_threshold_df)
}

cutoffs <- c(0:8,16,32,64,128,256)
Het_mOSN_receptor_threshold_df <- c()
for (i in 1:length(cutoffs)){
  Het_mOSN_receptor_threshold_df <- Het_mOSN_calculate_OSNs_percentage(cutoffs[i]) %>% rbind(Het_mOSN_receptor_threshold_df, .)
}
Het_mOSN_receptor_threshold_df$cutoff <- factor(Het_mOSN_receptor_threshold_df$cutoff,levels=cutoffs)
Het_mOSN_receptor_threshold_df$Olfr_N <- factor(Het_mOSN_receptor_threshold_df$Olfr_N,levels=c("0","1",">1"))
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])

pdf("seurat_V5_analysis/plots/Het_mOSNs_UMIs_cutoff_tunning.pdf", height = 8, width=9)
ggplot(data=Het_mOSN_receptor_threshold_df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of mOSNs",color="# of ORs")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=rel(2)),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_text(color="black",size=rel(1.8)),
        axis.line = element_line(colour="black",linewidth = 1))
dev.off()


#### tune the UMI threshold of receptor expression in Het_iOSN ####
Idents(CCAIntegration_obj) <- "celltype"
iOSN <- subset(CCAIntegration_obj, idents = "iOSN")
Idents(iOSN) <- "orig.ident"
Het_iOSN <- subset(iOSN, idents = "Het-1")
Het_iOSN_scRNA <- DietSeurat(
  object=Het_iOSN,
  layers = c("counts", "data"),
  assays="SCT"
)
# DietSeurat
# Save memory
# Reduce file size for sharing/saving
# Speed up processing
# Focus only on specific components of the data

Het_iOSN_receptor_counts <- FetchData(
  object = Het_iOSN_scRNA,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(Het_iOSN_scRNA)),
  layer = "counts"
)
dim(Het_iOSN_receptor_counts)
# 835 642
length(intersect(colnames(Het_iOSN_receptor_counts), ClassI_OR_name))
# 86
length(intersect(colnames(Het_iOSN_receptor_counts), ClassII_OR_name))
# 543
length(intersect(colnames(Het_iOSN_receptor_counts), Taar_name))
# 13

Het_iOSN_calculate_OSNs_percentage <- function(cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(Het_iOSN_receptor_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(Het_iOSN_receptor_counts>=cutoff)
  }
  
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  
  Het_iOSN_receptor_threshold_df <- data.frame(cutoff=cutoff,
                                               Olfr_N=c("0","1",">1"),
                                               percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(Het_iOSN_receptor_threshold_df)
}

cutoffs <- c(0:8,16,32,64,128,256)
Het_iOSN_receptor_threshold_df <- c()
for (i in 1:length(cutoffs)){
  Het_iOSN_receptor_threshold_df <- calculate_OSNs_percentage(cutoffs[i]) %>% rbind(Het_iOSN_receptor_threshold_df, .)
}
Het_iOSN_receptor_threshold_df$cutoff <- factor(Het_iOSN_receptor_threshold_df$cutoff,levels=cutoffs)
Het_iOSN_receptor_threshold_df$Olfr_N <- factor(Het_iOSN_receptor_threshold_df$Olfr_N,levels=c("0","1",">1"))
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])

pdf("seurat_V5_analysis/plots/Het_iOSNs_UMIs_cutoff_tunning.pdf", height = 8, width=9)
ggplot(data=Het_iOSN_receptor_threshold_df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of iOSNs",color="# of ORs")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=rel(2)),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_text(color="black",size=rel(1.8)),
        axis.line = element_line(colour="black",linewidth = 1))
dev.off()

#### Cell counts expressing different numbers of receptors (mOSN) (UMI > 1) ####
# 获取每个细胞所属组
Idents(CCAIntegration_obj) <- "celltype"
mOSN <- subset(CCAIntegration_obj, idents = "mOSN")
mOSN_receptor_counts <-  FetchData(
  object = mOSN,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(mOSN)),
  layer = "counts"
)
dim(mOSN_receptor_counts)
# 1544  642

cell_group <- mOSN$orig.ident

# 计算每个细胞表达的受体基因数量
mOSN_receptor_count_per_cell <- rowSums(mOSN_receptor_counts > 1)

#  合并成 data.frame
mOSN_receptor_expr_summary <- data.frame(
  cell = names(mOSN_receptor_count_per_cell),
  receptor_count = mOSN_receptor_count_per_cell,
  group = cell_group
)

# 分组统计
mOSN_receptor_summary_by_group <- table(mOSN_receptor_expr_summary$group, mOSN_receptor_expr_summary$receptor_count)
mOSN_receptor_summary_by_group
# 0   1   2   3   4
# Het-1 108 841  69   5   1
# Homo  206 276  37   1   0

# 转为长格式
mOSN_receptor_summary_by_group_df <- as.data.frame(mOSN_receptor_summary_by_group)
colnames(mOSN_receptor_summary_by_group_df) <- c("Group", "receptor_count", "Cell_number")

mOSN_receptor_summary_prop_df <- mOSN_receptor_summary_by_group_df %>%
  group_by(Group) %>%
  mutate(Percent = Cell_number / sum(Cell_number))
# Group receptor_count Cell_number  Percent
# <fct> <fct>                <int>    <dbl>
#   1 Het-1 0                      108 0.105   
# 2 Homo  0                      206 0.396   
# 3 Het-1 1                      841 0.821   
# 4 Homo  1                      276 0.531   
# 5 Het-1 2                       69 0.0674  
# 6 Homo  2                       37 0.0712  
# 7 Het-1 3                        5 0.00488 
# 8 Homo  3                        1 0.00192 
# 9 Het-1 4                        1 0.000977
# 10 Homo  4                        0 0 
# 绘制堆叠条形图
pdf("seurat_V5_analysis/plots/barplot_receptor_number_class_mOSN_UMI_2.pdf", width = 6, height = 6)
ggplot(mOSN_receptor_summary_prop_df, aes(x = Group, y = Cell_number, fill = as.factor(receptor_count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(mOSN_receptor_summary_prop_df$receptor_count)))) +
  labs(
    title = "The number of cells expressing \ndifferent receptor numbers",
    x = "Group",
    y = "Number of cells",
    fill = "Number of expressed receptor genes (UMI > 1)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggplot(mOSN_receptor_summary_prop_df, aes(x = Group, y = Percent, fill = as.factor(receptor_count))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(mOSN_receptor_summary_prop_df$receptor_count)))) +
  labs(
    title = "Proportion of cells expressing \ndifferent receptor numbers",
    x = "Group",
    y = "Percentage of cells (%)",
    fill = "Number of expressed receptor genes (UMI > 1)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()

#### Cell counts expressing different classes of receptors (mOSN expressing single receptor) (UMI > 1) ####
# 定义函数：提取表达单个受体mOSN并分类
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

# 调用函数：Het 与 Homo
single_Het <- get_single_receptor_mOSN(mOSN, "Het-1")
single_Homo <- get_single_receptor_mOSN(mOSN, "Homo")

# 合并结果
receptor_class_df <- bind_rows(
  data.frame(group = "Het", receptor_class = single_Het$receptor_class),
  data.frame(group = "Homo", receptor_class = single_Homo$receptor_class)
)

# 绘制比例图
receptor_class_summary <- receptor_class_df %>%
  group_by(group, receptor_class) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  # n() 统计每个小组中的行数
  # .groups = "drop_last": summarise 后保留第一层 group（这里是 group），但移除最后一层（receptor_class）。
  mutate(prop = count / sum(count))
# A tibble: 5 × 4
# Groups:   group [2]
# group receptor_class count  prop
# <chr> <chr>          <int> <dbl>
#   1 Het   Class_II_OR      588 0.699
# 2 Het   Class_I_OR       110 0.131
# 3 Het   Taar             143 0.170
# 4 Homo  Class_II_OR      162 0.587
# 5 Homo  Class_I_OR       114 0.413

pdf("seurat_V5_analysis/plots/barplot_receptor_class_in_mOSN_expressing_single_receptor_UMI_2.pdf",width = 6, height = 6)
ggplot(receptor_class_summary, aes(x = group, y = count, fill = receptor_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Taar" = "#E64B35", "Class_I_OR" = "#4DBBD5",
    "Class_II_OR" = "#00A087", "None" = "#999999"
  )) +
  labs(
    y = "Number of single-receptor mOSNs",
    fill = "Receptor class",
    title = "Receptor class count in Het vs Homo mOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(receptor_class_summary, aes(x = group, y = count, fill = receptor_class)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = c(
    "Taar" = "#E64B35", "Class_I_OR" = "#4DBBD5",
    "Class_II_OR" = "#00A087", "None" = "#999999"
  )) +
  labs(
    y = "Proportion of receptor classes",
    fill = "Receptor class",
    title = "Receptor class distribution in Het vs Homo mOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#### Cell counts expressing different numbers of receptors (EGFP_Cre_mOSN) (UMI > 1) ####
# 获取每个细胞所属组
EGFP_Cre_mOSN <- subset(mOSN, subset = EGFP_or_Cre_marker_sum > 0)
save(EGFP_Cre_mOSN, file = "seurat_V5_analysis/rdata/EGFP_Cre_mOSN.rdata")

EGFP_Cre_mOSN_receptor_counts <- FetchData(
  EGFP_Cre_mOSN,
  vars =  intersect(all_functional_olfactory_receptors$gene_name, rownames(EGFP_Cre_mOSN)),
  assay = "SCT",
  layer = "counts"
)

cell_group <- EGFP_Cre_mOSN$orig.ident

# 计算每个细胞表达的基因数量
receptor_count_per_cell <- rowSums(EGFP_Cre_mOSN_receptor_counts > 1)

#  合并成 data.frame
receptor_expr_summary <- data.frame(
  cell = names(receptor_count_per_cell),
  receptor_count = receptor_count_per_cell,
  group = cell_group
)

# 分组统计
receptor_summary_by_group <- table(receptor_expr_summary$group, receptor_expr_summary$receptor_count)
receptor_summary_by_group
# 0   1   2
# Het-1   8 101   4
# Homo   27  64  10

# 转为长格式
receptor_summary_by_group_df <- as.data.frame(receptor_summary_by_group)
colnames(receptor_summary_by_group_df) <- c("Group", "receptor_count", "Cell_number")

receptor_summary_prop_df <- receptor_summary_by_group_df %>%
  group_by(Group) %>%
  mutate(Percent = Cell_number / sum(Cell_number))
# # A tibble: 6 × 4
# # Groups:   Group [2]
# Group receptor_count Cell_number Percent
# <fct> <fct>                <int>   <dbl>
#   1 Het-1 0                        8  0.0708
# 2 Homo  0                       27  0.267 
# 3 Het-1 1                      101  0.894 
# 4 Homo  1                       64  0.634 
# 5 Het-1 2                        4  0.0354
# 6 Homo  2                       10  0.0990

pdf("seurat_V5_analysis/plots/barplot_receptor_number_class_EGFP_Cre_mOSN_UMI_2.pdf",width = 6, height = 6)

ggplot(receptor_summary_by_group_df, aes(x = Group, y = Cell_number, fill = as.factor(receptor_count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "Paired"))(length(unique(receptor_summary_by_group_df$receptor_count)))) +
  labs(
    title = "The number of cells \ndifferent receptor numbers",
    x = "Group",
    y = "Number of cells",
    fill = "receptor number (UMI > 1)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggplot(receptor_summary_prop_df, aes(x = Group, y = Percent, fill = as.factor(receptor_count))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "Paired"))(length(unique(receptor_summary_by_group_df$receptor_count)))) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(
    title = "Proportion of cells expressing \ndifferent receptor numbers",
    x = "Group",
    y = "Percentage of cells (%)",
    fill = "Number of expressed receptor genes (UMI > 1)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()

#### Cell counts expressing different classes of receptors (EGFP_Cre_mOSN expressing single receptor) (UMI > 1) ####
# 调用函数：Het 与 Homo
single_Het <- get_single_receptor_mOSN(EGFP_Cre_mOSN, "Het-1")
single_Homo <- get_single_receptor_mOSN(EGFP_Cre_mOSN, "Homo")

# 合并结果
receptor_class_df <- bind_rows(
  data.frame(group = "Het", receptor_class = single_Het$receptor_class),
  data.frame(group = "Homo", receptor_class = single_Homo$receptor_class)
)

# 绘制比例图
receptor_class_summary <- receptor_class_df %>%
  group_by(group, receptor_class) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(prop = count / sum(count))
# # A tibble: 5 × 4
# # Groups:   group [2]
# group receptor_class count  prop
# <chr> <chr>          <int> <dbl>
#   1 Het   Class_II_OR       11 0.109
# 2 Het   Class_I_OR        11 0.109
# 3 Het   Taar              79 0.782
# 4 Homo  Class_II_OR       18 0.281
# 5 Homo  Class_I_OR        46 0.719
write.csv(receptor_class_summary, file = "seurat_V5_analysis/files/receptor_class_percent_EGFP_Cre_mOSN_expressing_single_receptor.csv")

pdf("seurat_V5_analysis/plots/barplot_receptor_class_in_single_receptor_EGFP_Cre_mOSN_UMI_2.pdf",width = 6, height = 6)
ggplot(receptor_class_summary, aes(x = group, y = count, fill = receptor_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Taar" = "#E64B35", "Class_I_OR" = "#4DBBD5",
    "Class_II_OR" = "#00A087", "None" = "#999999"
  )) +
  labs(
    y = "Number of single-receptor EGFP_Cre_mOSNs",
    fill = "Receptor class",
    title = "Receptor class count in \nHet vs Homo EGFP_Cre_mOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(receptor_class_summary, aes(x = group, y = count, fill = receptor_class)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Taar" = "#E64B35", "Class_I_OR" = "#4DBBD5",
    "Class_II_OR" = "#00A087", "None" = "#999999"
  )) +
  labs(
    y = "Proportion of receptor classes",
    fill = "Receptor class",
    title = "Receptor class distribution in \nHet vs Homo EGFP_Cre_mOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


#### Cell counts expressing different numbers of receptors (iOSN) (UMI > 0) ####
Idents(CCAIntegration_obj) <- "celltype"
iOSN <- subset(CCAIntegration_obj, idents = "iOSN")

iOSN_receptor_counts <- FetchData(
  iOSN,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(iOSN)),
  assay = "SCT",
  layer = "counts"
)

dim(iOSN_receptor_counts)
# 1318  642

iOSN$receptor_number <- rowSums(iOSN_receptor_counts > 0)
table(iOSN$receptor_number)
# 0   1   2   3   4   5   6   7   8   9  12  18 
# 250 291 262 211 149  93  34  18   3   5   1   1 

table(iOSN$orig.ident)
# Het-1  Homo 
# 833   485  

df <- table(iOSN$orig.ident, iOSN$receptor_number) %>% 
  as.data.frame()
colnames(df) <- c("group", "receptor_number", "cell_count")

df_prop <- df %>%
  group_by(group) %>%
  mutate(percent = cell_count / sum(cell_count) * 100)

df_prop <- df_prop %>%
  mutate(percent = round(percent, 2))

df_prop$receptor_number <- factor(df_prop$receptor_number, levels = sort(unique(df_prop$receptor_number)))

pdf("seurat_V5_analysis/plots/barplot_receptor_number_class_iOSN.pdf", width = 6, height = 6)
ggplot(df, aes(x = factor(receptor_number), y = cell_count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Het-1" = "#4DBBD5", "Homo" = "#E64B35")) +
  labs(
    x = "Number of receptors per cell",
    y = "Cell count",
    fill = "Group",
    title = "Distribution of receptor numbers in iOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(df_prop, aes(x = group, y = percent, fill = factor(receptor_number))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df_prop$receptor_number)))) +
  labs(
    x = "",
    y = "Percentage of cells (%)",
    fill = "Number of receptors per cell",
    title = "receptor numbers in Het vs Homo iOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

dev.off()


#### Cell counts expressing different numbers of receptors (EGFP_Cre_iOSN) (UMI > 0) ####
EGFP_Cre_iOSN <- subset(iOSN, subset = EGFP_or_Cre_marker_sum > 0)

EGFP_Cre_iOSN_receptor_counts <- FetchData(
  EGFP_Cre_iOSN,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(EGFP_Cre_iOSN)),
  assay = "SCT",
  layer = "counts"
)

dim(EGFP_Cre_iOSN_receptor_counts)
# 571 642

EGFP_Cre_iOSN$receptor_number <- rowSums(EGFP_Cre_iOSN_receptor_counts > 0)
table(EGFP_Cre_iOSN$receptor_number)
# 0   1   2   3   4   5   6   7   9  18 
# 155 140 112  75  41  32   7   6   2   1 

table(EGFP_Cre_iOSN$orig.ident)
# Het-1  Homo 
# 269   302

df <- table(EGFP_Cre_iOSN$orig.ident, EGFP_Cre_iOSN$receptor_number) %>% 
  as.data.frame()
colnames(df) <- c("group", "receptor_number", "cell_count")

df_prop <- df %>%
  group_by(group) %>%
  mutate(percent = cell_count / sum(cell_count) * 100)

df_prop <- df_prop %>%
  mutate(percent = round(percent, 2))

df_prop$receptor_number <- factor(df_prop$receptor_number, levels = sort(unique(df_prop$receptor_number)))

pdf("seurat_V5_analysis/plots/barplot_cell_count_per_receptor_number_EGFP_Cre_iOSN.pdf")
# 绘图
ggplot(df, aes(x = factor(receptor_number), y = cell_count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Het-1" = "#4DBBD5", "Homo" = "#E64B35")) +
  labs(
    x = "Number of receptors per cell",
    y = "Cell count",
    fill = "Group",
    title = "Distribution of receptor numbers in iOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(df_prop, aes(x = group, y = percent, fill = factor(receptor_number))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df_prop$receptor_number)))) +
  labs(
    x = "",
    y = "Percentage of cells (%)",
    fill = "Number of receptors per cell",
    title = "receptor numbers in Het vs Homo iOSN"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()

#### Taar/OR coexpression in Het iOSNs expressing multiple receptors (Het iOSN) (UMI > 0) ####
Idents(CCAIntegration_obj) <- "celltype"
iOSN <- subset(CCAIntegration_obj, idents = "iOSN")
iOSN_receptor_counts <- FetchData(
  iOSN,
  vars = intersect(all_functional_olfactory_receptors$gene_name, rownames(iOSN)),
  assay = "SCT",
  layer = "counts"
)
iOSN$receptor_number <- rowSums(iOSN_receptor_counts > 0)
table(iOSN$orig.ident)
# Het-1  Homo 
# 833   485
iOSN_multiple_receptors <- subset(iOSN, subset = receptor_number > 1)
table(iOSN_multiple_receptors$orig.ident)
# Het-1  Homo 
# 667   110  

# Taar
iOSN_Taar_counts <- FetchData(
  iOSN_multiple_receptors,
  vars = grep("Taar", rownames(iOSN_multiple_receptors), value = T),
  assay = "SCT",
  layer = "counts"
)

iOSN_multiple_receptors$Taar_number <- rowSums(iOSN_Taar_counts > 0)

# Olfr
iOSN_Olfr_counts <- FetchData(
  iOSN_multiple_receptors,
  vars = grep("Olfr", rownames(iOSN_multiple_receptors), value = T),
  assay = "SCT",
  layer = "counts"
)

iOSN_multiple_receptors$Olfr_number <- rowSums(iOSN_Olfr_counts > 0)

Idents(iOSN_multiple_receptors) <- "orig.ident"
Het_iOSN_multiple_receptors <- subset(iOSN_multiple_receptors, idents = "Het-1")

Taar_and_Olfr_coexpression_df_upset <- data.frame(
  Taar = Het_iOSN_multiple_receptors$Taar_number > 0,
  Olfr = Het_iOSN_multiple_receptors$Olfr_number > 0
)

pdf("seurat_V5_analysis/plots/upset_plot_Taar_Olfr_coexpression_in_Het_multiple_receptors_iOSN.pdf",width = 8, height = 6)
upset(
  Taar_and_Olfr_coexpression_df_upset,
  c("Taar", "Olfr"),
  name = "Cells per combination"
) + ggtitle("Taar & Olfr Co-expression in Het iOSNs")
dev.off()

#### Taar/OR coexpression in EGFP_Cre_iOSN expressing multiple receptors (EGFP_Cre_iOSN) (UMI > 0) ####
EGFP_Cre_iOSN_multiple_receptors <- subset(EGFP_Cre_iOSN, subset = receptor_number > 1)
table(EGFP_Cre_iOSN_multiple_receptors$orig.ident)
# Het-1  Homo 
# 207    69 

# Taar
EGFP_Cre_iOSN_Taar_counts <- FetchData(
  EGFP_Cre_iOSN_multiple_receptors,
  vars = grep("Taar", rownames(EGFP_Cre_iOSN_multiple_receptors), value = T),
  assay = "SCT",
  layer = "counts"
)

EGFP_Cre_iOSN_multiple_receptors$Taar_number <- rowSums(EGFP_Cre_iOSN_Taar_counts > 0)

# Olfr
EGFP_Cre_iOSN_Olfr_counts <- FetchData(
  EGFP_Cre_iOSN_multiple_receptors,
  vars = grep("Olfr", rownames(EGFP_Cre_iOSN_multiple_receptors), value = T),
  assay = "SCT",
  layer = "counts"
)

EGFP_Cre_iOSN_multiple_receptors$Olfr_number <- rowSums(EGFP_Cre_iOSN_Olfr_counts > 0)

##### Het #####
Idents(EGFP_Cre_iOSN_multiple_receptors) <- "orig.ident"
Het_EGFP_Cre_iOSN_multiple_receptors <- subset(EGFP_Cre_iOSN_multiple_receptors, idents = "Het-1")

df_upset <- data.frame(
  Taar = Het_EGFP_Cre_iOSN_multiple_receptors$Taar_number > 0,
  Olfr = Het_EGFP_Cre_iOSN_multiple_receptors$Olfr_number > 0
)

table(df_upset$Taar,df_upset$Olfr)
# FALSE TRUE
# FALSE     0  126
# TRUE      1   80
pdf("seurat_V5_analysis/plots/upset_plot_Taar_Olfr_coexpression_in_Het_EGFP_Cre_multiple_receptors_iOSN.pdf",width = 8, height = 6)
upset(
  df_upset,
  c("Taar", "Olfr"),
  name = "Cells per combination"
) + ggtitle("Taar & Olfr Co-expression in Het EGFP/Cre+ iOSNs")
dev.off()

##### Homo #####
Idents(EGFP_Cre_iOSN_multiple_receptors) <- "orig.ident"
Homo_EGFP_Cre_iOSN_multiple_receptors <- subset(EGFP_Cre_iOSN_multiple_receptors, idents = "Homo")

df_upset <- data.frame(
  Taar = Homo_EGFP_Cre_iOSN_multiple_receptors$Taar_number > 0,
  Olfr = Homo_EGFP_Cre_iOSN_multiple_receptors$Olfr_number > 0
)

pdf("seurat_V5_analysis/plots/upset_plot_Taar_Olfr_coexpression_in_Homo_EGFP_Cre_multiple_receptors_iOSN.pdf",width = 8, height = 6)
upset(
  df_upset,
  c("Taar", "Olfr"),
  name = "Cells per combination"
) + ggtitle("Taar & Olfr Co-expression in Homo EGFP/Cre+ iOSNs")
dev.off()

#### Taar/Class I OR/Class II OR coexpression in EGFP_Cre_iOSN expressing multiple receptors (EGFP_Cre_iOSN) (UMI > 0) ####
# Taar
EGFP_Cre_iOSN_Taar_counts <- FetchData(
  EGFP_Cre_iOSN_multiple_receptors,
  vars = grep("Taar", rownames(EGFP_Cre_iOSN_multiple_receptors), value = T),
  assay = "SCT",
  layer = "counts"
)
EGFP_Cre_iOSN_multiple_receptors$Taar_number <- rowSums(EGFP_Cre_iOSN_Taar_counts > 0)

# Class I OR
EGFP_Cre_iOSN_Class_I_OR_counts <- FetchData(
  EGFP_Cre_iOSN_multiple_receptors,
  vars = intersect(ClassI_OR_name, rownames(EGFP_Cre_iOSN_multiple_receptors)),
  assay = "SCT",
  layer = "counts"
)
EGFP_Cre_iOSN_multiple_receptors$Class_I_OR <- rowSums(EGFP_Cre_iOSN_Class_I_OR_counts > 0)

# Class II OR
EGFP_Cre_iOSN_Class_II_OR_counts <- FetchData(
  EGFP_Cre_iOSN_multiple_receptors,
  vars = intersect(ClassII_OR_name, rownames(EGFP_Cre_iOSN_multiple_receptors)),
  assay = "SCT",
  layer = "counts"
)
EGFP_Cre_iOSN_multiple_receptors$Class_II_OR <- rowSums(EGFP_Cre_iOSN_Class_II_OR_counts > 0)

##### Het #####
Idents(EGFP_Cre_iOSN_multiple_receptors) <- "orig.ident"
Het_EGFP_Cre_iOSN_multiple_receptors <- subset(EGFP_Cre_iOSN_multiple_receptors, idents = "Het-1")
# An object of class Seurat 
# 53020 features across 207 samples within 2 assays

df_upset <- data.frame(
  Taar = Het_EGFP_Cre_iOSN_multiple_receptors$Taar_number > 0,
  Class_I_OR = Het_EGFP_Cre_iOSN_multiple_receptors$Class_I_OR > 0,
  Class_II_OR = Het_EGFP_Cre_iOSN_multiple_receptors$Class_II_OR > 0
)

Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors <- data.frame(
  Taar_coexpression = sum(df_upset$Taar & !(df_upset$Class_I_OR) & !(df_upset$Class_II_OR)),
  Class_I_OR_coexpression = sum(!(df_upset$Taar) & df_upset$Class_I_OR & !(df_upset$Class_II_OR)),
  Class_II_OR_coexpression = sum(!(df_upset$Taar) & !(df_upset$Class_I_OR) & df_upset$Class_II_OR),
  Taar_ClassIOR_coexpression = sum(df_upset$Taar & df_upset$Class_I_OR & !(df_upset$Class_II_OR)),
  Taar_ClassIIOR_coexpression = sum(df_upset$Taar & !(df_upset$Class_I_OR) & df_upset$Class_II_OR),
  ClassIOR_ClassIIOR_coexpression = sum(!(df_upset$Taar) & df_upset$Class_I_OR & df_upset$Class_II_OR),
  Taar_ClassIOR_ClassIIOR_coexpression = sum(df_upset$Taar & df_upset$Class_I_OR & df_upset$Class_II_OR)
) %>% t()
colnames(Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors) <- "cellNumber"
colSums(Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors)
Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors <- Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors %>% 
  as.data.frame() %>% 
  mutate(percent = cellNumber / nrow(df_upset))
write.csv(Taar_ClassIOR_ClassIIOR_Het_coexpression_EGFP_Cre_iOSN_multiple_receptors, 
          file = "seurat_V5_analysis/files/Taar_ClassIOR_ClassIIOR_coexpression_Het_EGFP_Cre_iOSN_multiple_receptors.csv")

pdf("seurat_V5_analysis/plots/upset_plot_Taar_classI_classII_OR_coexpression_in_Het_EGFP_Cre_multiple_receptors_iOSN.pdf",width = 8, height = 6)
upset(
  df_upset,
  c("Taar", "Class_I_OR", 'Class_II_OR'),
  name = "Cells per combination") + 
  ggtitle("Receptor Co-expression in Het EGFP/Cre+ iOSNs")
dev.off()

##### Homo #####
Idents(EGFP_Cre_iOSN_multiple_receptors) <- "orig.ident"
Homo_EGFP_Cre_iOSN_multiple_receptors <- subset(EGFP_Cre_iOSN_multiple_receptors, idents = "Homo")

df_upset <- data.frame(
  Taar = Homo_EGFP_Cre_iOSN_multiple_receptors$Taar_number > 0,
  Class_I_OR = Homo_EGFP_Cre_iOSN_multiple_receptors$Class_I_OR > 0,
  Class_II_OR = Homo_EGFP_Cre_iOSN_multiple_receptors$Class_II_OR > 0
)

Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors <- data.frame(
  Taar_coexpression = sum(df_upset$Taar & !(df_upset$Class_I_OR) & !(df_upset$Class_II_OR)),
  Class_I_OR_coexpression = sum(!(df_upset$Taar) & df_upset$Class_I_OR & !(df_upset$Class_II_OR)),
  Class_II_OR_coexpression = sum(!(df_upset$Taar) & !(df_upset$Class_I_OR) & df_upset$Class_II_OR),
  Taar_ClassIOR_coexpression = sum(df_upset$Taar & df_upset$Class_I_OR & !(df_upset$Class_II_OR)),
  Taar_ClassIIOR_coexpression = sum(df_upset$Taar & !(df_upset$Class_I_OR) & df_upset$Class_II_OR),
  ClassIOR_ClassIIOR_coexpression = sum(!(df_upset$Taar) & df_upset$Class_I_OR & df_upset$Class_II_OR),
  Taar_ClassIOR_ClassIIOR_coexpression = sum(df_upset$Taar & df_upset$Class_I_OR & df_upset$Class_II_OR)
) %>% t()
colnames(Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors) <- "cellNumber"
colSums(Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors)
Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors <- Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors %>% 
  as.data.frame() %>% 
  mutate(percent = cellNumber / nrow(df_upset))
write.csv(Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors, 
          file = "seurat_V5_analysis/files/Taar_ClassIOR_ClassIIOR_coexpression_Homo_EGFP_Cre_iOSN_multiple_receptors.csv")
sum(Taar_ClassIOR_ClassIIOR_Homo_coexpression_EGFP_Cre_iOSN_multiple_receptors$percent)

pdf("seurat_V5_analysis/plots/upset_plot_Taar_classI_classII_OR_coexpression_in_Homo_EGFP_Cre_multiple_receptors_iOSN.pdf",width = 8, height = 6)
upset(
  df_upset,
  c("Taar", "Class_I_OR", 'Class_II_OR'),
  name = "Cells per combination") + 
  ggtitle("Receptor Co-expression in Homo GFP/Cre+ iOSNs")
dev.off()

#### cell number per cell type (EGFP_Cre+ cells) ####
EGFP_Cre_cells <- subset(CCAIntegration_obj, subset = EGFP_or_Cre_marker_sum > 0)
Idents(EGFP_Cre_cells) <- "celltype"
celltype_summary_by_group <- table(EGFP_Cre_cells$celltype, EGFP_Cre_cells$orig.ident)
# Het-1 Homo
# mOSN                          113  101
# iOSN                          269  302
# late INP                       11   46
# early INP                       3    2
# GBC                             4    3
# HBC                             3    0
# SUS                             4    9
# MV                              1    1
# OEC                             0    1
# vasculature cell                3    3
# immune cell                     0    1
# respiratory epithelial cell     0    3
# immature excitatory neuron    175   83
# immature interneuron            2    0

# cell count
celltype_summary_by_group_df <- as.data.frame(celltype_summary_by_group)
colnames(celltype_summary_by_group_df) <- c("celltype", "group", "Cell_number")

pdf("seurat_V5_analysis/plots/barplot_celltype_class_in_EGFP_Cre_cells.pdf", width = 6, height = 6)
ggplot(celltype_summary_by_group_df, aes(x = group, y = Cell_number, fill = as.factor(celltype))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(celltype_summary_by_group_df$celltype)))) +
  labs(
    x = "Group",
    y = "Number of cells",
    fill = "celltype"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggplot(celltype_summary_by_group_df, aes(x = group, y = Cell_number, fill = as.factor(celltype))) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(celltype_summary_by_group_df$celltype)))) +
  labs(
    x = "Group",
    y = "Percentage of cells (%)",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# percent
celltype_summary_by_group_df <- celltype_summary_by_group_df %>%
  group_by(group) %>%
  mutate(Percent = Cell_number / sum(Cell_number) )

ggplot(celltype_summary_by_group_df, aes(x = celltype, y = Percent, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 并列显示不同 group
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(
    x = "Cell type",
    y = "Proportion of cells (%)",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

#### EGFP/Cre & Olfr or Taar coexpression (OSN-lineage) ####
gene_info <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T)
all_functional_olfactory_receptors <- gene_info[which(grepl("Taar|Olfr", gene_info$gene_name) & gene_info$gene_type == "protein_coding" & gene_info$annotation_source == "HAVANA"),]

mOSN <- subset(CCAIntegration_obj, idents = "mOSN")
receptor_genes <- intersect(
  all_functional_olfactory_receptors$gene_name,
  rownames(mOSN)
)

TAAR_info <- all_functional_olfactory_receptors[grep("Taar", all_functional_olfactory_receptors$gene_name),]
Taar_name <- TAAR_info$gene_name

ClassI_ClassII_OR_info <- read.csv("/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
rownames(ClassI_ClassII_OR_info) <- gsub("\\.ps", "-ps",  rownames(ClassI_ClassII_OR_info))

ClassI_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_I" & ClassI_ClassII_OR_info$gene_type == "protein_coding")]
length(ClassI_OR_name)

ClassII_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_II" & ClassI_ClassII_OR_info$gene_type == "protein_coding")]
length(ClassII_OR_name)

receptor_class_sets <- list(
  Taar_genes = Taar_name,
  ClassI_OR_genes = ClassI_OR_name,
  ClassII_OR_genes = ClassII_OR_name
)

### extract OSN lineage 
CCAIntegration_obj_data <- GetAssayData(CCAIntegration_obj, assay = "SCT", layer  = "data") 

for (name in names(receptor_class_sets)) {
  genes <- intersect(receptor_class_sets[[name]], rownames(CCAIntegration_obj_data))  # 确保存在
  CCAIntegration_obj[[paste0(name, "_data_sum")]] <- Matrix::colSums(CCAIntegration_obj_data[genes, , drop = FALSE])
}

Idents(CCAIntegration_obj) <- "celltype"

OSN_lineage <- subset(
  CCAIntegration_obj,
  idents = c("mOSN", "iOSN", "late INP", "early INP",
             "GBC", "HBC")
)

### receptor expression without filter
OSN_lineage$Taar_genes_filter_data_sum <-
  OSN_lineage$Taar_genes_data_sum

OSN_lineage$ClassI_OR_genes_filter_data_sum <-
  OSN_lineage$ClassI_OR_genes_data_sum

OSN_lineage$ClassII_OR_genes_filter_data_sum <-
  OSN_lineage$ClassII_OR_genes_data_sum

### extract mOSN (filter using UMI > 1)
mOSN <- subset(OSN_lineage, idents = "mOSN")

mOSN_counts <- FetchData(
  mOSN,
  vars = receptor_genes,
  assay = "SCT",
  layer = "counts"
)

mOSN_data <- FetchData(
  mOSN,
  vars = receptor_genes,
  assay = "SCT",
  layer = "data"
)

# counts ≤ 1 的基因置 0
mOSN_data[mOSN_counts <= 1] <- 0

receptor_class_sets <- list(
  Taar_genes = Taar_name,
  ClassI_OR_genes = ClassI_OR_name,
  ClassII_OR_genes = ClassII_OR_name
)

# 将UMI>1过滤后计算好的data sum写入mOSN
receptor_sums_mOSN <- lapply(receptor_class_sets, function(gene_set) {
  genes <- intersect(gene_set, colnames(mOSN_data))
  Matrix::rowSums(mOSN_data[, genes, drop = FALSE])
})

### change the mOSN data to filtered data
mOSN_cells <- colnames(mOSN)
all(mOSN_cells == names(receptor_sums_mOSN[["Taar_genes"]]))

for (name in names(receptor_sums_mOSN)) {
  col_name <- paste0(name, "_filter_data_sum")
  
  # 直接替换 mOSN_cells 的值，其它细胞保持原值
  OSN_lineage@meta.data[mOSN_cells, col_name] <- receptor_sums_mOSN[[name]]
}

# coexpression
OSN_lineage$EGFP_Taar_filter_coexpression <- OSN_lineage$Taar_genes_filter_data_sum > 0 & OSN_lineage$EGFP_or_Cre_marker_sum > 0

OSN_lineage$OR_genes_filter_data_sum <- OSN_lineage$ClassI_OR_genes_filter_data_sum + OSN_lineage$ClassII_OR_genes_filter_data_sum

OSN_lineage$EGFP_Olfr_filter_coexpression <-
  OSN_lineage$OR_genes_filter_data_sum > 0 & OSN_lineage$EGFP_or_Cre_marker_sum > 0

### checke receptor expression in EGFP/Cre+ OSNs
EGFP_Cre_OSN <- subset(OSN_lineage, subset = EGFP_or_Cre_marker_sum > 0)

table(EGFP_Cre_OSN$EGFP_Taar_filter_coexpression,EGFP_Cre_OSN$orig.ident)
# Het-1 Homo
# FALSE   228  447
# TRUE    175    7

table(EGFP_Cre_OSN$EGFP_Olfr_filter_coexpression,EGFP_Cre_OSN$orig.ident)
# Het-1 Homo
# FALSE   119  182
# TRUE    284  272

OSN_lineage$orig.ident <- factor(OSN_lineage$orig.ident, levels = c("Het-1", "Homo"))

save(OSN_lineage, file = "seurat_V5_analysis/rdata/OSN_lineage_coexpression.rdata")

### mOSN profile
# devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

df <- as.data.frame(Embeddings(OSN_lineage, "umap.CCA"))
df$cluster <- Idents(OSN_lineage)
# 如果用该阈值，最后看起来Homo组中mOSN内有Taar和Tbr1共表达的细胞
# core_df <- df %>%
#   filter(cluster == "mOSN") %>%
#   filter(
#     between(
#       umapCCA_1,
#       quantile(umapCCA_1, 0.025),
#       quantile(umapCCA_1, 0.975)
#     ),
#     between(
#       umapCCA_2,
#       quantile(umapCCA_2, 0.025),
#       quantile(umapCCA_2, 0.975)
#     )
#   )

# 改成该阈值
core_df <- df %>%
  filter(cluster == "mOSN") %>%
  filter(
    between(
      umapCCA_1,
      quantile(umapCCA_1, 0.2),
      quantile(umapCCA_1, 0.8)
    ),
    between(
      umapCCA_2,
      quantile(umapCCA_2, 0.2),
      quantile(umapCCA_2, 0.8)
    )
  )

core_df$orig.ident <- OSN_lineage$orig.ident[
  match(rownames(core_df), colnames(OSN_lineage))
]

core_df$orig.ident <- factor(
  core_df$orig.ident,
  levels = levels(OSN_lineage$orig.ident)
)

# color
barplot(rep(1,7), col = c("#33a02c","#4575b4","#023858","#dd3497","#7a0177","#ec7014","#993404"),
        legend.text = c("EGFP", "Taar", "EGFP_Taar_co",
                        "OR", "EGFP_OR_co",
                        "Gucy1b2", "EGFP_Gucy1b2_co"))
# plot

pdf("seurat_V5_analysis/plots/featurePlot_EGFP_Olfr_coexpression_filtered_data_with_mOSN_profile_in_OSN_lineage_order_small_point_0.1_max_cutoff_15_small_circle.pdf", width = 8, height = 4)

p1 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "EGFP_or_Cre_marker_sum",
  order = TRUE,
  cols = c("gray90","#33a02c"),
  pt.size = 0.1
) +
  labs(color = "EGFP_or_Cre_marker_sum (order = T)")

p1 +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05,
    colour = "black"
  )

p2 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "EGFP_or_Cre_marker_sum",
  order = TRUE,
  cols = c("gray90","#33a02c"),
  pt.size = 0.1,
  split.by = "orig.ident",
  keep.scale = "all"   # 加这一行！
) +
  labs(color = "EGFP_or_Cre_marker_sum (order = T)")

p2[[1]] <- p2[[1]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p2[[2]] <- p2[[2]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p2

p3 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "OR_genes_filter_data_sum",
  order = TRUE,
  cols = c("gray90","#dd3497"),
  pt.size = 0.1,
  max.cutoff = 15
) +
  labs(color = "OR_genes_filter_data_sum (order = T)")

p3 +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05,
    colour = "black"
  )

p4 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "OR_genes_filter_data_sum",
  order = TRUE,
  cols = c("gray90","#dd3497"),
  pt.size = 0.1,
  max.cutoff = 15,
  split.by = "orig.ident",
  keep.scale = "all"   # 加这一行！
) +
  labs(color = "OR_genes_filter_data_sum (order = T)")

p4[[1]] <- p4[[1]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p4[[2]] <- p4[[2]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p4

p5 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "EGFP_Olfr_filter_coexpression",
  order = TRUE,
  cols = c("gray90","#7a0177"),
  pt.size = 0.1
) +
  labs(color = "EGFP_Olfr_filter_coexpression (order = T)")

p5 +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05,
    colour = "black"
  )

p6 <- FeaturePlot(
  OSN_lineage,
  reduction = "umap.CCA",
  features = "EGFP_Olfr_filter_coexpression",
  order = TRUE,
  cols = c("gray90","#7a0177"),
  pt.size = 0.1,
  split.by = "orig.ident",
  keep.scale = "all"   # 加这一行！
) +
  labs(color = "EGFP_Olfr_filter_coexpression (order = T)")

p6[[1]] <- p6[[1]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p6[[2]] <- p6[[2]] +
  geom_encircle(
    data = subset(core_df, cluster == "mOSN"),
    aes(umapCCA_1, umapCCA_2),
    linetype = "dashed",
    expand = 0.05
  )

p6

dev.off()






