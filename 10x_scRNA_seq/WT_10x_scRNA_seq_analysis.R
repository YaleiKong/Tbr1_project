rm(list=ls())

## Load required packages
.libPaths(new = "/lustre/home/acct-medlqian/medlqian-loop3/R/x86_64-redhat-linux-gnu-library/4.2.2")
library(Seurat) # 文章中用的v4.1.0
library(ggplot2)
library(dplyr)
library(purrr) 
library(tidyverse) 
library(stringr)
library(RColorBrewer)
library(patchwork)
library(colorspace) # darken
library(harmony)
library(ggvenn)
library(ggplot2)
library(Matrix)
library(clustree)
library(ComplexUpset)
library(openxlsx)

set.seed(24)
setwd("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_10X_scRNA_Seq_seurat_analysis")

## Creat output directory 
out_dir <- "/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_10X_scRNA_Seq_seurat_analysis/integrated_analysis_seuratV5/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load the dataset
GSE120199_samples <- c("GSE120199_P0_rep1", "GSE120199_P0_rep2", "GSE120199_P3_rep1", "GSE120199_P3_rep2", 
                       "GSE120199_P7_rep1", "GSE120199_P7_rep2", "GSE120199_P21_rep1", "GSE120199_P21_rep2")
GSE120199_scRNA_dirs <- str_c("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_sc_seurat_analysis/GSE120199_C.RonYu_Neuron_2018/raw_data/", GSE120199_samples)

GSE147459_samples <- c("GSE147459_E18.5","GSE147459_P14","GSE147459_adult")
GSE147459_scRNA_dirs <- str_c("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_sc_seurat_analysis/GSE147459_C.RonYu_bioRxiv_2022/raw_data/", GSE147459_samples)

GSE173947_samples <- paste("GSE173947_homecage",1:6,sep="_")
GSE173947_scRNA_dirs <- str_c("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_sc_seurat_analysis/GSE173947_SandeepRobertDatta_Cell_2021/raw_data/", GSE173947_samples)

GSE224604_samples <- c("GSE224604_19d_rep1","GSE224604_19d_rep2","GSE224604_8w_rep1","GSE224604_8w_rep2")
GSE224604_scRNA_dirs <- str_c("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/20250518_WT_sc_seurat_analysis/GSE224604_JinXu_cellReports_2024/scRNA/cellranger/", GSE224604_samples)

scRNA_dirs <- c(GSE120199_scRNA_dirs, GSE147459_scRNA_dirs, GSE173947_scRNA_dirs, GSE224604_scRNA_dirs)
samples <- c(GSE120199_samples, GSE147459_samples, GSE173947_samples, GSE224604_samples)

objList <- lapply(1:length(scRNA_dirs),function(i){
  CreateSeuratObject(counts = Read10X(scRNA_dirs[i]),
                     project = samples[i],assay = "RNA") 
})

for ( i in 1:length(objList)){
  objList[[i]][["percent.mt"]] = PercentageFeatureSet(objList[[i]], pattern="^mt-")
  pdf(str_c(out_dir,"plots/", samples[i],"_qc_plot.pdf"),width=9)
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

cell_counts_by_condition <- lapply(objList, function(obj) {
  obj[["percent.hb"]] = PercentageFeatureSet(obj, pattern="^Hb[ab]")
  data.frame(
    Sample = obj@project.name,
    Total = ncol(obj),
    nFeature_200_plus     = sum(obj$nFeature_RNA >= 200),
    nFeature_400_plus     = sum(obj$nFeature_RNA >= 400),
    nFeature_5000_minus   = sum(obj$nFeature_RNA <= 5000),
    nCount_700_plus       = sum(obj$nCount_RNA >= 700),
    nCount_15000_minus    = sum(obj$nCount_RNA <= 15000),
    percent_mt_under_10   = sum(obj$percent.mt < 10),
    percent_hb_under_1 = sum(obj$percent.hb < 1),
    percent_hb_under_5 = sum(obj$percent.hb < 5)
  )
})

cell_counts_by_condition_df <- do.call(rbind, cell_counts_by_condition)
print(cell_counts_by_condition_df)
# Sample Total nFeature_200_plus nFeature_400_plus nFeature_5000_minus nCount_700_plus
# 1     GSE120199_P0_rep1  2761              2761              2704                2760            2623
# 2     GSE120199_P0_rep2  2769              2769              2709                2768            2636
# 3     GSE120199_P3_rep1  3923              3923              3787                3923            3659
# 4     GSE120199_P3_rep2  3951              3951              3811                3951            3684
# 5     GSE120199_P7_rep1  2061              2061              1987                2060            1951
# 6     GSE120199_P7_rep2  2068              2068              1997                2067            1960
# 7    GSE120199_P21_rep1  1574              1573              1530                1574            1515
# 8    GSE120199_P21_rep2  1580              1579              1540                1580            1522
# 9       GSE147459_E18.5  3790              3785              3737                3770            3703
# 10        GSE147459_P14  3532              3532              3497                3528            3280
# 11      GSE147459_adult  3861              3861              3383                3861            2992
# 12 GSE173947_homecage_1 13702             13398             13189               11445           13424
# 13 GSE173947_homecage_2 10900             10703             10462                8759           10766
# 14 GSE173947_homecage_3 11053             10791             10490                8913           10884
# 15 GSE173947_homecage_4 10349             10230              9962                8121           10248
# 16 GSE173947_homecage_5 12998             12837             12652               11167           12853
# 17 GSE173947_homecage_6 12293             12117             11957               10273           12171
# 18   GSE224604_19d_rep1 10450              9731              9334                8861           10369
# 19   GSE224604_19d_rep2 16270             15571             15026               15985           15800
# 20    GSE224604_8w_rep1  5351              5343              5246                5039            5169
# 21    GSE224604_8w_rep2 12183             12165             11986               12054           11917
# nCount_15000_minus percent_mt_under_10 percent_hb_under_1 percent_hb_under_5
# 1                2760                2618               2759               2759
# 2                2768                2625               2767               2767
# 3                3921                3545               3923               3923
# 4                3949                3553               3951               3951
# 5                2057                1757               2061               2061
# 6                2064                1758               2068               2068
# 7                1569                1435               1572               1572
# 8                1574                1439               1578               1578
# 9                3635                3415               3784               3786
# 10               3506                2914               3532               3532
# 11               3851                3276               3861               3861
# 12               8070               12675              13691              13691
# 13               5329                9955              10877              10877
# 14               5141                9914              11027              11028
# 15               3961                9451              10340              10340
# 16               8668               11975              12986              12987
# 17               7425               11388              12280              12280
# 18               8105                8030              10179              10188
# 19              15542               15519              13687              14299
# 20               4662                5179               5172               5288
# 21              11855               11430              12146              12156

write.csv(cell_counts_by_condition_df, file = "integrated_analysis_seuratV5/files/cell_counts_by_condition_df.csv")

for (i in 1:length(objList)){
  objList[[i]] <- subset(objList[[i]],subset = nFeature_RNA >= 400 & 
                           nCount_RNA >= 700  & percent.mt < 10)
}

merged_obj <- merge(objList[[1]],
                    y=c(objList[[2]],objList[[3]],objList[[4]],objList[[5]],objList[[6]],objList[[7]],
                        objList[[8]],objList[[9]],objList[[10]],objList[[11]],objList[[12]],objList[[13]],
                        objList[[14]],objList[[15]],objList[[16]],objList[[17]],objList[[18]],objList[[19]],
                        objList[[20]],objList[[21]]),
                    add.cell.ids = samples)

merged_obj <- SCTransform(merged_obj, vars.to.regress = c("percent.mt")) %>% 
  RunPCA(verbose = FALSE)

# HarmonyIntegration
HarmonyIntegration_obj <- IntegrateLayers(
  object = merged_obj,
  method = HarmonyIntegration,
  assay = "SCT",
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.Harmony",
  verbose = TRUE,
  .options = harmony_options(
    theta = 2,  # Diversity penalty (adjust if over/under-correction observed)
    mmax_iter = 50,  # Increase iterations for large datasets
    nclust = NULL,  # Let Harmony estimate number of clusters
    early_stop = TRUE  # Allow early stopping for convergence
  )
)

HarmonyIntegration_obj <- RunUMAP(HarmonyIntegration_obj,  reduction = "integrated.Harmony", reduction.name = "umap.Harmony", dims=1:30)

#### tune the PCs ####
# # tune the number of PCs
# markers <- c( # OE
#   "Syt1", "Syp", "Snap25", "Vamp2", # neurons
#   "Nqo1","Ncam2", # dorsal/ventral OSNs
#   "Omp", "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1", # Mature OSNs
#   "Gap43","Gng8", "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", "Dpysl3", "Stmn1", "Stmn2", "Ebf2", "Cbx8", "Trib3", #Immature OSNs
#   "Sox11","Neurog1", "Neurod1","Neurod2", "Top2a", "Mki67", "Lhx2", "Ebf1", #INPs, immediate neuronal precursors 
#   "Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a",#GBCs
#   "Krt5","Krt14", "Trp63", "Cxcl14", "Meg3",#HBCs
#   "Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3", # Ms4a-expressing chemosensory receptor cell (Ms4)
#   "Sox2","Ermn","Cyp2g1","Cyp1a2", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Hey1", "Sult1c1", # Sustentacular cells
#   "Coch", "Ascl3", "Cftr", "Foxi1","Hepacam2", # olfactory microvillar cells (Trpm5-)
#   "Sox9", "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Ltc4s", # mv (Trpm5+)
#   
#   # 固有层 (lamina propria)
#   "Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3", # Ensheathing glia
#   "Sox9","Sox10", "Krt18", "Muc5ac", "Gpx3", # Bowman's gland
#   "Eng","Sox17",# Pericytes # Vascular wall cells
#   "Tagln", "Myh11", # vascular smooth muscle cell # Vascular wall cells
#   "Pecam1", "Kdr", # Endothelial cells # Vascular wall cells
#   
#   # respiratory cells
#   "Foxj1", "Cfap126", "Stoml3", # respiratory ciliated cells 
#   "Sox9",  # respiratory gland progenitor cells 
#   "Cyp4b1", "Muc5ac", "Tff3", # respiratory secretory cells
#   
#   # blood and immune cells
#   "Cd3d", "Cd3e", "Cd8a", # cd8+ t cell
#   "Cd3d", "Cd3e", "Cd4", "Il7r", # cd4+ t cell
#   "Cx3cr1", # natural killer cell
#   "Cd37", "Cd19", "Cd79a", "Ms4a4a", # b cell
#   "Mzb1", "Sdc1", "Cd79a", # plasma cell
#   "Hbb-bs","Hbb-bt",# Erythrocytes
#   "S100a9","S100a8",# Neutrophils
#   "Cd14", "Clec10a", "Lyz2", "S100a4", # monocyte
#   "C1qa", "C1qb", "C1qc", "Ms4a7", # macrophage
#   "Hba-a1", #blood cell
#   
#   "Col1a1","Bglap" #Osteogenic cells
# )
# 
# markers <- markers[markers %in% rownames(HarmonyIntegration_obj)]
# all(markers %in% rownames(HarmonyIntegration_obj))
# unique_markers <- unique(markers)
# 
# marker_sets <- list(
#   mOSN = c("Omp", "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1"),
#   imOSN = c("Gap43","Gng8", "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", "Dpysl3", "Stmn1", "Stmn2", "Ebf2", "Cbx8", "Trib3"),
#   INP = c("Sox11","Neurog1", "Neurod1","Neurod2", "Top2a", "Mki67", "Lhx2", "Ebf1"),
#   GBC = c("Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a"),
#   HBC = c("Krt5","Krt14", "Trp63", "Cxcl14", "Meg3"),
#   Ms4a = c("Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3"),
#   SUS = c("Sox2","Ermn","Cyp2g1","Cyp1a2", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Hey1", "Sult1c1"),
#   Ensheathing = c("Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3"),
#   mv1 = c("Coch", "Ascl3", "Cftr", "Hepacam2"),
#   mv2 = c("Sox9", "Trpm5"),
#   blood_immune = c(
#     "Cd3d", "Cd3e", "Cd8a", "Cd4", "Il7r", "Cx3cr1",
#     "Cd37", "Cd19", "Cd79a", "Ms4a4a", "Mzb1", "Sdc1",
#     "Hbb-bs","Hbb-bt", "S100a9","S100a8", "Cd14", "Clec10a",
#     "Lyz2", "S100a4", "C1qa", "C1qb", "C1qc", "Ms4a7", "Hba-a1"
#   )
# )
# 
# mat <- GetAssayData(HarmonyIntegration_obj, assay = "SCT", layer = "data")  # 通常用 "data" slot
# 
# # 加入 meta.data
# for (name in names(marker_sets)) {
#   genes <- intersect(marker_sets[[name]], rownames(mat))  # 确保存在
#   HarmonyIntegration_obj[[paste0(name, "_marker_sum")]] <- Matrix::colSums(mat[genes, , drop = FALSE])
# }
# 
# for (nPCs in seq(20,40,5)){
#   HarmonyIntegration_obj <- FindNeighbors(HarmonyIntegration_obj, dims = 1:nPCs, reduction = "integrated.Harmony")
#   HarmonyIntegration_obj <- FindClusters(HarmonyIntegration_obj, resolution = 0.2, reduction = "integrated.Harmony")
#   HarmonyIntegration_obj <- RunUMAP(HarmonyIntegration_obj, dims = 1:nPCs, reduction = "integrated.Harmony", reduction.name = "umap")
#   
#   pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
#   p1 <- DimPlot(HarmonyIntegration_obj, reduction = "umap", label=TRUE, raster=FALSE) + theme(legend.position = "bottom")
#   p2 <- DimPlot(HarmonyIntegration_obj, reduction = "umap", group.by="orig.ident", raster=FALSE) + theme(legend.position = "bottom")
#   print(p1)
#   print(p2)
#   dev.off()
#   
#   DefaultAssay(HarmonyIntegration_obj) <- "SCT"
#   pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
#   p1 <- DimPlot(HarmonyIntegration_obj, reduction = "umap", label=TRUE, raster=FALSE) + theme(legend.position = "bottom")
#   p2 <- DimPlot(HarmonyIntegration_obj, reduction = "umap", group.by="orig.ident", raster=FALSE) + theme(legend.position = "bottom")
#   print(p1)
#   print(p2)
#   dev.off()
#   
#   DefaultAssay(HarmonyIntegration_obj) <- "SCT"
#   pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",nPCs,"_markers_sum_UMAP.pdf"))
#   p <- FeaturePlot(HarmonyIntegration_obj, reduction = "umap", features=paste0(names(marker_sets),"_marker_sum"), combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
#   print(p)
#   dev.off()
# }

#### tune the resolution (40 PCs) ####
HarmonyIntegration_obj <- FindNeighbors(HarmonyIntegration_obj, reduction = "integrated.Harmony", dims = 1:40)
HarmonyIntegration_obj <- RunUMAP(HarmonyIntegration_obj, dims = 1:40, reduction = "integrated.Harmony", reduction.name = "umap.Harmony")

# marker
markers <- c( # OE
  "Syt1", "Syp", "Snap25", "Vamp2", # neurons
  "Nqo1","Ncam2", # dorsal/ventral OSNs
  "Omp", "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1", # Mature OSNs
  "Gap43","Gng8", "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", "Dpysl3", "Stmn1", "Stmn2", "Ebf2", "Cbx8", "Trib3", #Immature OSNs
  "Neurog1", "Neurod1", "Tex15", "Sox11","Neurod2", "Top2a", "Mki67", "Lhx2", "Ebf1", #INPs, immediate neuronal precursors
  "Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a",#GBCs
  "Krt5","Krt14", "Trp63", "Cxcl14", "Meg3",#HBCs
  "Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3", # Ms4a-expressing chemosensory receptor cell (Ms4)
  "Sox2","Ermn","Cyp2g1","Cyp1a2", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Hey1", "Sult1c1", # Sustentacular cells
  "Coch", "Ascl3", "Cftr", "Foxi1","Hepacam2", # olfactory microvillar cells (Trpm5-)
  "Sox9", "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Ltc4s", # mv (Trpm5+)
  
  # 固有层 (lamina propria)
  "Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3", # Ensheathing glia
  "Sox9","Sox10", "Krt18", "Muc5ac", "Gpx3", # Bowman's gland
  "Eng","Sox17",# Pericytes # Vascular wall cells
  "Tagln", "Myh11", # vascular smooth muscle cell # Vascular wall cells
  "Pecam1", "Kdr", # Endothelial cells # Vascular wall cells
  
  # respiratory cells
  "Foxj1", "Cfap126", "Stoml3", # respiratory ciliated cells
  "Sox9",  # respiratory gland progenitor cells
  "Cyp4b1", "Muc5ac", "Tff3", # respiratory secretory cells
  
  # blood and immune cells
  "Cd3d", "Cd3e", "Cd8a", # cd8+ t cell
  "Cd3d", "Cd3e", "Cd4", "Il7r", # cd4+ t cell
  "Cx3cr1", # natural killer cell
  "Cd37", "Cd19", "Cd79a", "Ms4a4a", # b cell
  "Mzb1", "Sdc1", "Cd79a", # plasma cell
  "Hbb-bs","Hbb-bt",# Erythrocytes
  "S100a9","S100a8",# Neutrophils
  "Cd14", "Clec10a", "Lyz2", "S100a4", # monocyte
  "C1qa", "C1qb", "C1qc", "Ms4a7", # macrophage
  "Hba-a1", #blood cell
  
  "Col1a1","Bglap" #Osteogenic cells
)
markers <- markers[markers %in% rownames(HarmonyIntegration_obj)]
all(markers %in% rownames(HarmonyIntegration_obj))
unique_markers <- unique(markers)

pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",40,"_markers_UMAP_order.pdf"))
p <- FeaturePlot(HarmonyIntegration_obj, reduction = "umap.Harmony", features=unique_markers, combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
print(p)
dev.off()

# marker sum
marker_sets <- list(
  mOSN = c("Omp", "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1"),
  imOSN = c("Gap43","Gng8", "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", "Dpysl3", "Stmn1", "Stmn2", "Ebf2", "Cbx8", "Trib3"),
  INP = c("Neurog1", "Neurod1", "Tex15","Sox11","Neurod2", "Top2a", "Mki67", "Lhx2", "Ebf1"),
  GBC = c("Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a"),
  HBC = c("Krt5","Krt14", "Trp63", "Cxcl14", "Meg3"),
  Ms4a = c("Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3"),
  SUS = c("Sox2","Ermn","Cyp2g1","Cyp1a2", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Hey1", "Sult1c1"),
  Ensheathing = c("Atp1a2","Fabp7", "S100b", "Plp1", "Mpz", "Alx3"),
  mv1 = c("Coch", "Ascl3", "Cftr", "Hepacam2"),
  mv2 = c("Sox9", "Trpm5"),
  blood_immune = c(
    "Cd3d", "Cd3e", "Cd8a", "Cd4", "Il7r", "Cx3cr1",
    "Cd37", "Cd19", "Cd79a", "Ms4a4a", "Mzb1", "Sdc1",
    "Hbb-bs","Hbb-bt", "S100a9","S100a8", "Cd14", "Clec10a",
    "Lyz2", "S100a4", "C1qa", "C1qb", "C1qc", "Ms4a7", "Hba-a1"
  )
)

# 加入 meta.data
mat <- GetAssayData(HarmonyIntegration_obj, assay = "SCT", layer = "data")  # 通常用 "data" slot
for (name in names(marker_sets)) {
  genes <- intersect(marker_sets[[name]], rownames(mat))  # 确保存在
  HarmonyIntegration_obj[[paste0(name, "_marker_sum")]] <- Matrix::colSums(mat[genes, , drop = FALSE])
}

pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",40,"_markers_sum_UMAP_order.pdf"))
p <- FeaturePlot(HarmonyIntegration_obj, reduction = "umap.Harmony", features=paste0(names(marker_sets),"_marker_sum"), combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
print(p)
dev.off()

pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC",40,"_markers_sum_UMAP.pdf"))
p <- FeaturePlot(HarmonyIntegration_obj, reduction = "umap.Harmony", features=paste0(names(marker_sets),"_marker_sum"), 
                 combine=FALSE, cols=c("lightgrey", "red"), order=F, raster=FALSE)
print(p)
dev.off()

# 运行多个 resolution
resolutions <- seq(0.1, 1.5, by = 0.1)
for (res in resolutions) {
  HarmonyIntegration_obj <- FindClusters(HarmonyIntegration_obj, resolution = res)
  
  pdf(str_c(out_dir,"plots/integrated_OE_scRNA_PC40_resolution",res,"_UMAP.pdf"))
  p1 <- DimPlot(HarmonyIntegration_obj, reduction = "umap.Harmony", label=TRUE, group.by=str_c("SCT_snn_res.", res))
  p2 <- DimPlot(HarmonyIntegration_obj, reduction = "umap.Harmony", group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  
  pdf(width = 15, height = 80, str_c(out_dir,"plots/integrated_OE_scRNA_PC40_resolution", res,"_vlnplot_marker.pdf"))
  p1 <- VlnPlot(HarmonyIntegration_obj, group.by=str_c("SCT_snn_res.", res), features = unique_markers, pt.size = 0) & xlab(NULL) & ylab(NULL)
  print(p1)
  dev.off()
  
  pdf(width = 24, height = 20, str_c(out_dir,"plots/integrated_OE_scRNA_PC40_resolution",res,"_dotplot_marker.pdf"))
  p2 <- DotPlot(HarmonyIntegration_obj, group.by=str_c("SCT_snn_res.", res), features = rev(unique_markers), cols = c("lightgrey", "dodgerblue3")) + 
    RotatedAxis()+
    coord_flip () 
  print(p2)
  dev.off()
}

save(HarmonyIntegration_obj, file = paste0(out_dir, "rdata/HarmonyIntegration_obj_40PCs_different_resolutions.rdata"))

#### subset OE celltypes and reintegrate ####
load(paste0(out_dir, "rdata/HarmonyIntegration_obj_40PCs_different_resolutions.rdata"))
Idents(HarmonyIntegration_obj) <- "SCT_snn_res.0.1"
OE_HarmonyIntegration_obj <- subset(HarmonyIntegration_obj, idents = c(0, 1, 2, 3, 4, 6, 7, 10, 14, 16, 17))

DefaultAssay(OE_HarmonyIntegration_obj) <- "RNA"

OE_HarmonyIntegration_obj <- SCTransform(OE_HarmonyIntegration_obj, vars.to.regress = c("percent.mt")) %>% 
  RunPCA(verbose = FALSE)
save(OE_HarmonyIntegration_obj, file = "integrated_analysis_seuratV5/rdata/OE_HarmonyIntegration_obj_sctransform.rdata")

# HarmonyIntegration
OE_Harmony_reintegration_obj <- IntegrateLayers(object = OE_HarmonyIntegration_obj, method = HarmonyIntegration, 
                                                assay = "SCT", normalization.method = "SCT",
                                                orig.reduction = "pca",
                                                new.reduction = "integrated.Harmony",
                                                verbose = FALSE)
OE_Harmony_reintegration_obj <- RunUMAP(OE_Harmony_reintegration_obj,  reduction = "integrated.Harmony", reduction.name = "umap.Harmony", dims=1:30)
DimPlot(OE_Harmony_reintegration_obj, reduction = "umap.Harmony", group.by = "orig.ident") + ggtitle("Harmony Integration")
save(OE_HarmonyIntegration_obj, file = "integrated_analysis_seuratV5/rdata/OE_HarmonyIntegration_obj_reintegrate.rdata")

##### tune the PC number ####
HarmonyIntegration_emb <- Embeddings(OE_Harmony_reintegration_obj, "integrated.Harmony")  
# 计算每列（每个维度）的标准差
HarmonyIntegration_stdevs <- apply(HarmonyIntegration_emb, 2, sd)

# 画 Elbow Plot
plot(x = 1:length(HarmonyIntegration_stdevs), y = HarmonyIntegration_stdevs, type = "b", pch = 20,
     xlab = "Dimensions", ylab = "Standard Deviation",
     main = "Elbow Plot for integrated.Harmony")

OE_HarmonyIntegration_obj_marker <- c(
  "Nqo1","Ncam2", # dorsal/ventral OSNs
  "Omp", "Adcy3", "Gnal", "Cnga4", "Cnga2",  "Gng13", "Stoml3", "Rtp1", # Mature OSNs
  "Gap43", "Gnas", "Gng8", "Ablim1", "Dbn1", "Dpysl3", "Dpysl5", "Crmp1", "Marcksl1",  "Stmn1", "Stmn2", "Ebf2", "Atf5",  "Hdac2", "Cbx8",  #Immature OSNs
  "Neurog1", "Neurod1","Neurod2", "Tex15", "Top2a", "Mki67", #INPs, immediate neuronal precursors 
  "Ascl1","Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4", "Top2a",#GBCs
  "Krt5","Krt14", "Trp63", #HBCs
  "Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3", # Ms4a-expressing chemosensory receptor cell (Ms4)
  "Ermn","Cyp2g1","Cyp1a2", "Sult1c1", "Gpx6", "Sox2", "Notch2", "Notch3", "Hes1", "Hes5",  "Hey1",  # Sustentacular cells
  "Coch", "Ascl3", "Cftr", "Foxi1","Hepacam2", # olfactory microvillar cells (Trpm5-)
  "Sox9", "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Ltc4s" # mv (Trpm5+)
)

OE_HarmonyIntegration_obj_marker <- OE_HarmonyIntegration_obj_marker[OE_HarmonyIntegration_obj_marker %in% rownames(OE_Harmony_reintegration_obj)]
all(OE_HarmonyIntegration_obj_marker %in% rownames(OE_Harmony_reintegration_obj))
unique_OE_HarmonyIntegration_obj_marker <- unique(OE_HarmonyIntegration_obj_marker)
length(unique_OE_HarmonyIntegration_obj_marker)
# 70

for (nPCs in seq(20,40,5)){
  OE_Harmony_reintegration_obj <- FindNeighbors(OE_Harmony_reintegration_obj, dims = 1:nPCs, reduction = "integrated.Harmony")
  OE_Harmony_reintegration_obj <- FindClusters(OE_Harmony_reintegration_obj, resolution = 0.2, reduction = "integrated.Harmony")
  OE_Harmony_reintegration_obj <- RunUMAP(OE_Harmony_reintegration_obj, dims = 1:nPCs, reduction = "integrated.Harmony", reduction.name = "umap")
  
  pdf(str_c(out_dir,"plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(OE_Harmony_reintegration_obj, reduction = "umap", label=TRUE, raster=FALSE) + 
    theme(legend.position = "bottom")
  p2 <- DimPlot(OE_Harmony_reintegration_obj, reduction = "umap", group.by="orig.ident", raster=FALSE) + 
    theme(legend.position = "bottom")
  print(p1)
  print(p2)
  dev.off()
  
  DefaultAssay(OE_Harmony_reintegration_obj) <- "SCT"
  pdf(str_c(out_dir,"plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(OE_Harmony_reintegration_obj, reduction = "umap", features=unique_OE_HarmonyIntegration_obj_marker, combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)
  print(p)
  dev.off()
}
FeaturePlot(OE_Harmony_reintegration_obj, reduction = "umap", features="Cd36", combine=FALSE, cols=c("lightgrey", "red"), order=TRUE, raster=FALSE)

##### tune the resolution (40 PCs) #### 
OE_Harmony_reintegration_obj <- FindNeighbors(OE_Harmony_reintegration_obj, reduction = "integrated.Harmony", dims = 1:40)
OE_Harmony_reintegration_obj <- RunUMAP(OE_Harmony_reintegration_obj, dims = 1:40, reduction = "integrated.Harmony", reduction.name = "umap.Harmony")

OE_Harmony_reintegration_obj_mat <- GetAssayData(OE_Harmony_reintegration_obj, assay = "SCT", layer = "data")  # 通常用 "data" slot

pdf(str_c(out_dir, "plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC40_marker_UMAP.pdf"), width = 6, height = 6)
for (name in unique_OE_HarmonyIntegration_obj_marker) {
  print(FeaturePlot(OE_Harmony_reintegration_obj, features = name, raster=FALSE) +
          ggtitle(name))  # 推荐加 print()，确保图输出到 PDF
}
dev.off()

marker_sets <- list(
  mOSN = c("Omp", "Adcy3", "Gnal", "Cnga4", "Cnga2",  "Gng13", "Stoml3", "Rtp1"), # Mature OSNs
  imOSN = c("Gap43", "Gnas", "Gng8", "Ablim1", "Dbn1", "Dpysl3", "Dpysl5", "Crmp1", "Marcksl1",  "Stmn1", "Stmn2", "Ebf2", "Atf5"),  #Immature OSNs
  INP = c("Neurog1", "Neurod1","Neurod2", "Tex15"), #INPs, immediate neuronal precursors 
  GBC = c("Ascl1", "Kit", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Ezh2", "Tmprss4"),#GBCs
  HBC = c("Krt5","Krt14", "Trp63"), #HBCs
  Ms4a = c("Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3"), # Ms4a-expressing chemosensory receptor cell (Ms4)
  SUS = c("Ermn","Cyp2g1","Cyp1a2", "Sult1c1", "Gpx6", "Sox2", "Notch2", "Notch3", "Hes1", "Hes5",  "Hey1"),  # Sustentacular cells
  MV1 = c("Coch", "Ascl3", "Cftr", "Foxi1","Hepacam2"), # olfactory microvillar cells (Trpm5-)
  MV2 = c("Sox9", "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Ltc4s")
)

for (name in names(marker_sets)) {
  genes <- intersect(marker_sets[[name]], rownames(OE_Harmony_reintegration_obj_mat))  # 确保存在
  OE_Harmony_reintegration_obj[[paste0(name, "_marker_sum")]] <- Matrix::colSums(OE_Harmony_reintegration_obj_mat[genes, , drop = FALSE])
}

pdf(str_c(out_dir, "plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC40_marker_sum_UMAP.pdf"), width = 6, height = 6)
for (name in names(marker_sets)) {
  plot_feature <- paste0(name, "_marker_sum")  # 构建 meta.data 中的列名
  print(FeaturePlot(OE_Harmony_reintegration_obj, features = plot_feature,raster=FALSE) +
          ggtitle(plot_feature))  # 推荐加 print()，确保图输出到 PDF
}
dev.off()

resolutions <- seq(0.1, 1.5, by = 0.1)
for (res in resolutions) {
  OE_Harmony_reintegration_obj <- FindClusters(OE_Harmony_reintegration_obj, resolution = res)
  
  pdf(str_c(out_dir,"plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC40_resolution",res,"_UMAP.pdf"))
  p1 <- DimPlot(OE_Harmony_reintegration_obj, reduction = "umap.Harmony", label=TRUE, group.by=str_c("SCT_snn_res.", res))
  p2 <- DimPlot(OE_Harmony_reintegration_obj, reduction = "umap.Harmony", group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  
  pdf(width = 15, height = 80, str_c(out_dir,"plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC40_resolution", res,"_vlnplot_marker.pdf"))
  p1 <- VlnPlot(OE_Harmony_reintegration_obj, group.by=str_c("SCT_snn_res.", res), features = unique_OE_HarmonyIntegration_obj_marker, pt.size = 0) & xlab(NULL) & ylab(NULL)
  print(p1)
  dev.off()
  
  pdf(width = 24, height = 20, str_c(out_dir,"plots/integrated_OE_HarmonyIntegration_obj_scRNA_PC40_resolution",res,"_dotplot_marker.pdf"))
  p2 <- DotPlot(OE_Harmony_reintegration_obj, group.by=str_c("SCT_snn_res.", res), features = rev(unique_OE_HarmonyIntegration_obj_marker), cols = c("lightgrey", "dodgerblue3")) + 
    RotatedAxis()+
    coord_flip () 
  print(p2)
  dev.off()
}

save(OE_Harmony_reintegration_obj, file = paste0(out_dir, "rdata/OE_Harmony_reintegration_obj_different_resolutions.rdata"))

##### Findmarker (resolution = 0.6) #####
load("integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_different_resolutions.rdata")
Idents(OE_Harmony_reintegration_obj) <- "SCT_snn_res.0.6"
OE_Harmony_reintegration_obj <- PrepSCTFindMarkers(OE_Harmony_reintegration_obj)

clusters <- levels(Idents(OE_Harmony_reintegration_obj))

markers_list <- map(clusters, ~ FindMarkers(
  OE_Harmony_reintegration_obj,
  ident.1 = .x,
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
))

# # DE analysis RNA assay data layer
DefaultAssay(OE_Harmony_reintegration_obj) <- "RNA"
OE_Harmony_reintegration_obj <- JoinLayers(OE_Harmony_reintegration_obj, assay = "RNA")
OE_Harmony_reintegration_obj <- NormalizeData(
  OE_Harmony_reintegration_obj,
  assay = "RNA",
  layer = "counts"
)
save(OE_Harmony_reintegration_obj, file = "integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_with_RNA_assay.rdata")

markers_list <- purrr::map(clusters, ~ FindMarkers(
  OE_Harmony_reintegration_obj,
  ident.1 = .x,
  assay = "RNA",
  layer = "data",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox",
  only.pos = TRUE,
  max.cells.per.ident = 5000
))

names(markers_list) <- clusters
save(markers_list, file = "integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_res0.6_cluster_markers_list.rdata")
load("integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_res0.6_cluster_markers_list.rdata")
markers_df <- bind_rows(
  lapply(clusters, function(cl) {
    df <- markers_list[[cl]]
    df$cluster <- cl
    df$gene <- rownames(df)
    df
  }),
  .id = "cluster_id"
)

write.csv(markers_df, file = "integrated_analysis_seuratV5/files/OE_Harmony_reintegration_obj_res0.6_clusters_markers.csv")

table(OE_Harmony_reintegration_obj$orig.ident, OE_Harmony_reintegration_obj$SCT_snn_res.0.6)
# 0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22   23   24   25   26   27    3    4    5    6    7    8    9
# GSE120199_P0_rep1     107  124   13  107    7  109  102   25   60   78    0    0  212   10    3    8   26   16   27    2    1   57  366   40  255  299  110   61
# GSE120199_P0_rep2      99  119   16  100    7  109  105   24   60   79    0    0  217   13    2    6   26   16   28    2    0   58  380   46  253  298  106   60
# GSE120199_P21_rep1    195  237   35   32   26   77   22    7   10   17    0    4  176   26    3    4   26    6    0    0    0  141   87   45   36   42   70   23
# GSE120199_P21_rep2    198  237   32   33   24   79   21    8   12   17    0    4  182   26    3    4   26    6    0    0    0  141   85   44   37   42   75   19
# GSE120199_P3_rep1      44  103   12   93    3  264  254   37  156  192    0    0  184   22    2    7   26   18   67    0    7  112  281   29  204  273  117   51
# GSE120199_P3_rep2      38  108    8   93    8  255  255   40  152  191    0    0  178   19    2    8   28   17   66    1    7  114  292   29  205  277  118   52
# GSE120199_P7_rep1      55  132    8   69    7   86   42   42   37   38    0    1  117    7    2   12   22    3   24    0    0  103  228   20  154  155   66   33
# GSE120199_P7_rep2      55  131    9   68    8   85   42   42   39   37    0    1  110    7    2   11   21    4   25    0    0  105  227   23  155  152   70   32
# GSE147459_adult       591  454  104   38   69   30   26   12   16   24    0    5  379   81   12   17    7    1    0    0    0  111  102   73   49   47  131   28
# GSE147459_E18.5        19   17    2   69    0  501  172  499   99  112    0    0   68    9    1   13   54   25   68    2    1   86  229   16  216  192   86   26
# GSE147459_P14         493  354   44   91   47   42   33   28   38   41    1    2  307    9   10   21   28   11    4    5    0  111  372   49  113   94  108   44
# GSE173947_homecage_1 2397 2246  498  520  320   81    6   14   19   13  290   98 1389  109   83    9   36   11    0   94    0  544  631  591  594  553  497  700
# GSE173947_homecage_2 2209 1895  338  349  338   56    5   21    9   18  194   96  745   56   59    5   32    5    0   35    0  500  479  473  471  408  315  504
# GSE173947_homecage_3 2197 1879  319  327  307   50    5    7    5    6  198  131 1049   74   41    4   26   19    0   43    0  359  464  391  411  416  388  514
# GSE173947_homecage_4 2196 1832  403  238  313   52    4   11    4    0  166   81 1247   42   63    3   17   17    0   43    0  514  361  476  284  240  216  389
# GSE173947_homecage_5 2292 2390  480  311  312  109    7   15    8    6  192   81 2049   71   74    4   42   16    0   80    0  600  448  515  419  347  316  469
# GSE173947_homecage_6 2171 2246  435  345  336   92    2   11    1    3  182   92 2015   67   60    4   34   19    0   29    0  572  438  518  399  315  290  418
# GSE224604_19d_rep1   1172 1349  398  123  160   20   10  239   10    5   22   42 1721   18   39   10   25   33    0    9    0  450  210  340  124  104  388  202
# GSE224604_19d_rep2    442  844  223  220   81  298  859  464  686  503   55   15  538   34   45   62   79  116   73   21   13  585  477  116  208  217  274  164
# GSE224604_8w_rep1    1115  696  155    9  145  118    5  113    8    0    0  107  745    8   66  178    0   43    0    0  272  424    9  241    0    8  115   70
# GSE224604_8w_rep2    2721 2132  371   31  411   32    2  173    0    0    0  144 1992   27  147  321    2   31    0    0    0 1077   39  650    6   11  347  147

OE_Harmony_reintegration_obj$celltype <- case_when(
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(0, 1, 2, 3, 5, 10, 12, 19, 21, 25, 26) ~ "mOSN",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(4, 8, 9, 11) ~ "iOSN",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(6, 7, 17) ~ "INP",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(14, 16, 18, 22, 27) ~ "GBC",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(13) ~ "HBC",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(15, 24) ~ "SUS",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(20) ~ "Trpm5_MV",
  OE_Harmony_reintegration_obj$SCT_snn_res.0.6 %in% c(23) ~ "MV",
  TRUE ~ NA_character_
)
OE_Harmony_reintegration_obj$celltype <- factor(OE_Harmony_reintegration_obj$celltype, levels = c("mOSN", "iOSN", "INP", "GBC", "HBC", "SUS", "Trpm5_MV", "MV"))
save(OE_Harmony_reintegration_obj, file = "integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_with_celltype.rdata")

#### receptor expression #### 
# identiﬁcation of the OR expressed in each OSN
gene_info <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T)
receptor_info <- gene_info[grepl("Taar|Olfr", gene_info$gene_name) & grepl("protein_coding", gene_info$gene_type),]
receptor_name <- receptor_info$gene_name
length(receptor_name)
# 1155

##### tune the UMI threshold of receptor expression in mOSN ####
load("integrated_analysis_seuratV5/rdata/OE_Harmony_reintegration_obj_with_celltype.rdata")
Idents(OE_Harmony_reintegration_obj) <- "celltype"
DimPlot(OE_Harmony_reintegration_obj, group.by = "celltype")
mOSN_diet_scRNA <- DietSeurat(
  object = subset(OE_Harmony_reintegration_obj, idents = "mOSN"),
  assays = "SCT",
  layers = c("counts", "data")   
)
# DietSeurat
# Save memory
# Reduce file size for sharing/saving
# Speed up processing
# Focus only on specific components of the data

mOSN_olfr_counts <- FetchData(
  object = mOSN_diet_scRNA,
  vars = receptor_name,
  layer = "counts"
)
# Warning message:
#   The following requested variables were not found (10 out of 65 shown): Olfr248, Olfr340, Olfr345, Olfr1008, Olfr1022, Olfr1040, Olfr1053, Olfr1061, Olfr1066, Olfr1128 
dim(mOSN_olfr_counts)
# [1] 56536  1087
# 76643  1087

calculate_mOSNs_percentage <- function(cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(mOSN_olfr_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(mOSN_olfr_counts>=cutoff)
  }
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  mOSN_df <- data.frame(cutoff=cutoff,Olfr_N=c("0","1",">1"),percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(mOSN_df)
}

cutoffs <- c(0:8,16,32,64,128,256)
mOSN_df <- c()
for (i in 1:length(cutoffs)){
  mOSN_df <- calculate_mOSNs_percentage(cutoffs[i]) %>% rbind(mOSN_df, .)
}
mOSN_df$cutoff <- factor(mOSN_df$cutoff,levels=cutoffs)
mOSN_df$Olfr_N <- factor(mOSN_df$Olfr_N,levels=c("0","1",">1"))
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])

pdf(str_c(out_dir,"plots/mOSNs_UMIs_cutoff_threshold_tunning.pdf"), height = 8, width=9)
ggplot(data=mOSN_df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
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
        axis.line = element_line(colour="black",linewidth=1))
dev.off()

##### mOSN (UMI ≥ 3) #####
# 设置阈值，小于 3 的置零
mOSN_olfr_counts[mOSN_olfr_counts < 3] <- 0

# 选出单受体细胞（只表达一个受体）
single_cells <- rownames(mOSN_olfr_counts)[rowSums(mOSN_olfr_counts > 0) == 1]

mOSN_scRNA <- subset(OE_Harmony_reintegration_obj, idents = "mOSN")

# subset 得到单受体细胞对象
mOSN_scRNA_single_receptor <- subset(mOSN_scRNA, cells = single_cells)

# 给每个 cell 标记 receptor
mOSN_scRNA_single_receptor$receptor <- apply(
  mOSN_olfr_counts[single_cells, ], 1, function(x) {
    colnames(mOSN_olfr_counts)[which(x > 0)]
  }
)

###### add class info (functional receptor) #####
gene_info <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T)
TAAR_info <- gene_info[grepl("Taar", gene_info$gene_name) & grepl("protein_coding", gene_info$gene_type),]
Taar_name <- TAAR_info$gene_name

ClassI_ClassII_OR_info <- read.csv("/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
ClassI_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_I" & ClassI_ClassII_OR_info$gene_type == "protein_coding")]%>% 
  gsub("\\.ps", "-ps",  .)
length(ClassI_OR_name)
# 131

ClassII_OR_name <- rownames(ClassI_ClassII_OR_info)[which(ClassI_ClassII_OR_info$OR_Class == "Class_II" & ClassI_ClassII_OR_info$gene_type == "protein_coding")] %>% 
  gsub("\\.ps", "-ps",  .)
length(ClassII_OR_name)
# 1006

mOSN_scRNA_single_receptor$receptor_class <- ifelse(mOSN_scRNA_single_receptor$receptor %in% Taar_name ,"Taar",
                                                    ifelse(mOSN_scRNA_single_receptor$receptor %in% ClassI_OR_name, "class_I_OR",
                                                           ifelse(mOSN_scRNA_single_receptor$receptor %in% ClassII_OR_name,"class_II_OR",
                                                                  "undefined")))

Idents(mOSN_scRNA_single_receptor) <- "receptor_class"
table(mOSN_scRNA_single_receptor$receptor_class)
# class_I_OR class_II_OR        Taar 
# 7142       51607         724 
save(mOSN_scRNA_single_receptor,file = "integrated_analysis_seuratV5/rdata/mOSN_expressing_single_receptor.rdata")

##### tune the UMI threshold of receptor expression in iOSN ####
Idents(OE_Harmony_reintegration_obj) <- "celltype"
DefaultAssay(OE_Harmony_reintegration_obj) <- "SCT"
iOSN_diet_scRNA <- DietSeurat(
  object=subset(OE_Harmony_reintegration_obj, idents = "iOSN"),
  layers = c("counts", "data"),
  assays="SCT"
)
# DietSeurat
# Save memory
# Reduce file size for sharing/saving
# Speed up processing
# Focus only on specific components of the data

iOSN_olfr_counts <- FetchData(
  object = iOSN_diet_scRNA,
  vars = receptor_name,
  layer = "counts"
)
# Warning message:
#   The following requested variables were not found (10 out of 65 shown): Olfr248, Olfr340, Olfr345, Olfr1008, Olfr1022, Olfr1040, Olfr1053, Olfr1061, Olfr1066, Olfr1128 
dim(iOSN_olfr_counts)
# 17680  1087

calculate_iOSNs_percentage <- function(cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(iOSN_olfr_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(iOSN_olfr_counts>=cutoff)
  }
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  iOSN_df <- data.frame(cutoff=cutoff,Olfr_N=c("0","1",">1"),percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(iOSN_df)
}

cutoffs <- c(0:8,16,32,64,128,256)
iOSN_df <- c()
for (i in 1:length(cutoffs)){
  iOSN_df <- calculate_iOSNs_percentage(cutoffs[i]) %>% rbind(iOSN_df, .)
}
iOSN_df$cutoff <- factor(iOSN_df$cutoff,levels=cutoffs)
iOSN_df$Olfr_N <- factor(iOSN_df$Olfr_N,levels=c("0","1",">1"))
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])

pdf(str_c(out_dir,"plots/iOSNs_UMIs_cutoff_threshold_tunning.pdf"), height = 8, width=9)
ggplot(data=iOSN_df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
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
        axis.line = element_line(colour="black",linewidth=1))
dev.off()

##### iOSN without threshold #####
iOSN_scRNA <- subset(OE_Harmony_reintegration_obj, idents = "iOSN")
iOSN_scRNA$receptor_number <- rowSums(iOSN_olfr_counts>0)
table(iOSN_scRNA$receptor_number)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   28   29   30 
# 3296 4017 3853 2737 1587  839  455  272  172  131   87   63   42   44   29   18   11    8    3    5    2    3    1    1    1    1    1    1 
table(iOSN_scRNA$receptor_number,iOSN_scRNA$orig.ident)
FeaturePlot(iOSN_scRNA, features = "receptor_number",order = T)
DimPlot(OE_Harmony_reintegration_obj, reduction = "umap.Harmony")

Idents(iOSN_scRNA) <- "orig.ident"
iOSN_scRNA_filter <- subset(iOSN_scRNA, idents = unique(grep("homecage|GSE224604", iOSN_scRNA$orig.ident, value = T)))
table(iOSN_scRNA_filter$receptor_number)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   28   29   30 
# 322 2602 3479 2619 1552  833  452  270  172  131   86   63   42   44   29   18   11    8    3    5    2    3    1    1    1    1    1    1 
unique(iOSN_scRNA$orig.ident)
# [1] "GSE120199_P0_rep1"    "GSE120199_P0_rep2"    "GSE120199_P3_rep1"    "GSE120199_P3_rep2"    "GSE120199_P7_rep1"    "GSE120199_P7_rep2"    "GSE120199_P21_rep1"   "GSE120199_P21_rep2"   "GSE147459_E18.5"     
# [10] "GSE147459_P14"        "GSE147459_adult"      "GSE173947_homecage_1" "GSE173947_homecage_2" "GSE173947_homecage_3" "GSE173947_homecage_4" "GSE173947_homecage_5" "GSE173947_homecage_6" "GSE224604_19d_rep1"  
# [19] "GSE224604_19d_rep2"   "GSE224604_8w_rep1"    "GSE224604_8w_rep2" 

iOSN_scRNA_filter$receptor_number_class <- factor(
  ifelse(iOSN_scRNA_filter$receptor_number == 0, "0",
         ifelse(iOSN_scRNA_filter$receptor_number == 1, "1", "multiple")),
  levels = c("0", "1", "multiple")
)

table(iOSN_scRNA_filter$receptor_number_class)
# 0        1 multiple 
# 322     2602     9828 

#### plot: receptor coexpression in iOSN (expressing multiple receptors) ####
Idents(iOSN_scRNA_filter) <- "receptor_number_class"
iOSN_scRNA_filter_multiple_receptors <- subset(iOSN_scRNA_filter, idents = "multiple")

common_ClassI_OR_genes <- intersect(ClassI_OR_name, colnames(iOSN_olfr_counts))
iOSN_scRNA_filter_multiple_receptors$class_I_OR <- rowSums(
  iOSN_olfr_counts[, common_ClassI_OR_genes, drop = FALSE]
) > 0

common_ClassII_OR_genes <- intersect(ClassII_OR_name, colnames(iOSN_olfr_counts))
iOSN_scRNA_filter_multiple_receptors$class_II_OR <- rowSums(
  iOSN_olfr_counts[, common_ClassII_OR_genes, drop = FALSE]
) > 0

common_Taar_genes <- intersect(Taar_name, colnames(iOSN_olfr_counts))
iOSN_scRNA_filter_multiple_receptors$Taar <- rowSums(
  iOSN_olfr_counts[, common_Taar_genes, drop = FALSE]
) > 0

# 三列逻辑值
coexp_df <- iOSN_scRNA_filter_multiple_receptors@meta.data[, c("class_I_OR", "class_II_OR", "Taar")]

# 统计所有组合的细胞数
coexp_table <- as.data.frame(table(
  Class_I_OR  = coexp_df$class_I_OR,
  Class_II_OR = coexp_df$class_II_OR,
  Taar        = coexp_df$Taar
))

coexp_table
write.csv(coexp_table, file = "integrated_analysis_seuratV5/files/recepotr_class_coexpression_in_iOSNs_expressing_multiple_receptors.csv")

pdf("integrated_analysis_seuratV5/plots/iOSN_multiple_receptors_coexpression_cell_number.pdf", width = 10, height = 6)
upset(
  iOSN_scRNA_filter_multiple_receptors@meta.data,
  c("class_I_OR", "class_II_OR", "Taar"),
  name = "Cells per combination"
) 
dev.off()

