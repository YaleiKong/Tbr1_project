Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list=ls())
getwd()  

## set seed for repeat
set.seed(42)

## library 
library_Ls <- c("DESeq2","LSD","gProfileR",
                "stringr","readxl",
                "dplyr","tidyr","tidyverse",
                "ggplot2","pheatmap","ggrepel","gg.gap","RColorBrewer","grid",
                "gridExtra","gplots","reshape2","ggpubr","ggbreak","gtools","textshape",
                "circlize","ComplexHeatmap","hutils"
)
lapply(library_Ls,function(x){suppressPackageStartupMessages(library(x,character.only = T))})

## set the workpath
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/6_lsl_Tbr1_Pax6_Cre_bulk_RNA_Seq/20240516_lsl_Tbr1_Pax6_Cre_bulk_RNA_Seq_202405/DE_analysis_20240516")
setwd(workpath)

## create a result folder
folder <- "Results/"
if (!file.exists(folder)) { dir.create(folder) }

folder_plot <- "Plots/"
if (!file.exists(folder_plot)) { dir.create(folder_plot) }

#### data preparation ####
## Load data
dat <- read.table("../featureCounts/featureCounts_multiOverlap_geneId.txt", header=T, row.names=1, as.is=T)
# as.is=T means keep all character vectors, do not convert to factors
dat_Cre_tdTomato <- read.table("../featureCounts/Cre_tdTomato_featureCounts.txt", header=T, row.names=1, as.is=T) %>% 
  t()

## load ensembl gene name and gene description file
# import ensembl info and remove the duplicated gene_name
# unique(ensembl_gene_list_good[grep("Olfr", ensembl_gene_list_good$gene_name),"gene_type"])
# [1] "protein_coding"                     "unprocessed_pseudogene"             "transcribed_unprocessed_pseudogene" "polymorphic_pseudogene" 
ensembl_gene_list_good <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) %>% 
  mutate(gene_id_Ensembl_remove_version_meta = str_remove(gene_id_Ensembl, "\\.[0-9]")) %>% 
  filter(gene_type %in%  c("protein_coding", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "polymorphic_pseudogene")) %>% 
  filter(grepl("^chr[1-9XY]", chr_name)) %>% 
  filter(!duplicated(gene_name))
dim(ensembl_gene_list_good)
# [1] 24899    15

## match the ensembl ID of the data to the ensembl ID from ensembl_gene_list_good
good_dat <- dat[rownames(dat) %in% ensembl_gene_list_good$gene_id_Ensembl, ]

## obtain the gene names
good_dat_matched_gene_name_with_gene_length <- good_dat %>%
  rownames_to_column(var = "gene_id") %>% 
  mutate(gene_name = ensembl_gene_list_good[match(gene_id, ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name) %>% 
  filter(!duplicated(gene_name)) %>%
  dplyr::rename(gene_length = Length) 
good_dat_matched_gene_name_with_gene_length <- column_to_rownames(good_dat_matched_gene_name_with_gene_length, "gene_name")

###### sample info #####
meta <- data.frame(
  sample = paste0("sample_", 1:18),
  sample_name = c(paste0("lsl_Tbr1_tdTomato_het_Pax6_Cre_het_", 1:9), paste0("lsl_Tbr1_tdTomato_het_", 1:9)),
  genotype = rep(c("lsl_Tbr1_tdTomato_het_Pax6_Cre_het", "lsl_Tbr1_tdTomato_het"), each = 9),
  mice_info = c(rep("group1", 2), rep("group2",7),
                rep("group1", 6), rep("group2",2), "group3")
)

## change the colnames of data
# 提取原始列名
original_colnames <- colnames(good_dat_matched_gene_name_with_gene_length)

# 创建一个映射表，将 sample 列名映射到 sample_name
colname_mapping <- setNames(meta$sample_name, meta$sample)

# 使用 ifelse 和 %in% 语句替换列名
new_colnames <- ifelse(original_colnames %in% meta$sample, colname_mapping[original_colnames], original_colnames)

# 赋值新的列名给数据框
colnames(good_dat_matched_gene_name_with_gene_length) <- new_colnames

## export the clean data
write.csv(good_dat_matched_gene_name_with_gene_length,"Results/good_dat_matched_gene_name_with_gene_length.csv")

good_dat_matched_gene_name <- good_dat_matched_gene_name_with_gene_length[,-grep("gene", colnames(good_dat_matched_gene_name_with_gene_length))] %>% 
  dplyr::select(10:18, 1:9)

write.csv(good_dat_matched_gene_name,"Results/good_dat_matched_gene_name.csv")

#### calculate TPM ####
good_dat_matched_gene_name_with_gene_length <- read.csv("Results/good_dat_matched_gene_name_with_gene_length.csv", header = T, row.names = 1)
gene_length_kb <- good_dat_matched_gene_name_with_gene_length$gene_length / 1000
rpk_raw_count <- good_dat_matched_gene_name_with_gene_length[, -grep("gene", colnames(good_dat_matched_gene_name_with_gene_length))] / gene_length_kb 
tpm_raw_count <- t(t(rpk_raw_count)/colSums(rpk_raw_count) * 1000000) %>% as.data.frame()
colnames(tpm_raw_count) <- paste0(colnames(good_dat_matched_gene_name_with_gene_length)[-grep("gene", colnames(good_dat_matched_gene_name_with_gene_length))], "_TPM")
write.csv(tpm_raw_count,file="Results/tpm.csv")

#### DESeq2 ####
good_dat_matched_gene_name <- read.csv("Results/good_dat_matched_gene_name.csv", header = T, row.names = 1)

###load metadata which is the sample categories
meta <- meta[,c("sample_name", "genotype", "mice_info")]
row.names(meta) <- meta$sample_name
dim(meta)
summary(meta)

### Check that sample names match in both files (dat and meta)
# "all" function to check whether all the logic are true
which(row.names(meta) %in% colnames(good_dat_matched_gene_name))  # return the position of the former
match(row.names(meta) , colnames(good_dat_matched_gene_name))     # return the position of the latter
all(row.names(meta) == colnames(good_dat_matched_gene_name))

### create a DESeqDataSet object
good_dat_matched_gene_name[] <- lapply(good_dat_matched_gene_name, function(x) as.integer(round(x)))
dds_dat <- DESeqDataSetFromMatrix(countData = good_dat_matched_gene_name, colData = meta, design= ~ mice_info + genotype)
dds_dat$genotype <- relevel(dds_dat$genotype, "lsl_Tbr1_tdTomato_het")  

##### PCA #####
vsd <- vst(dds_dat, blind = FALSE)
DESeq2::plotPCA(vsd, intgroup = "genotype")
pca_data <- DESeq2::plotPCA(vsd, intgroup = "genotype", returnData = TRUE)

# ggplot2
percentVar <- round(100 * attr(pca_data, "percentVar"))
p.pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(lsl_Tbr1_tdTomato_het_Pax6_Cre_het = "orangered",
                                lsl_Tbr1_tdTomato_het = "springgreen"))+
  geom_text_repel(aes(label = pca_data$name), color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 15),
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title=element_blank(),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm'),
        panel.grid = element_blank())

ggsave("Plots/PCA.pdf", p.pca, width = 15, height = 10)

#### DESeq2 (remove lsl_Tbr1_tdTomato_het_9) ####
good_dat_matched_gene_name <- read.csv("Results/good_dat_matched_gene_name.csv", header = T, row.names = 1)
good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9 <- good_dat_matched_gene_name %>% 
  dplyr::select(-"lsl_Tbr1_tdTomato_het_Pax6_Cre_het_9")

### load metadata which is the sample categories
meta <- meta[,c("sample_name", "genotype","mice_info")]
meta_exclude_lsl_Tbr1_tdTomato_het_9 <- meta %>% 
  dplyr::filter(sample_name!="lsl_Tbr1_tdTomato_het_Pax6_Cre_het_9")

row.names(meta_exclude_lsl_Tbr1_tdTomato_het_9) <- meta_exclude_lsl_Tbr1_tdTomato_het_9$sample_name
dim(meta_exclude_lsl_Tbr1_tdTomato_het_9)
summary(meta_exclude_lsl_Tbr1_tdTomato_het_9)

### Check that sample names match in both files (dat and meta_exclude_lsl_Tbr1_tdTomato_het_9)
# "all" function to check whether all the logic are true
which(row.names(meta_exclude_lsl_Tbr1_tdTomato_het_9) %in% colnames(good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9))
match(row.names(meta_exclude_lsl_Tbr1_tdTomato_het_9) , colnames(good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9))
all(row.names(meta_exclude_lsl_Tbr1_tdTomato_het_9) == colnames(good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9))

### create a DESeqDataSet object
good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9[] <- lapply(good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9, function(x) as.integer(round(x)))
dds_dat_exclude_lsl_Tbr1_tdTomato_het_9 <- DESeqDataSetFromMatrix(countData = good_dat_matched_gene_name_exclude_lsl_Tbr1_tdTomato_het_9, colData = meta_exclude_lsl_Tbr1_tdTomato_het_9, design = ~ mice_info + genotype)
# design = ~Strain + Time means that deseq2 will test the effect of the Time (the last factor), controlling for the effect of Strain (the first factor), so that the algorithm returns the fold change result only from the effect of time. 
dds_dat_exclude_lsl_Tbr1_tdTomato_het_9$genotype <- relevel(dds_dat_exclude_lsl_Tbr1_tdTomato_het_9$genotype, "lsl_Tbr1_tdTomato_het")  

##### DEseq2 #####
dds_DE_exclude_lsl_Tbr1_tdTomato_het_9 <- DESeq(dds_dat_exclude_lsl_Tbr1_tdTomato_het_9)

##### export results #####
## export the DESeq2 results
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- results(dds_DE_exclude_lsl_Tbr1_tdTomato_het_9, contrast = c("genotype", "lsl_Tbr1_tdTomato_het_Pax6_Cre_het", "lsl_Tbr1_tdTomato_het"), test="Wald", independentFiltering = FALSE) %>% 
  as.data.frame() %>% 
  mutate(gene_type = ensembl_gene_list_good$gene_type,
         gene_id_Ensembl = ensembl_gene_list_good$gene_id_Ensembl,
         gene_id_Ensembl_remove_version_info = ensembl_gene_list_good$gene_id_Ensembl_remove_version_meta,
         chr_name = ensembl_gene_list_good$chr_name,
         genomic_strand = ensembl_gene_list_good$genomic_strand,
         genomic_start = ensembl_gene_list_good$genomic_start,
         genomic_end = ensembl_gene_list_good$genomic_end,
         gene_name = rownames(.)) %>% 
  na.omit()
# set **independentFiltering to FALSE** so that there will be less NA values of padj (see below for padj). 
# Especially for lowly expressed genes such as olfactory receptors. 
# The independent filtering is designed only to filter out low count genes to the extent that they are not enriched with small p-values.
write.csv(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, "Results/DEseq2_results_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9.csv")

## export the normalized counts
normalized_count_exclude_lsl_Tbr1_tdTomato_het_9 <- DESeq2::counts(dds_DE_exclude_lsl_Tbr1_tdTomato_het_9, normalized = TRUE) %>% 
  as.data.frame() %>% 
  mutate(gene_type = ensembl_gene_list_good$gene_type,
         gene_id_Ensembl = ensembl_gene_list_good$gene_id_Ensembl,
         gene_id_Ensembl_remove_version_info = ensembl_gene_list_good$gene_id_Ensembl_remove_version_meta,
         chr_name = ensembl_gene_list_good$chr_name,
         genomic_strand = ensembl_gene_list_good$genomic_strand,
         genomic_start = ensembl_gene_list_good$genomic_start,
         genomic_end = ensembl_gene_list_good$genomic_end,
         gene_name = rownames(.)) %>% 
  na.omit()
write.csv(normalized_count_exclude_lsl_Tbr1_tdTomato_het_9, "Results/DEseq2_normalized_count_exclude_lsl_Tbr1_tdTomato_het_9.csv")

##### check the results of olfactory receptor genes #####
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- read.csv("Results/DEseq2_results_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9.csv", header = T, row.names = 1)
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "^Olfr", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                         grepl(pattern = "protein_coding", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ] 
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$padj < 0.05 & 
                                                                                                                                                  abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$log2FoldChange) > 0.585), ]
dim(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff)
# [1] 118  14
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "^Olfr", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                         grepl(pattern = "pseudogene", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$padj < 0.05 & 
                                                                                                                                                  abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$log2FoldChange) > 0.585), ]
dim(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR_diff)
# [1]  2 14
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                           grepl(pattern = "protein_coding", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$padj < 0.05 & 
                                                                                                                                                      abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$log2FoldChange) > 0.585), ]
dim(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff)
# [1]  1 14
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "^Taar", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                           grepl(pattern = "pseudogene", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR$padj < 0.05 & 
                                                                                                                                                      abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR$log2FoldChange) > 0.585), ]
dim(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR_diff)
# [1]  0 14

sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj < 0.05 & abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange) > 0.585, na.rm = T)
# 574
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange > 0.585, na.rm = T)
# 269
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange < -0.585, na.rm = T)
# 305

sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$padj < 0.05 & abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$log2FoldChange) > 0.585, na.rm = T)
# 118
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$log2FoldChange > 0.585, na.rm = T)
# 4
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$log2FoldChange < -0.585, na.rm = T)
# 114

sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$padj < 0.05 & abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$log2FoldChange) > 0.585, na.rm = T)
# 2
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$log2FoldChange > 0.585, na.rm = T)
# 0
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_OR$log2FoldChange < -0.585, na.rm = T)
# 2

sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$padj < 0.05 & abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$log2FoldChange) > 0.585, na.rm = T)
# 1
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$log2FoldChange > 0.585, na.rm = T)
# 1
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$padj < 0.05 & res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$log2FoldChange < -0.585, na.rm = T)
# 0
sum(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR$padj < 0.05 & abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_pseudogene_TAAR$log2FoldChange) > 0.585, na.rm = T) 
# 0

##### check the change of marker genes #####
markers <- c("Krt5","Krt14","Trp63", # HBC
             "Hes6","Sox2","Kit","Ascl1","Neurog1","Neurod1", # GBC
             "Gap43","Gnas","Gng8","Hdac2","Stmn4", # imOSN
             "Omp","Gnal","Gng13","Adcy3","Cnga4","Cnga2","Cngb1","Slc17a6","Stoml3", # mOSN
             "Cyp2g1","Cyp1a2","Hes1", #SUS
             "Sox9", "Sox10", # bowman gland
             "Car2", "Gucy2d", "Cnga3", "Emx1", "Ms4a3", "Ms4a4a", "Ms4a4b", "Ms4a4c", "Ms4a6b", "Ms4a6c", "Ms4a6d", "Ms4a7", "Ms4a8a", "Ms4a10", "Ms4a15", # GC-D/Ms4a
             "Trpm5", "Pou2f3", "Chat", "Avil", "Il25", "Il17rb", "Ltc4s", # Trpm5+ MV
             "Itpr3", "Nt5e", "Ascl3", "Npy", "Foxi1", "Coch", "Cftr", "Cd24a", "Hepacam2",  # Trpm5- MV
             "Trpc2", # Trpc2+ 
             "Gucy1b2" # Trpc3+ typeB
)

res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_markers <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[markers,] 
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_markers_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_markers %>% 
  dplyr::filter(padj < 0.05)
dim(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_markers_diff)
# [1] 20 14

#### plot ####
##### plot 1: volcano plot #####
## check the range of log2FC and padj
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- read.csv("Results/DEseq2_results_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9.csv", header = T, row.names = 1) %>% 
  as.data.frame()
range(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange, na.rm = T) 
# [1] -4.317617  4.071962
range(-log10(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj), na.rm = T)
# 2.375320e-05 5.728133e+01

## pvalue & padj
# https://www.biostars.org/p/415023/
# the p-value in DESeq2 is calculated using the wald test. 
# The null hypothesis of the wald test is that: for each gene, there is no differential expression across two sample groups (e.g., treated vs control). 
# If the p-value is small (e.g., p<0.05), the null hypothesis is rejected, as there is only 5% of chance that the null hypothesis is true. 
# However, when you have many genes being tested, by chance (5%), there is a number of genes that are not significantly expressed, but obtained significant p-values. 
# Therefore, we need to correct this problem caused by multiple testing. DESeq2 adjust the p value from wald test using Benjamini and Hochberg method (BH-adjusted p values), which is presented in the column of padj in the results object.

# check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj == 0), ] 
# 0
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[which(grepl("Olfr", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                               res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type %in% "protein_coding"), ]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR<- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[which(grepl("Taar[2-9]", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & 
                                                                                                                                res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type %in% "protein_coding"), ]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$padj < 0.05  & 
                                                                                                                                                  abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR$log2FoldChange) > 0.585) ,]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR[which(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$padj < 0.05  & 
                                                                                                                                                      abs(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR$log2FoldChange) > 0.585) ,]
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_rm_functional_OR_diff_functional_TAAR_diff <- subset(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, !(gene_name %in% c(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff$gene_name, 
                                                                                                                                                                               res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff$gene_name)))

pdf("./Plots/volcano_plot_DEGs.pdf", height = 8, width = 8)
ggplot(data = res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_rm_functional_OR_diff_functional_TAAR_diff, aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point(aes(colour = padj < 0.05  & abs(log2FoldChange) > 0.585), alpha = 0.75, pch = 16, size = 1) + 
  scale_color_manual(values = c("grey", "black")) +
  # The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  
  geom_point(data = res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff, alpha = 0.75, pch = 16, size = 2, color = "blue") + 
  geom_text_repel(data = res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_OR_diff, aes(label = gene_name), color = "blue", size = 3) + 
  
  geom_point(data = res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff, alpha = 0.75, pch = 16, size = 2, color = "red") + 
  geom_text_repel(data = res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR_diff, aes(label = gene_name), color = "red", size = 3) + 
  
  # geom_point(data = subset(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, gene_name %in% "Gucy1b2"), alpha = 0.75, pch = 16, size = 2, color = "green") + 
  # geom_text(data = subset(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, gene_name %in% "Gucy1b2"), aes(label = gene_name)) + 
  
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  xlab("log2 fold change") + 
  ylab("-log10 padj") + 
  xlim(floor(min(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange)), ceiling(max(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,  ceiling(max(-log10(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$padj)))))+
  ggtitle("res_OverExpression_Pax6")
dev.off()

##### plot 2: Taar normalized expression counts bar plot #####
normalized_count_exclude_lsl_Tbr1_tdTomato_het_9_protein_coding_and_pseudogene <- read.csv("./Results/DEseq2_normalized_count_exclude_lsl_Tbr1_tdTomato_het_9.csv", header = T, row.names = 1)
Taar_normalized_count <- normalized_count_exclude_lsl_Tbr1_tdTomato_het_9_protein_coding_and_pseudogene[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count_exclude_lsl_Tbr1_tdTomato_het_9_protein_coding_and_pseudogene)) & normalized_count_exclude_lsl_Tbr1_tdTomato_het_9_protein_coding_and_pseudogene$gene_type %in% "protein_coding"), grepl("Tbr1", colnames(normalized_count_exclude_lsl_Tbr1_tdTomato_het_9_protein_coding_and_pseudogene))]

Taar_normalized_count_for_plot <- Taar_normalized_count %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count)[grep("Tbr1", colnames(Taar_normalized_count))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_het_Pax6_Cre_het", sampletype), 
                           "lsl_Tbr1_tdTomato_het_Pax6_Cre_het", "lsl_Tbr1_tdTomato_het")) %>%  
  left_join(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR[grepl("Tbr1|gene_name|padj", colnames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_het","lsl_Tbr1_tdTomato_het_Pax6_Cre_het")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))) 


pdf("./Plots/bar_plot_taar_normalized_count.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_for_plot,  sampletype == "lsl_Tbr1_tdTomato_het_Pax6_Cre_het_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_for_plot$normalized_count))*0.97), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_for_plot$normalized_count))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts")

dev.off()

Taar_normalized_count_compared_to_control_for_plot <- as.data.frame(Taar_normalized_count/rowMeans(Taar_normalized_count[, grep("lsl_Tbr1_tdTomato_het_\\d",colnames(Taar_normalized_count))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count)[grep("Tbr1", colnames(Taar_normalized_count))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_het_Pax6_Cre_het", sampletype), 
                           "lsl_Tbr1_tdTomato_het_Pax6_Cre_het", "lsl_Tbr1_tdTomato_het")) %>%  
  left_join(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR[grepl("Tbr1|gene_name|padj", colnames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9_functional_TAAR))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_het","lsl_Tbr1_tdTomato_het_Pax6_Cre_het")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))) 

pdf("./Plots/bar_plot_taar_normalized_count_compared_to_control.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_compared_to_control_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + 
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_compared_to_control_for_plot,  sampletype == "lsl_Tbr1_tdTomato_het_Pax6_Cre_het_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_compared_to_control_for_plot$normalized_count))*0.97), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_compared_to_control_for_plot$normalized_count))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts compared to lsl_Tbr1_tdTomato_het")
dev.off()

##### plot 4: log2FC and genome location (Taar gene cluster) ####
taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "^Taar[2-9]|Slc18b1|Stx7|Moxd1", row.names(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & grepl("protein_coding", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ]

pdf("./Plots/taar_log2FC_with_genome_location.pdf", width = 7, height = 7)
ggplot(data = taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, aes(genomic_start, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = T, color = "blue")+ 
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, (padj < 0.05 & abs(log2FoldChange) > 0.585) == F), size = 2, pch = 1, na.rm = T, color = "blue")+ 
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*0.9, ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.3))+
  scale_x_continuous(expand = c(0, 0), limits = c(23790000, 24310000), breaks = seq(23800000, 24300000, by = 100000) , labels = seq(23800000, 24300000, by = 100000)/1000000)+
  
  geom_text_repel(size = 4, colour = 'blue', nudge_x = -0.25, nudge_y = -0.25)+
  
  theme_classic()+ 
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Position on chromosome 10(MB)", y = "log2 Gene expression changes") +
  
  geom_rect(aes(xmin = genomic_start, xmax = genomic_end, ymin = ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange)), ymax = ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.1), colour = "black", fill = "black")+
  
  geom_text(data = data.frame(label = c("Slc18b1", "Taar", "Stx7", "Moxd1"), 
                              x = c(23815000, 23900000, 24168000, 24260000), 
                              y = rep(ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.2, 4)), 
            mapping = aes(x = x, y = y, label = label))+
  geom_text(data = data.frame(label = gsub("Taar", "", taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name[grepl("Taar",taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name)]), 
                              x = ((taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_start+taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_end)/2)[grepl("Taar",taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name)], 
                              y = rep(ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.2, 14)), 
            mapping = aes(x = x, y = y, label = label), 
            size = 2
  )

dev.off()

pdf("./Plots/taar_log2FC_with_genome_location_change_color.pdf", width = 7, height = 7)
ggplot(data = taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, aes(genomic_start, log2FoldChange, label = gene_name)) + 
  # up taar red
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, padj < 0.05 & log2FoldChange > 0.585), size = 2, pch = 16, na.rm = T, color = "red")+ 
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, log2FoldChange >= 0), size = 2, pch = 1, na.rm = T, color = "red") + 
  geom_text_repel(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, padj < 0.05 & log2FoldChange > 0.585), size = 4, colour = 'red', nudge_x = -0.25, nudge_y = -0.25)+ 
  geom_text_repel(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, log2FoldChange >= 0), size = 4, colour = 'red', nudge_x = -0.25, nudge_y = -0.25) + 
  
  # down taar blue
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, padj < 0.05 & log2FoldChange < -0.585), size = 2, pch = 16, na.rm = T, color = "blue")+ 
  geom_point(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, log2FoldChange < 0), size = 2, pch = 1, na.rm = T, color = "blue")+ 
  geom_text_repel(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, padj < 0.05 & log2FoldChange < -0.585), size = 4, colour = 'blue', nudge_x = -0.25, nudge_y = -0.25)+ 
  geom_text_repel(data = subset(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, log2FoldChange < 0), size = 4, colour = 'blue', nudge_x = -0.25, nudge_y = -0.25) + 
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*0.9, ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.3))+
  scale_x_continuous(expand = c(0, 0), limits = c(23790000, 24310000), breaks = seq(23800000, 24300000, by = 100000) , labels = seq(23800000, 24300000, by = 100000)/1000000)+
  
  
  theme_classic()+ 
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Position on chromosome 10(MB)", y = "log2 Gene expression changes") +
  
  geom_rect(aes(xmin = genomic_start, xmax = genomic_end, ymin = ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange)), ymax = ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.1), colour = "black", fill = "black")+
  
  geom_text(data = data.frame(label = c("Slc18b1", "Taar", "Stx7", "Moxd1"), 
                              x = c(23815000, 23900000, 24168000, 24260000), 
                              y = rep(ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.2, 4)), 
            mapping = aes(x = x, y = y, label = label))+
  geom_text(data = data.frame(label = gsub("Taar", "", taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name[grepl("Taar",taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name)]), 
                              x = ((taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_start+taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_end)/2)[grepl("Taar",taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_name)], 
                              y = rep(ceiling(max(taar_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*1.2, 14)), 
            mapping = aes(x = x, y = y, label = label), 
            size = 2
  )

dev.off()
##### plot 5: log2FC with genome location (OR and Taar cluster)####
res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- read.csv("Results/DEseq2_results_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9.csv", row.names = 1, header = T)
receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9[grepl(pattern = "Olfr|Taar", rownames(res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9)) & grepl("protein_coding", res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$gene_type), ]
receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 <- receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9 %>% 
  mutate(genomic_start_relative = 1:nrow(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9))
write.csv(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, file = "Results/olfactory_receptor_res_OverExpression.csv")

# Generate genome location data
genome_location <- as.data.frame(table(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$chr_name))
genome_location <- genome_location[match(unique(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$chr_name), genome_location$Var1), ]
genome_location$start <- 0
genome_location <- genome_location %>%
  mutate(end = cumsum(Freq), 
         start = if_else(row_number() == 1, start, start + lag(end)),
         Var1 = factor(Var1, levels = paste0("chr", c(1:19, "X", "Y")))) 

col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)
col_for_chr <- col_for_chr_universal 

# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y")) 
col_for_chr <- col_for_chr[match(chr_all, unique(genome_location$Var1))]
col_for_chr <- col_for_chr[!is.na(col_for_chr)]

pdf("./Plots/receptors_log2FC_with_genome_location.pdf", width = 8, height = 7)
ggplot(data = receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  # OR red
  geom_point(data = subset(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, grepl("^Olfr", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  
  # Taar blue
  geom_point(data = subset(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, grepl("^Taar", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)) , size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, grepl("^Taar", gene_name) & padj < 0.05 &  abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") +
  
  geom_text_repel(data = subset(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, 
            aes(xmin = start, xmax = end, ymin =  ceiling(max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*0.95 , ymax =  ceiling(max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))*0.98, fill = Var1), 
            inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange))-1, ceiling(max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$log2FoldChange)))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_end_relative)), breaks = seq(0, max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_end_relative), by = 200), labels = seq(0, max(receptor_res_OverExpression_exclude_lsl_Tbr1_tdTomato_het_9$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes")
dev.off()

