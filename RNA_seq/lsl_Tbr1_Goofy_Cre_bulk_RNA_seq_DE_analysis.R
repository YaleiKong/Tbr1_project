Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list=ls())
getwd()  

## set seed for repeat
set.seed(1)

## library 
library_Ls <- c("DESeq2","LSD","gProfileR",
                "stringr","readxl",
                "dplyr","tidyr","tidyverse",
                "ggplot2","pheatmap","ggrepel","gg.gap","RColorBrewer","grid",
                "gridExtra","gplots","reshape2","ggpubr","ggbreak","gtools","textshape",
                "circlize","ComplexHeatmap","hutils"
)

lapply(library_Ls,function(x){suppressPackageStartupMessages(library(x,character.only = T))})
library(fdrtool)
# if too few (or none) differential expressed genes, fdrtools can fix the null hypothesis in DESeq2.
# https://seqqc.wordpress.com/2016/01/27/too-few-or-none-differential-expressed-genes-a-way-to-fix-the-null-hypothesis-in-deseq2/


## set the workpath
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/7_lsl_Tbr1_Goofy_Cre_bulk_RNA_Seq/20250705_lsl_Tbr1_Goofy_Cre_bulk_RNA_Seq_2507/DE_analysis_20250705/")
setwd(workpath)

## create a result folder
folder <- "Results/"
if (!file.exists(folder)) { dir.create(folder) }

folder_plot <- "Plots/"
if (!file.exists(folder_plot)) { dir.create(folder_plot) }

#### data preparation ####
## Load data
dat <- read.table("../featureCounts/featureCounts_geneId.txt", header=T, row.names=1, as.is=T)
# as.is=T means keep all character vectors, do not convert to factors
# dat_Cre_tdTomato <- read.table("featureCounts/Cre_tdTomato_featureCounts.txt", header=T, row.names=1, as.is=T) %>% 
#   t()

## load ensembl gene name and gene description file
# import ensembl info and remove the duplicated gene_name
ensembl_gene_list_good <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) %>% 
  mutate(gene_id_Ensembl_remove_version_meta = str_remove(gene_id_Ensembl, "\\.[0-9]")) %>% 
  filter(gene_type %in%  c("protein_coding", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "polymorphic_pseudogene")) %>% 
  filter(grepl("^chr[1-9XY]", chr_name)) %>% 
  filter(!duplicated(gene_name))
dim(ensembl_gene_list_good)
# [1] 24899    15

unique(ensembl_gene_list_good[grep("Olfr", ensembl_gene_list_good$gene_name),"gene_type"])
# [1] "protein_coding"                     "unprocessed_pseudogene"             "transcribed_unprocessed_pseudogene" "polymorphic_pseudogene" 

## match the ensembl ID of the data to the ensembl ID from ensembl_gene_list_good
good_dat <- dat[rownames(dat) %in% ensembl_gene_list_good$gene_id_Ensembl, ]

## obtain the gene names
good_dat_matched_gene_name_with_gene_length <- good_dat %>%
  rownames_to_column(var = "gene_id") %>% 
  mutate(gene_name = ensembl_gene_list_good[match(gene_id, ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name) %>% 
  filter(!duplicated(gene_name)) %>%
  dplyr::rename(gene_length = Length) 
good_dat_matched_gene_name_with_gene_length <- column_to_rownames(good_dat_matched_gene_name_with_gene_length, "gene_name")
dim(good_dat_matched_gene_name_with_gene_length)

###### sample info #####
meta <- data.frame(
  sample = paste0("X",1:16),
  sample_name = c(paste0("lsl_Tbr1_tdTomato_Het_", 1:4), paste0("lsl_Tbr1_tdTomato_Het_Goofy_cre_", 1:4),
                  paste0("lsl_Tbr1_tdTomato_Homo_", 1:4), paste0("lsl_Tbr1_tdTomato_Homo_Goofy_cre_",1:4))
)
meta <- meta %>% 
  mutate(genotype = str_remove(sample_name, "_[1-4]$"))

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
  select(meta$sample_name)

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
dds_dat <- DESeqDataSetFromMatrix(countData = good_dat_matched_gene_name, colData = meta, design= ~ genotype)
dds_dat$genotype <- relevel(dds_dat$genotype, "lsl_Tbr1_tdTomato_Het")  

##### PCA #####
vsd <- vst(dds_dat, blind = FALSE)
DESeq2::plotPCA(vsd, intgroup = "genotype")

# ggplot2
pca_data <- DESeq2::plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p.pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(lsl_Tbr1_tdTomato_Het_Goofy_cre = "orangered",
                                lsl_Tbr1_tdTomato_Het = "springgreen",
                                lsl_Tbr1_tdTomato_Homo_Goofy_cre = "deeppink2",
                                lsl_Tbr1_tdTomato_Homo = "aquamarine"))+
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

##### DEseq2 #####
dds_dat <- DESeq(dds_dat)

##### export the DESeq2 results #####
###### LRT #####
res_lrt <- results(dds_dat, independentFiltering = FALSE) %>% 
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
hist(res_lrt$pvalue)
write.csv(res_lrt, "Results/DEseq2_results_lrt.csv")

###### Het #####
res_OverExpression_Het <- results(dds_dat, contrast = c("genotype", "lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Het"), test="Wald", independentFiltering = FALSE) %>% 
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
write.csv(res_OverExpression_Het, "Results/DEseq2_results_OverExpression_Het.csv")
hist(res_OverExpression_Het$pvalue)
res_OverExpression_Het[res_OverExpression_Het$padj < 0.05 & abs(res_OverExpression_Het$log2FoldChange) > 0.585,]
nrow(res_OverExpression_Het[res_OverExpression_Het$padj < 0.05 & abs(res_OverExpression_Het$log2FoldChange) > 0.585,]) 
# 16
table(res_OverExpression_Het$padj < 0.05)
# FALSE  TRUE 
# 21730    24 

###### Homo #####
res_OverExpression_Homo <- results(dds_dat, contrast = c("genotype", "lsl_Tbr1_tdTomato_Homo_Goofy_cre", "lsl_Tbr1_tdTomato_Homo"), test="Wald", independentFiltering = FALSE) %>% 
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
write.csv(res_OverExpression_Homo, "Results/DEseq2_results_OverExpression_Homo.csv")
hist(res_OverExpression_Homo$pvalue)
res_OverExpression_Homo[res_OverExpression_Homo$padj < 0.05 & abs(res_OverExpression_Homo$log2FoldChange) > 0.585,]
nrow(res_OverExpression_Homo[res_OverExpression_Homo$padj < 0.05 & abs(res_OverExpression_Homo$log2FoldChange) > 0.585,]) 
# 579
table(res_OverExpression_Homo$padj < 0.05)
# FALSE  TRUE 
# 18596  3158

###### Homo vs Het #####
res_OverExpression_Homo_Het <- results(dds_dat, contrast = c("genotype", "lsl_Tbr1_tdTomato_Homo_Goofy_cre", "lsl_Tbr1_tdTomato_Het_Goofy_cre"), test="Wald", independentFiltering = FALSE) %>% 
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
write.csv(res_OverExpression_Homo_Het, "Results/DEseq2_results_OverExpression_Homo_Het.csv")

###### Homo vs Het-control #####
res_OverExpression_Homo_vs_Het_control <- results(dds_dat, contrast = c("genotype", "lsl_Tbr1_tdTomato_Homo_Goofy_cre", "lsl_Tbr1_tdTomato_Het"), test="Wald", independentFiltering = FALSE) %>% 
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
write.csv(res_OverExpression_Homo_vs_Het_control, "Results/DEseq2_results_OverExpression_Homo_vs_Het_control.csv")

##### export the normalized counts #####
normalized_count <- DESeq2::counts(dds_dat, normalized = TRUE) %>% 
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
write.csv(normalized_count, "Results/DEseq2_normalized_count.csv")

#### plot: volcano plot (padj) #####
##### res_OverExpression_Het (padj) #####
res_OverExpression_Het <- read.csv("Results/DEseq2_results_OverExpression_Het.csv", header = T, row.names = 1) %>% 
  as.data.frame()
range(res_OverExpression_Het$log2FoldChange, na.rm = T) 
# [1] -23.790721   4.860511
range(-log10(res_OverExpression_Het$padj), na.rm = T)
#  0.00000 59.41825
sum(-log10(res_OverExpression_Het$padj)>30)
# 1
sum(res_OverExpression_Het$log2FoldChange < (-11))
# 1 
res_OverExpression_Het <- res_OverExpression_Het %>%
  mutate(
    # 如果 padj 小于 1e-30（包括0），都设为 1e-30
    padj = ifelse(padj <= 1e-30, 1e-30, padj),
    # 限制 log2FoldChange 在 [-11, 11] 之间
    log2FoldChange = pmin(pmax(log2FoldChange, -11), 11)
  )

res_OverExpression_Het_functional_OR <- res_OverExpression_Het[which(grepl("Olfr", rownames(res_OverExpression_Het)) &
                                                                       res_OverExpression_Het$gene_type %in% "protein_coding"), ]
res_OverExpression_Het_functional_Taar<- res_OverExpression_Het[which(grepl("Taar[2-9]", rownames(res_OverExpression_Het)) & 
                                                                        res_OverExpression_Het$gene_type %in% "protein_coding"), ]
res_OverExpression_Het_functional_OR_diff <- res_OverExpression_Het_functional_OR[which(res_OverExpression_Het_functional_OR$padj < 0.05  & 
                                                                                          abs(res_OverExpression_Het_functional_OR$log2FoldChange) > 0.585) ,]
dim(res_OverExpression_Het_functional_OR_diff)

res_OverExpression_Het_functional_Taar_diff <- res_OverExpression_Het_functional_Taar[which(res_OverExpression_Het_functional_Taar$padj < 0.05  & 
                                                                                              abs(res_OverExpression_Het_functional_Taar$log2FoldChange) > 0.585) ,]
res_OverExpression_Het_rm_functional_OR_diff_functional_Taar_diff <- subset(res_OverExpression_Het, !(gene_name %in% c(res_OverExpression_Het_functional_OR_diff$gene_name, res_OverExpression_Het_functional_Taar_diff$gene_name)))

res_OverExpression_Het_volcano_plot <- res_OverExpression_Het %>%
  mutate(group = case_when(
    gene_name %in% res_OverExpression_Het_functional_OR_diff$gene_name ~ "DE ORs",
    gene_name %in% res_OverExpression_Het_functional_Taar_diff$gene_name ~ "DE Taars",
    padj < 0.05 & abs(log2FoldChange) > 0.585 ~ "Significant DEGs",
    TRUE ~ "Other"
  ))

write.csv(res_OverExpression_Het_volcano_plot, file = "Results/res_OverExpression_Het_with_label.csv")

pdf("./Plots/volcano_plot_OverExpression_Het_DEGs_padj.pdf", height = 8, width = 8)
ggplot() +
  # 先画 "Other" 和 "Significant DEGs" (灰色 + 黑色)
  geom_point(
    data = subset(res_OverExpression_Het_volcano_plot, group %in% c("Other", "Significant DEGs")),
    aes(x = log2FoldChange, y = -log10(padj), color = group, fill = group),
    alpha = 0.75, size = 1
  ) +
  # 再画重点 (蓝色 ORs, 红色 Taars, 绿色 Gucy1b2)
  geom_point(
    data = subset(res_OverExpression_Het_volcano_plot, group %in% c("DE ORs")),
    aes(x = log2FoldChange, y = -log10(padj), color = group),
    size = 2
  ) +
  geom_point(
    data = subset(res_OverExpression_Het_volcano_plot, group %in% "DE Taars"),
    aes(x = log2FoldChange, y = -log10(padj), color = group),
    size = 2
  ) +
  # 标签
  geom_text_repel(
    data = subset(res_OverExpression_Het_volcano_plot, group %in% c("DE ORs", "DE Taars")),
    aes(x = log2FoldChange, y = -log10(padj), label = gene_name, color = group),
    size = 3, box.padding = 0.3, max.overlaps = 30, show.legend = FALSE
  ) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c(
      "Other" = "grey",
      "Significant DEGs" = "black",
      "DE ORs" = "blue",
      "DE Taars" = "red",
      "DE Gucy1b2" = "green"
    )
  ) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(-11, 11), oob = scales::squish) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30), oob = scales::squish) +
  ggtitle("res_OverExpression_Het padj")
dev.off()

##### res_OverExpression_Homo (padj) #####
res_OverExpression_Homo <- read.csv("Results/DEseq2_results_OverExpression_Homo.csv", header = T, row.names = 1) %>% 
  as.data.frame()
range(res_OverExpression_Homo$log2FoldChange, na.rm = T) 
# -4.348601 10.255184
range(-log10(res_OverExpression_Homo$padj), na.rm = T)
# 0.00000 65.87927

sum(-log10(res_OverExpression_Homo$padj)>30)
# 2
sum(res_OverExpression_Homo$log2FoldChange < (-11))
# 0

res_OverExpression_Homo <- res_OverExpression_Homo %>%
  mutate(
    # 如果 padj 小于 1e-30（包括0），都设为 1e-30
    padj = ifelse(padj <= 1e-30, 1e-30, padj),
    # 限制 log2FoldChange 在 [-11, 11] 之间
    log2FoldChange = pmin(pmax(log2FoldChange, -11), 11)
  )

res_OverExpression_Homo_functional_OR <- res_OverExpression_Homo[which(grepl("Olfr", rownames(res_OverExpression_Homo)) & res_OverExpression_Homo$gene_type %in% "protein_coding"), ]
res_OverExpression_Homo_functional_Taar<- res_OverExpression_Homo[which(grepl("Taar[2-9]", rownames(res_OverExpression_Homo)) & res_OverExpression_Homo$gene_type %in% "protein_coding"), ]
res_OverExpression_Homo_functional_OR_diff <- res_OverExpression_Homo_functional_OR[which(res_OverExpression_Homo_functional_OR$padj < 0.05  & abs(res_OverExpression_Homo_functional_OR$log2FoldChange) > 0.585) ,]
dim(res_OverExpression_Homo_functional_OR_diff)

res_OverExpression_Homo_functional_Taar_diff <- res_OverExpression_Homo_functional_Taar[which(res_OverExpression_Homo_functional_Taar$padj < 0.05  & abs(res_OverExpression_Homo_functional_Taar$log2FoldChange) > 0.585) ,]
res_OverExpression_Homo_rm_functional_OR_diff_functional_Taar_diff <- subset(res_OverExpression_Homo, !(gene_name %in% c(res_OverExpression_Homo_functional_OR_diff$gene_name, res_OverExpression_Homo_functional_Taar_diff$gene_name)))

res_OverExpression_Homo_volcano_plot <- res_OverExpression_Homo %>%
  mutate(group = case_when(
    gene_name %in% res_OverExpression_Homo_functional_OR_diff$gene_name ~ "DE ORs",
    gene_name %in% res_OverExpression_Homo_functional_Taar_diff$gene_name ~ "DE Taars",
    padj < 0.05 & abs(log2FoldChange) > 0.585 ~ "Significant DEGs",
    TRUE ~ "Other"
  ))
unique(res_OverExpression_Homo_volcano_plot$group)
write.csv(res_OverExpression_Homo_volcano_plot, file = "Results/res_OverExpression_Homo_with_label.csv")

pdf("./Plots/volcano_plot_OverExpression_Homo_DEGs_padj.pdf", height = 8, width = 8)
ggplot() +
  # 先画 "Other" 和 "Significant DEGs" (灰色 + 黑色)
  geom_point(
    data = subset(res_OverExpression_Homo_volcano_plot, group %in% c("Other", "Significant DEGs")),
    aes(x = log2FoldChange, y = -log10(padj), color = group, fill = group),
    alpha = 0.75, size = 1
  ) +
  # 再画重点 (蓝色 ORs, 红色 Taars, 绿色 Gucy1b2)
  geom_point(
    data = subset(res_OverExpression_Homo_volcano_plot, group %in% c("DE ORs")),
    aes(x = log2FoldChange, y = -log10(padj), color = group),
    size = 2
  ) +
  geom_point(
    data = subset(res_OverExpression_Homo_volcano_plot, group %in% "DE Taars"),
    aes(x = log2FoldChange, y = -log10(padj), color = group),
    size = 2
  ) +
  # 标签
  geom_text_repel(
    data = subset(res_OverExpression_Homo_volcano_plot, group %in% c("DE ORs", "DE Taars")),
    aes(x = log2FoldChange, y = -log10(padj), label = gene_name, color = group),
    size = 3, box.padding = 0.3, max.overlaps = 30, show.legend = FALSE
  ) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c(
      "Other" = "grey",
      "Significant DEGs" = "black",
      "DE ORs" = "blue",
      "DE Taars" = "red",
      "DE Gucy1b2" = "green"
    )
  ) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(-11, 11), oob = scales::squish) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30), oob = scales::squish) +
  ggtitle("res_OverExpression_Homo padj")
dev.off()

#### plot: Taar normalized expression counts bar plot #####
normalized_count <- read.csv("./Results/DEseq2_normalized_count.csv", header = T, row.names = 1)
##### Het padj #####
Taar_normalized_count_Het <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                              grepl("Het", colnames(normalized_count))]

Taar_normalized_count_Het_for_plot <- Taar_normalized_count_Het %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count_Het)[grep("Het", colnames(Taar_normalized_count_Het))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_Het_Goofy_cre", sampletype), 
                           "lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Het")) %>%  
  left_join(res_OverExpression_Het_functional_Taar[grepl("gene_name|padj|fdrtool_qval", colnames(res_OverExpression_Het_functional_Taar))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Het","lsl_Tbr1_tdTomato_Het_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count_Het))) 

pdf("./Plots/bar_plot_Taar_normalized_count_Het_padj.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_Het_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_Het_for_plot,  sampletype == "lsl_Tbr1_tdTomato_Het_Goofy_cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_Het_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_Het_for_plot$normalized_count))*1.05))+ 
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

##### Homo padj #####
Taar_normalized_count_Homo <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                               grepl("Homo", colnames(normalized_count))]

Taar_normalized_count_Homo_for_plot <- Taar_normalized_count_Homo %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count_Homo)[grep("Homo", colnames(Taar_normalized_count_Homo))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_Homo_Goofy_cre", sampletype), 
                           "lsl_Tbr1_tdTomato_Homo_Goofy_cre", "lsl_Tbr1_tdTomato_Homo")) %>%  
  left_join(res_OverExpression_Homo_functional_Taar[grepl("gene_name|padj|fdrtool_qval", colnames(res_OverExpression_Homo_functional_Taar))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Homo","lsl_Tbr1_tdTomato_Homo_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count_Homo))) 

pdf("./Plots/bar_plot_Taar_normalized_count_Homo_padj.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_Homo_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_Homo_for_plot,  sampletype == "lsl_Tbr1_tdTomato_Homo_Goofy_cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_Homo_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_Homo_for_plot$normalized_count))*1.05))+ 
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

##### Het padj scale #####
Taar_normalized_count_Het <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                              grepl("Het", colnames(normalized_count))]

Taar_normalized_count_Het_scale_for_plot <- as.data.frame(Taar_normalized_count_Het/rowMeans(Taar_normalized_count_Het[, grep("lsl_Tbr1_tdTomato_Het_\\d",colnames(Taar_normalized_count_Het))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count_Het)[grep("Het", colnames(Taar_normalized_count_Het))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_Het_Goofy_cre", sampletype), 
                           "lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Het")) %>%  
  left_join(res_OverExpression_Het_functional_Taar[grepl("gene_name|padj|fdrtool_qval", colnames(res_OverExpression_Het_functional_Taar))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Het","lsl_Tbr1_tdTomato_Het_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count_Het))) 

pdf("./Plots/bar_plot_Taar_normalized_count_Het_scale_padj.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_Het_scale_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_Het_scale_for_plot,  sampletype == "lsl_Tbr1_tdTomato_Het_Goofy_cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_Het_scale_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_Het_scale_for_plot$normalized_count))*1.05))+ 
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
  labs(x = "", y = "normalized counts compared to lsl_Tbr1_tdTomato_Het") +
  ggtitle("bar_plot_Taar_normalized_count_Het_padj")
dev.off()

##### Homo padj scale #####
Taar_normalized_count_Homo <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                               grepl("Homo", colnames(normalized_count))]

Taar_normalized_count_Homo_scale_for_plot <- as.data.frame(Taar_normalized_count_Homo/rowMeans(Taar_normalized_count_Homo[, grep("lsl_Tbr1_tdTomato_Homo_\\d",colnames(Taar_normalized_count_Homo))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count_Homo)[grep("Homo", colnames(Taar_normalized_count_Homo))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = ifelse(grepl("lsl_Tbr1_tdTomato_Homo_Goofy_cre", sampletype), 
                           "lsl_Tbr1_tdTomato_Homo_Goofy_cre", "lsl_Tbr1_tdTomato_Homo")) %>%  
  left_join(res_OverExpression_Homo_functional_Taar[grepl("gene_name|padj|fdrtool_qval", colnames(res_OverExpression_Homo_functional_Taar))], by = "gene_name")  %>%
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Homo","lsl_Tbr1_tdTomato_Homo_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count_Homo))) 

pdf("./Plots/bar_plot_Taar_normalized_count_Homo_scale_padj.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_Homo_scale_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_Homo_scale_for_plot,  sampletype == "lsl_Tbr1_tdTomato_Homo_Goofy_cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_Homo_scale_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_Homo_scale_for_plot$normalized_count))*1.05))+ 
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
  labs(x = "", y = "normalized counts compared to lsl_Tbr1_tdTomato_Homo") +
  ggtitle("bar_plot_Taar_normalized_count_Homo_padj")
dev.off()

#### plot: log2FC with genome location (OR and Taar)####
##### Het #####
res_OverExpression_Het <- read.csv("Results/DEseq2_results_OverExpression_Het.csv", row.names = 1, header = T)
receptor_res_OverExpression_Het <- res_OverExpression_Het[grepl(pattern = "Olfr|Taar", rownames(res_OverExpression_Het)) & grepl("protein_coding", res_OverExpression_Het$gene_type), ]
range(receptor_res_OverExpression_Het$log2FoldChange)
# [1] -1.805102  2.854058
receptor_res_OverExpression_Het <- receptor_res_OverExpression_Het %>% 
  mutate(genomic_start_relative = 1:nrow(receptor_res_OverExpression_Het),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(receptor_res_OverExpression_Het))
write.csv(receptor_res_OverExpression_Het, file = "Results/olfactory_receptor_res_OverExpression_Het.csv")

# Generate genome location data
genome_location <- as.data.frame(table(receptor_res_OverExpression_Het$chr_name))
genome_location <- genome_location[match(unique(receptor_res_OverExpression_Het$chr_name), genome_location$Var1), ]
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

pdf("./Plots/receptors_log2FC_with_genome_location_OverExpression_Het.pdf", width = 8, height = 7)
ggplot(data = receptor_res_OverExpression_Het, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  # OR red
  geom_point(data = subset(receptor_res_OverExpression_Het, grepl("^Olfr", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(receptor_res_OverExpression_Het, grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  
  # Taar OR blue
  geom_point(data = subset(receptor_res_OverExpression_Het, grepl("^Taar", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(receptor_res_OverExpression_Het, grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") + 
  
  geom_text_repel(data = subset(receptor_res_OverExpression_Het, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin = 3.5, ymax = 4, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(-4, 4)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(receptor_res_OverExpression_Het$genomic_end_relative)), breaks = seq(0, max(receptor_res_OverExpression_Het$genomic_end_relative), by = 200), labels = seq(0, max(receptor_res_OverExpression_Het$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes") +
  ggtitle("OverExpression_Het")
dev.off()

##### Homo #####
res_OverExpression_Homo <- read.csv("Results/DEseq2_results_OverExpression_Homo.csv", row.names = 1, header = T)
receptor_res_OverExpression_Homo <- res_OverExpression_Homo[grepl(pattern = "Olfr|Taar", rownames(res_OverExpression_Homo)) & grepl("protein_coding", res_OverExpression_Homo$gene_type), ]
range(receptor_res_OverExpression_Homo$log2FoldChange)
# [1] -2.889217  3.584960
receptor_res_OverExpression_Homo <- receptor_res_OverExpression_Homo %>% 
  mutate(genomic_start_relative = 1:nrow(receptor_res_OverExpression_Homo),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(receptor_res_OverExpression_Homo))
write.csv(receptor_res_OverExpression_Homo, file = "Results/olfactory_receptor_res_OverExpression_Homo.csv")

# Generate genome location data
genome_location <- as.data.frame(table(receptor_res_OverExpression_Homo$chr_name))
genome_location <- genome_location[match(unique(receptor_res_OverExpression_Homo$chr_name), genome_location$Var1), ]
genome_location$start <- 0
genome_location <- genome_location %>%
  mutate(end = cumsum(Freq), 
         start = if_else(row_number() == 1, start, start + lag(end)),
         Var1 = factor(Var1, levels = paste0("chr", c(1:19, "X", "Y")) )) 
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)
col_for_chr <- col_for_chr_universal 

# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y")) 
col_for_chr <- col_for_chr[match(chr_all, unique(genome_location$Var1))]
col_for_chr <- col_for_chr[!is.na(col_for_chr)]

pdf("./Plots/receptors_log2FC_with_genome_location_OverExpression_Homo.pdf", width = 8, height = 7)
ggplot(data = receptor_res_OverExpression_Homo, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  # OR red
  geom_point(data = subset(receptor_res_OverExpression_Homo, grepl("^Olfr", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(receptor_res_OverExpression_Homo, grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  # Taar OR blue
  geom_point(data = subset(receptor_res_OverExpression_Homo, grepl("^Taar", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(receptor_res_OverExpression_Homo, grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") + 
  
  geom_text_repel(data = subset(receptor_res_OverExpression_Homo, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin = 3.5, ymax = 4, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(-4, 4)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(receptor_res_OverExpression_Homo$genomic_end_relative)), breaks = seq(0, max(receptor_res_OverExpression_Homo$genomic_end_relative), by = 200), labels = seq(0, max(receptor_res_OverExpression_Homo$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes") +
  ggtitle("OverExpression_Homo")
dev.off()

#### plot: Taar normalized count barplot ####
normalized_count <- read.csv("./Results/DEseq2_normalized_count.csv", header = T, row.names = 1)
Taar_normalized_count <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                          grepl("lsl", colnames(normalized_count))]
res_lrt <- read.csv("Results/DEseq2_results_lrt.csv", header = T, row.names = 1)
res_lrt_functional_Taar <- res_lrt[which(grepl(pattern = "^Taar[2-9]", rownames(res_lrt)) & res_lrt$gene_type %in% "protein_coding"),]

res_OverExpression_Het <- read.csv("Results/DEseq2_results_OverExpression_Het.csv", header = T, row.names = 1)
res_Het_functional_Taar <- res_OverExpression_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Het)) & res_OverExpression_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo <- read.csv("Results/DEseq2_results_OverExpression_Homo.csv", header = T, row.names = 1)
res_Homo_functional_Taar <- res_OverExpression_Homo[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo)) & res_OverExpression_Homo$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_Het <- read.csv("Results/DEseq2_results_OverExpression_Homo_Het.csv", header = T, row.names = 1)
res_Homo_Het_functional_Taar <- res_OverExpression_Homo_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_Het)) & res_OverExpression_Homo_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_vs_Het_control <- read.csv("Results/DEseq2_results_OverExpression_Homo_vs_Het_control.csv", header = T, row.names = 1)
res_Homo_vs_Het_control_functional_Taar <- res_OverExpression_Homo_vs_Het_control[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_vs_Het_control)) & res_OverExpression_Homo_vs_Het_control$gene_type %in% "protein_coding"),]

Taar_normalized_count_scale_for_plot <- as.data.frame(Taar_normalized_count/rowMeans(Taar_normalized_count[, grep("lsl_Tbr1_tdTomato_Het_\\d",colnames(Taar_normalized_count))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count)[grep("lsl", colnames(Taar_normalized_count))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = sub("_[0-9]+$", "", sampletype)) %>%  
  left_join(res_lrt_functional_Taar[grepl("gene_name|padj", colnames(res_lrt_functional_Taar))], by = "gene_name")  %>%
  mutate(Het_padj=res_Het_functional_Taar$padj[match(gene_name, rownames(res_Het_functional_Taar))],
         Homo_padj=res_Homo_functional_Taar$padj[match(gene_name, rownames(res_Homo_functional_Taar))],
         Homo_Het_padj=res_Homo_Het_functional_Taar$padj[match(gene_name, rownames(res_Homo_Het_functional_Taar))],
         Homo_vs_Het_control_padj=res_Homo_vs_Het_control_functional_Taar$padj[match(gene_name, rownames(res_Homo_vs_Het_control_functional_Taar))]) %>% 
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Het","lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Homo", "lsl_Tbr1_tdTomato_Homo_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))) 
write.csv(Taar_normalized_count_scale_for_plot, file = "Results/Taar_normalized_counts_with_stats.csv")

# Generate comparisons_all with dynamic y-position
comparisons_all <- Taar_normalized_count_scale_for_plot %>%
  group_by(gene_name) %>%
  summarise(max_val = max(normalized_count, na.rm = TRUE)) %>%
  rowwise() %>%
  do({
    gene <- .$gene_name
    maxy <- .$max_val
    
    # Extract p-values for the gene with error handling
    gene_data <- Taar_normalized_count_scale_for_plot %>% 
      filter(gene_name == gene) %>% 
      slice_head(n = 1)  # Use slice_head for clarity
    
    # Define comparisons and corresponding p-values
    data.frame(
      gene_name = gene,
      group1 = c("lsl_Tbr1_tdTomato_Het",
                 "lsl_Tbr1_tdTomato_Homo",
                 "lsl_Tbr1_tdTomato_Het_Goofy_cre",
                 "lsl_Tbr1_tdTomato_Het"),
      group2 = c("lsl_Tbr1_tdTomato_Het_Goofy_cre",
                 "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
                 "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
                 "lsl_Tbr1_tdTomato_Homo_Goofy_cre"),
      p.adj = c(
        ifelse(nrow(gene_data) > 0 && !is.na(gene_data$Het_padj), gene_data$Het_padj, NA),
        ifelse(nrow(gene_data) > 0 && !is.na(gene_data$Homo_padj), gene_data$Homo_padj, NA),
        ifelse(nrow(gene_data) > 0 && !is.na(gene_data$Homo_Het_padj), gene_data$Homo_Het_padj, NA),
        ifelse(nrow(gene_data) > 0 && !is.na(gene_data$Homo_vs_Het_control_padj), gene_data$Homo_vs_Het_control_padj, NA)
      ),
      y.position = seq(maxy, maxy * 1.4, length.out = 4)  # Dynamic spacing
    )
  }) %>% 
  ungroup() %>%
  filter(!is.na(p.adj)) %>%  # Remove invalid p-values
  mutate(p.adj = format(p.adj, scientific = TRUE, digits = 3))  # Format p-values for readability


##### Create a list to store individual plots ####
plot_list <- list()

for(gene in unique(Taar_normalized_count_scale_for_plot$gene_name)) {
  
  plot_df <- Taar_normalized_count_scale_for_plot %>% filter(gene_name == gene)
  comp_df <- comparisons_all %>% filter(gene_name == gene)
  
  p <- ggplot(plot_df, aes(x = genotype, y = normalized_count)) +
    
    # Bar plot for mean
    stat_summary(
      aes(fill = genotype),
      fun = mean,
      geom = "bar",
      color = "black",
      linewidth = 0.4,
      width = 0.8,
      alpha = 0.85
    ) +
    
    # Error bars
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.25
    ) +
    
    # Jittered sample points
    geom_point(
      position = position_jitter(width = 0.1),
      shape = 21,
      fill = "white",
      color = "black",
      size = 2,
      stroke = 0.6
    ) +
    
    # Add p-values
    stat_pvalue_manual(
      comp_df,
      label = "p.adj",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      inherit.aes = FALSE,
      tip.length = 0.02,
      bracket.size = 0.5,
      size = 3.5
    ) +
    
    scale_fill_manual(values = c(
      "lsl_Tbr1_tdTomato_Het" = "#66C2A5", 
      "lsl_Tbr1_tdTomato_Het_Goofy_cre" = "#1B9E77",
      "lsl_Tbr1_tdTomato_Homo" =  "#8DA0CB",
      "lsl_Tbr1_tdTomato_Homo_Goofy_cre" = "#7570B3"
    )) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_df$normalized_count))*1.5) +
    
    
    theme_classic() +
    theme(axis.text.x = element_text(size = 3)) +
    
    labs(
      title = gene,
      x = "Genotype",
      y = "Normalized Counts"
    )
  
  # Save to list
  plot_list[[gene]] <- p
}

# Display a single plot as example
print(plot_list[[1]])

pdf("Plots/Taars_normalized_count_barplots.pdf", width = 6, height = 5)

for(gene in names(plot_list)){
  print(plot_list[[gene]])
}

dev.off()

##### Taar normalized count barplot (ylim = 16 exclude Taar6) ####
normalized_count <- read.csv("./Results/DEseq2_normalized_count.csv", header = T, row.names = 1)
Taar_normalized_count <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                          grepl("lsl", colnames(normalized_count))]
res_lrt <- read.csv("Results/DEseq2_results_lrt.csv", header = T, row.names = 1)
res_lrt_functional_Taar <- res_lrt[which(grepl(pattern = "^Taar[2-9]", rownames(res_lrt)) & res_lrt$gene_type %in% "protein_coding"),]

res_OverExpression_Het <- read.csv("Results/DEseq2_results_OverExpression_Het.csv", header = T, row.names = 1)
res_Het_functional_Taar <- res_OverExpression_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Het)) & res_OverExpression_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo <- read.csv("Results/DEseq2_results_OverExpression_Homo.csv", header = T, row.names = 1)
res_Homo_functional_Taar <- res_OverExpression_Homo[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo)) & res_OverExpression_Homo$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_Het <- read.csv("Results/DEseq2_results_OverExpression_Homo_Het.csv", header = T, row.names = 1)
res_Homo_Het_functional_Taar <- res_OverExpression_Homo_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_Het)) & res_OverExpression_Homo_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_vs_Het_control <- read.csv("Results/DEseq2_results_OverExpression_Homo_vs_Het_control.csv", header = T, row.names = 1)
res_Homo_vs_Het_control_functional_Taar <- res_OverExpression_Homo_vs_Het_control[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_vs_Het_control)) & res_OverExpression_Homo_vs_Het_control$gene_type %in% "protein_coding"),]

Taar_normalized_count_scale_for_plot <- as.data.frame(Taar_normalized_count/rowMeans(Taar_normalized_count[, grep("lsl_Tbr1_tdTomato_Het_\\d",colnames(Taar_normalized_count))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count)[grep("lsl", colnames(Taar_normalized_count))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = sub("_[0-9]+$", "", sampletype)) %>%  
  left_join(res_lrt_functional_Taar[grepl("gene_name|padj", colnames(res_lrt_functional_Taar))], by = "gene_name")  %>%
  mutate(Het_padj=res_Het_functional_Taar$padj[match(gene_name, rownames(res_Het_functional_Taar))],
         Homo_padj=res_Homo_functional_Taar$padj[match(gene_name, rownames(res_Homo_functional_Taar))],
         Homo_Het_padj=res_Homo_Het_functional_Taar$padj[match(gene_name, rownames(res_Homo_Het_functional_Taar))],
         Homo_vs_Het_control_padj=res_Homo_vs_Het_control_functional_Taar$padj[match(gene_name, rownames(res_Homo_vs_Het_control_functional_Taar))]) %>% 
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Het","lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Homo", "lsl_Tbr1_tdTomato_Homo_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))) 

comparisons_all_exclude_Taar6 <- Taar_normalized_count_scale_for_plot %>%
  group_by(gene_name) %>%
  summarise(max_val = max(normalized_count, na.rm = TRUE)) %>%
  rowwise() %>%
  do({
    gene <- .$gene_name
    
    gene_data <- Taar_normalized_count_scale_for_plot %>% 
      filter(gene_name == gene) %>% 
      slice_head(n = 1)
    
    pvals <- c(
      gene_data$Het_padj,
      gene_data$Homo_padj,
      gene_data$Homo_Het_padj,
      gene_data$Homo_vs_Het_control_padj
    )
    
    # 转换成星号格式
    p_stars <- case_when(
      pvals < 0.001 ~ "***",
      pvals < 0.01  ~ "**",
      pvals < 0.05  ~ "*",
      TRUE ~ "ns"
    )
    
    data.frame(
      gene_name = gene,
      group1 = c(
        "lsl_Tbr1_tdTomato_Het",
        "lsl_Tbr1_tdTomato_Homo",
        "lsl_Tbr1_tdTomato_Het_Goofy_cre",
        "lsl_Tbr1_tdTomato_Het"
      ),
      group2 = c(
        "lsl_Tbr1_tdTomato_Het_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre"
      ),
      pval = pvals,
      p.adj = p_stars,
      y.position = seq(12, 15, length.out = 4) # 固定在y轴之外
    )
  }) %>% 
  ungroup() 

comparisons_all_include_Taar6 <- Taar_normalized_count_scale_for_plot %>%
  group_by(gene_name) %>%
  summarise(max_val = max(normalized_count, na.rm = TRUE)) %>%
  rowwise() %>%
  do({
    gene <- .$gene_name
    
    gene_data <- Taar_normalized_count_scale_for_plot %>% 
      filter(gene_name == gene) %>% 
      slice_head(n = 1)
    
    pvals <- c(
      gene_data$Het_padj,
      gene_data$Homo_padj,
      gene_data$Homo_Het_padj,
      gene_data$Homo_vs_Het_control_padj
    )
    
    # 转换成星号格式
    p_stars <- case_when(
      pvals < 0.001 ~ "***",
      pvals < 0.01  ~ "**",
      pvals < 0.05  ~ "*",
      TRUE ~ "ns"
    )
    
    data.frame(
      gene_name = gene,
      group1 = c(
        "lsl_Tbr1_tdTomato_Het",
        "lsl_Tbr1_tdTomato_Homo",
        "lsl_Tbr1_tdTomato_Het_Goofy_cre",
        "lsl_Tbr1_tdTomato_Het"
      ),
      group2 = c(
        "lsl_Tbr1_tdTomato_Het_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre",
        "lsl_Tbr1_tdTomato_Homo_Goofy_cre"
      ),
      pval = pvals,
      p.adj = p_stars,
      y.position = seq(22, 25, length.out = 4) # 固定在y轴之外
    )
  }) %>% 
  ungroup() 

plot_list <- list()

for(gene in grep("Taar6", unique(Taar_normalized_count_scale_for_plot$gene_name), value = TRUE, invert = TRUE)) {
  
  plot_df <- Taar_normalized_count_scale_for_plot %>% filter(gene_name == gene)
  comp_df <- comparisons_all_exclude_Taar6 %>% filter(gene_name == gene)
  
  p <- ggplot(plot_df, aes(x = genotype, y = normalized_count)) +
    stat_summary(
      aes(fill = genotype),
      fun = mean,
      geom = "bar",
      color = "black",
      linewidth = 0.4,
      width = 0.8,
      alpha = 0.85
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
    geom_point(
      position = position_jitter(width = 0.3),
      shape = 21,
      fill = "white",
      color = "black",
      size = 4.74,
      stroke = 0.5
    ) +
    stat_pvalue_manual(
      comp_df,
      label = "p.adj",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0,
      bracket.size = 0.5,
      size = 4
    ) +
    scale_fill_manual(values = c(
      "lsl_Tbr1_tdTomato_Het" = "#66C2A5", 
      "lsl_Tbr1_tdTomato_Het_Goofy_cre" = "#1B9E77",
      "lsl_Tbr1_tdTomato_Homo" =  "#8DA0CB",
      "lsl_Tbr1_tdTomato_Homo_Goofy_cre" = "#7570B3"
    )) +
    scale_y_continuous(
      limits = c(0, 16),
      expand = expansion(mult = c(0, 0))  # 稍微留点空
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = gene,
      x = "Genotype",
      y = "Normalized Counts"
    )
  
  plot_list[[gene]] <- p
}

# 显示一个示例
print(plot_list[["Taar7b"]])

# 导出PDF
pdf("Plots/Taars_normalized_count_barplots_exclude_Taar6.pdf", width = 6.5, height = 5)
for(gene in names(plot_list)){
  print(plot_list[[gene]])
}
dev.off()

plot_list <- list()
for(gene in "Taar6") {
  plot_df <- Taar_normalized_count_scale_for_plot %>% filter(gene_name == gene)
  comp_df <- comparisons_all_include_Taar6 %>% filter(gene_name == gene)
  
  p <- ggplot(plot_df, aes(x = genotype, y = normalized_count)) +
    stat_summary(
      aes(fill = genotype),
      fun = mean,
      geom = "bar",
      color = "black",
      linewidth = 0.4,
      width = 0.8,
      alpha = 0.85
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
    geom_point(
      position = position_jitter(width = 0.3),
      shape = 21,
      fill = "white",
      color = "black",
      size = 4.74,
      stroke = 0.5
    ) +
    stat_pvalue_manual(
      comp_df,
      label = "p.adj",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0,
      bracket.size = 0.5,
      size = 4
    ) +
    scale_fill_manual(values = c(
      "lsl_Tbr1_tdTomato_Het" = "#66C2A5", 
      "lsl_Tbr1_tdTomato_Het_Goofy_cre" = "#1B9E77",
      "lsl_Tbr1_tdTomato_Homo" =  "#8DA0CB",
      "lsl_Tbr1_tdTomato_Homo_Goofy_cre" = "#7570B3"
    )) +
    scale_y_continuous(
      limits = c(0, 26),
      expand = expansion(mult = c(0, 0))  # 稍微留点空
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = gene,
      x = "Genotype",
      y = "Normalized Counts"
    )
  
  plot_list[[gene]] <- p
}

pdf("Plots/Taars_normalized_count_barplots_Taar6.pdf", width = 6.5, height = 5)
for(gene in names(plot_list)){
  print(plot_list[[gene]])
}
dev.off()

##### Taar normalized count barplot (ylim = 12 exclude Taar6) ####
normalized_count <- read.csv("./Results/DEseq2_normalized_count.csv", header = T, row.names = 1)
Taar_normalized_count <- normalized_count[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count)) & normalized_count$gene_type %in% "protein_coding"), 
                                          grepl("lsl", colnames(normalized_count))]
res_lrt <- read.csv("Results/DEseq2_results_lrt.csv", header = T, row.names = 1)
res_lrt_functional_Taar <- res_lrt[which(grepl(pattern = "^Taar[2-9]", rownames(res_lrt)) & res_lrt$gene_type %in% "protein_coding"),]

res_OverExpression_Het <- read.csv("Results/DEseq2_results_OverExpression_Het.csv", header = T, row.names = 1)
res_Het_functional_Taar <- res_OverExpression_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Het)) & res_OverExpression_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo <- read.csv("Results/DEseq2_results_OverExpression_Homo.csv", header = T, row.names = 1)
res_Homo_functional_Taar <- res_OverExpression_Homo[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo)) & res_OverExpression_Homo$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_Het <- read.csv("Results/DEseq2_results_OverExpression_Homo_Het.csv", header = T, row.names = 1)
res_Homo_Het_functional_Taar <- res_OverExpression_Homo_Het[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_Het)) & res_OverExpression_Homo_Het$gene_type %in% "protein_coding"),]

res_OverExpression_Homo_vs_Het_control <- read.csv("Results/DEseq2_results_OverExpression_Homo_vs_Het_control.csv", header = T, row.names = 1)
res_Homo_vs_Het_control_functional_Taar <- res_OverExpression_Homo_vs_Het_control[which(grepl(pattern = "^Taar[2-9]", rownames(res_OverExpression_Homo_vs_Het_control)) & res_OverExpression_Homo_vs_Het_control$gene_type %in% "protein_coding"),]

Taar_normalized_count_scale_for_plot <- as.data.frame(Taar_normalized_count/rowMeans(Taar_normalized_count[, grep("lsl_Tbr1_tdTomato_Het_\\d",colnames(Taar_normalized_count))]))  %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols =  colnames(Taar_normalized_count)[grep("lsl", colnames(Taar_normalized_count))],
                names_to = 'sampletype',
                values_to = 'normalized_count') %>% 
  mutate(genotype = sub("_[0-9]+$", "", sampletype)) %>%  
  left_join(res_lrt_functional_Taar[grepl("gene_name|padj", colnames(res_lrt_functional_Taar))], by = "gene_name")  %>%
  mutate(Het_padj=res_Het_functional_Taar$padj[match(gene_name, rownames(res_Het_functional_Taar))],
         Homo_padj=res_Homo_functional_Taar$padj[match(gene_name, rownames(res_Homo_functional_Taar))],
         Homo_Het_padj=res_Homo_Het_functional_Taar$padj[match(gene_name, rownames(res_Homo_Het_functional_Taar))],
         Homo_vs_Het_control_padj=res_Homo_vs_Het_control_functional_Taar$padj[match(gene_name, rownames(res_Homo_vs_Het_control_functional_Taar))]) %>% 
  mutate(genotype = factor(genotype, level=c("lsl_Tbr1_tdTomato_Het","lsl_Tbr1_tdTomato_Het_Goofy_cre", "lsl_Tbr1_tdTomato_Homo", "lsl_Tbr1_tdTomato_Homo_Goofy_cre")),
         gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))) 


plot_list <- list()

for(gene in grep("Taar6", unique(Taar_normalized_count_scale_for_plot$gene_name), value = TRUE, invert = TRUE)) {
  
  plot_df <- Taar_normalized_count_scale_for_plot %>% filter(gene_name == gene)
  
  p <- ggplot(plot_df, aes(x = genotype, y = normalized_count)) +
    stat_summary(
      aes(fill = genotype),
      fun = mean,
      geom = "bar",
      color = "black",
      linewidth = 0.4,
      width = 0.8,
      alpha = 0.85
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
    geom_point(
      position = position_jitter(width = 0.3),
      shape = 21,
      fill = "white",
      color = "black",
      size = 4.74,
      stroke = 0.5
    ) +
    scale_fill_manual(values = c(
      "lsl_Tbr1_tdTomato_Het" = "#66C2A5", 
      "lsl_Tbr1_tdTomato_Het_Goofy_cre" = "#1B9E77",
      "lsl_Tbr1_tdTomato_Homo" =  "#8DA0CB",
      "lsl_Tbr1_tdTomato_Homo_Goofy_cre" = "#7570B3"
    )) +
    scale_y_continuous(
      limits = c(0, 12), expand = c(0,0)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = gene,
      x = "Genotype",
      y = "Normalized Counts"
    )
  
  plot_list[[gene]] <- p
}

# 显示一个示例
print(plot_list[["Taar7b"]])

# 导出PDF
pdf("Plots/Taars_normalized_count_barplots_exclude_Taar6_ylim_12.pdf", width = 6.5, height = 5)
for(gene in names(plot_list)){
  print(plot_list[[gene]])
}
dev.off()

##### Taar normalized count barplot (ylim=5) ("Taar2", "Taar3", "Taar4", "Taar5", "Taar7b", "Taar7e", "Taar8a", "Taar9") ####
plot_list <- list()

for(gene in c("Taar2", "Taar3", "Taar4", "Taar5", "Taar7b", "Taar7e", "Taar8a", "Taar9")) {
  
  plot_df <- Taar_normalized_count_scale_for_plot %>% filter(gene_name == gene)
  
  p <- ggplot(plot_df, aes(x = genotype, y = normalized_count)) +
    stat_summary(
      aes(fill = genotype),
      fun = mean,
      geom = "bar",
      color = "black",
      linewidth = 0.4,
      width = 0.8,
      alpha = 0.85
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
    geom_point(
      position = position_jitter(width = 0.3),
      shape = 21,
      fill = "white",
      color = "black",
      size = 4.74,
      stroke = 0.5
    ) +
    scale_fill_manual(values = c(
      "lsl_Tbr1_tdTomato_Het" = "#66C2A5", 
      "lsl_Tbr1_tdTomato_Het_Goofy_cre" = "#1B9E77",
      "lsl_Tbr1_tdTomato_Homo" =  "#8DA0CB",
      "lsl_Tbr1_tdTomato_Homo_Goofy_cre" = "#7570B3"
    )) +
    scale_y_continuous(
      limits = c(0, 5), expand = c(0,0)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = gene,
      x = "Genotype",
      y = "Normalized Counts"
    )
  
  plot_list[[gene]] <- p
}

# 显示一个示例
print(plot_list[["Taar2"]])

# 导出PDF
pdf("Plots/Taars_normalized_count_barplots_ylim_5.pdf", width = 6.5, height = 5)
for(gene in names(plot_list)){
  print(plot_list[[gene]])
}
dev.off()
