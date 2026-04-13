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

## set the workpath
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/8_Tbr1_rescue_Bulk_RNA_seq/20250724_Tbr1_rescue_Bulk_RNA_seq")
setwd(workpath)

## create a result folder
folder <- "Results/"
if (!file.exists(folder)) { dir.create(folder) }

folder_plot <- "Plots/"
if (!file.exists(folder_plot)) { dir.create(folder_plot) }

#### data preparation ####
## Load data
dat <- read.table("featureCounts/featureCounts_geneId.txt", header=T, row.names=1, as.is=T)
# as.is=T means keep all character vectors, do not convert to factors

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

## export the clean data
write.csv(good_dat_matched_gene_name_with_gene_length,"Results/good_dat_matched_gene_name_with_gene_length.csv")

###### sample info #####
sample_name <- colnames(dat)[colnames(dat) != "Length"]

meta <- data.frame(
  sample_name = sample_name,
  genotype = stringr::str_remove(sample_name, "_[1-3]$")
)

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
dds_dat$genotype <- relevel(dds_dat$genotype, "WT_Rsc")  

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
  scale_colour_manual(values =c(TKO = "orangered",
                                TKO_Rsc = "springgreen",
                                WT_Rsc = "deeppink2"))+
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
dds_lrt <- DESeq(dds_dat, test = "LRT", reduced = ~ 1)

##### export lrt results #####
## pvalue & padj
# https://www.biostars.org/p/415023/
# the p-value in DESeq2 is calculated using the wald test. 
# The null hypothesis of the wald test is that: for each gene, there is no differential expression across two sample groups (e.g., treated vs control). 
# If the p-value is small (e.g., p<0.05), the null hypothesis is rejected, as there is only 5% of chance that the null hypothesis is true. 
# However, when you have many genes being tested, by chance (5%), there is a number of genes that are not significantly expressed, but obtained significant p-values. 
# Therefore, we need to correct this problem caused by multiple testing. DESeq2 adjust the p value from wald test using Benjamini and Hochberg method (BH-adjusted p values), which is presented in the column of padj in the results object.

res_lrt <- results(dds_lrt, independentFiltering = FALSE) %>% 
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

##### export wald results #####
## export the DESeq2 results
###### res_TKO_Rsc_vs_TKO #####
res_TKO_Rsc_vs_TKO <- results(dds_dat, contrast = c("genotype", "TKO_Rsc", "TKO"), test="Wald", independentFiltering = FALSE) %>% 
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
hist(res_TKO_Rsc_vs_TKO$pvalue)
write.csv(res_TKO_Rsc_vs_TKO, "Results/DEseq2_results_TKO_Rsc_vs_TKO.csv")

###### res_TKO_Rsc_vs_WT #####
res_TKO_Rsc_vs_WT_Rsc <- results(dds_dat, contrast = c("genotype", "TKO_Rsc", "WT_Rsc"), test="Wald", independentFiltering = FALSE) %>% 
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
hist(res_TKO_Rsc_vs_WT_Rsc$pvalue)
write.csv(res_TKO_Rsc_vs_WT_Rsc, "Results/DEseq2_results_TKO_Rsc_vs_WT_Rsc.csv")

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

#### plot: Taar normalized expression counts bar plot #####
normalized_count <- read.csv("./Results/DEseq2_normalized_count.csv", header = T, row.names = 1)
res_TKO_Rsc_vs_WT_Rsc <- read.csv("Results/DEseq2_results_TKO_Rsc_vs_WT_Rsc.csv", header = T, row.names = 1)
res_TKO_Rsc_vs_TKO <- read.csv("Results/DEseq2_results_TKO_Rsc_vs_TKO.csv", header = T, row.names = 1)
res_lrt <- read.csv("Results/DEseq2_results_lrt.csv", header = T, row.names = 1)

Taar_normalized_count <- normalized_count %>%
  dplyr::filter(grepl("^Taar[2-9]", rownames(.)) & gene_type == "protein_coding")

Taar_normalized_count_for_plot <- Taar_normalized_count %>%
  select(grep("TKO|Rsc", colnames(.))) %>%
  rownames_to_column("gene_name") %>% 
  pivot_longer(
    cols = setdiff(grep("TKO|Rsc", colnames(.), value = T), "gene_type"),
    names_to = "sampletype",
    values_to = "normalized_count"
  ) %>%
  mutate(genotype = case_when(
    grepl("TKO_Rsc", sampletype) ~ "TKO_Rsc",
    grepl("WT_Rsc", sampletype) ~ "WT_Rsc",
    TRUE ~ "TKO"
  )) %>%
  left_join(
    res_lrt %>%
      as.data.frame() %>%
      dplyr::select(gene_name, padj),
    by = "gene_name"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT_Rsc", "TKO", "TKO_Rsc")),
    gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))
  )

pdf("./Plots/bar_plot_Taar_normalized_count_padj.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_for_plot,  sampletype == "TKO_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_for_plot$normalized_count))*1.05))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("grey40","springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts") +
  ggtitle("padj")
dev.off()

Taar_normalized_count_scale_for_plot <- Taar_normalized_count %>%
  # 1. 选出 TKO 或 Rsc 样本列
  dplyr::select(dplyr::matches("TKO|Rsc")) %>%
  # 2. 每行除以 WT_Rsc 样本的均值
  `/`(rowMeans(dplyr::select(Taar_normalized_count, dplyr::matches("WT_Rsc")))) %>% 
  rownames_to_column("gene_name") %>% 
  pivot_longer(
    cols = setdiff(grep("TKO|Rsc", colnames(.), value = T), "gene_type"),
    names_to = "sampletype",
    values_to = "normalized_count"
  ) %>%
  mutate(genotype = case_when(
    grepl("TKO_Rsc", sampletype) ~ "TKO_Rsc",
    grepl("WT_Rsc", sampletype) ~ "WT_Rsc",
    TRUE ~ "TKO"
  )) %>%
  left_join(
    res_lrt %>%
      as.data.frame() %>%
      dplyr::select(gene_name, padj),
    by = "gene_name"
  ) %>%
  mutate(
    padj_TKO_Rsc_vs_TKO = res_TKO_Rsc_vs_TKO$padj[match(gene_name, res_TKO_Rsc_vs_TKO$gene_name)],
    padj_TKO_Rsc_vs_WT_Rsc = res_TKO_Rsc_vs_WT_Rsc$padj[match(gene_name, res_TKO_Rsc_vs_WT_Rsc$gene_name)],
    genotype = factor(genotype, levels = c("WT_Rsc", "TKO", "TKO_Rsc")),
    gene_name = factor(gene_name, levels = rownames(Taar_normalized_count))
  ) %>% 
  mutate(
    star_TKO_Rsc_vs_TKO = case_when(
      padj_TKO_Rsc_vs_TKO < 0.001 ~ "***",
      padj_TKO_Rsc_vs_TKO < 0.01  ~ "**",
      padj_TKO_Rsc_vs_TKO < 0.05  ~ "*",
      TRUE                   ~ "n.s."
    ),
    star_TKO_Rsc_vs_WT_Rsc = case_when(
      padj_TKO_Rsc_vs_WT_Rsc < 0.001 ~ "***",
      padj_TKO_Rsc_vs_WT_Rsc < 0.01  ~ "**",
      padj_TKO_Rsc_vs_WT_Rsc < 0.05  ~ "*",
      TRUE                    ~ "n.s."
    ),
    star_all = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE                   ~ "n.s."
    )
  )

pdf("./Plots/bar_plot_Taar_normalized_count_padj_scale.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_scale_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(Taar_normalized_count_scale_for_plot$normalized_count))*1.03), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_scale_for_plot$normalized_count))*1.05))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("grey40","springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts") +
  ggtitle("padj")
dev.off()

pdf("./Plots/bar_plot_Taar_normalized_count_padj_scale_pairwise_test.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_scale_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_all, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.4), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_TKO_Rsc_vs_WT_Rsc, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.3), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_TKO_Rsc_vs_TKO, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.2), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_count_scale_for_plot$normalized_count))*1.2))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("grey40","springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts") +
  ggtitle("padj") +
  ggtitle("all\nTKO_Rsc_vs_WT_Rsc\nTKO_Rsc_vs_TKO")

dev.off()

pdf("./Plots/bar_plot_Taar_normalized_count_padj_scale_pairwise_test_errorbar.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_count_scale_for_plot, aes(gene_name, normalized_count, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1.67, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_all, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.4), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_TKO_Rsc_vs_WT_Rsc, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.3), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_count_scale_for_plot,  sampletype == "TKO_1"), 
            aes(label = star_TKO_Rsc_vs_TKO, 
                y = max(Taar_normalized_count_scale_for_plot$normalized_count)*1.2), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(-0.05,ceiling(max(Taar_normalized_count_scale_for_plot$normalized_count))*1.2))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("grey40","springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts") +
  ggtitle("padj") +
  ggtitle("all\nTKO_Rsc_vs_WT_Rsc\nTKO_Rsc_vs_TKO")

dev.off()
#### plot: receptor log2FC ####
##### TKO_Rsc_vs_TKO #####
res_TKO_Rsc_vs_TKO <- read.csv("Results/DEseq2_results_TKO_Rsc_vs_TKO.csv", header = T, row.names = 1)
receptor_res_TKO_Rsc_vs_TKO <- res_TKO_Rsc_vs_TKO[grepl(pattern = "Olfr|Taar", rownames(res_TKO_Rsc_vs_TKO)) & grepl("protein_coding", res_TKO_Rsc_vs_TKO$gene_type), ]
range(receptor_res_TKO_Rsc_vs_TKO$log2FoldChange)
# -3.705742  7.468644
receptor_res_TKO_Rsc_vs_TKO <- receptor_res_TKO_Rsc_vs_TKO %>% 
  mutate(genomic_start_relative = 1:nrow(receptor_res_TKO_Rsc_vs_TKO),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(receptor_res_TKO_Rsc_vs_TKO))

# Generate genome location data
genome_location <- as.data.frame(table(receptor_res_TKO_Rsc_vs_TKO$chr_name))
genome_location <- genome_location[match(unique(receptor_res_TKO_Rsc_vs_TKO$chr_name), genome_location$Var1), ]
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

pdf("./Plots/receptors_log2FC_with_genome_location_TKO_Rsc_vs_TKO.pdf", width = 8, height = 7)
ggplot(data = receptor_res_TKO_Rsc_vs_TKO, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  # OR red
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_TKO, grepl("^Olfr", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_TKO, grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  # Taar OR blue
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_TKO, grepl("^Taar", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_TKO, grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") + 
  
  geom_text_repel(data = subset(receptor_res_TKO_Rsc_vs_TKO, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin = 9.5, ymax = 10, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(-7.5, 10)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(receptor_res_TKO_Rsc_vs_TKO$genomic_end_relative)), 
                     breaks = seq(0, max(receptor_res_TKO_Rsc_vs_TKO$genomic_end_relative), by = 200), 
                     labels = seq(0, max(receptor_res_TKO_Rsc_vs_TKO$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes") +
  ggtitle("TKO_Rsc_vs_TKO")
dev.off()

##### TKO_Rsc_vs_WT_Rsc #####
res_TKO_Rsc_vs_WT_Rsc <- read.csv("Results/DEseq2_results_TKO_Rsc_vs_WT_Rsc.csv", header = T, row.names = 1)

receptor_res_TKO_Rsc_vs_WT_Rsc <- res_TKO_Rsc_vs_WT_Rsc[grepl(pattern = "Olfr|Taar", rownames(res_TKO_Rsc_vs_WT_Rsc)) & grepl("protein_coding", res_TKO_Rsc_vs_WT_Rsc$gene_type), ]
range(receptor_res_TKO_Rsc_vs_WT_Rsc$log2FoldChange)
# -6.104453  1.803270
receptor_res_TKO_Rsc_vs_WT_Rsc <- receptor_res_TKO_Rsc_vs_WT_Rsc %>% 
  mutate(genomic_start_relative = 1:nrow(receptor_res_TKO_Rsc_vs_WT_Rsc),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(receptor_res_TKO_Rsc_vs_WT_Rsc))

# Generate genome location data
genome_location <- as.data.frame(table(receptor_res_TKO_Rsc_vs_WT_Rsc$chr_name))
genome_location <- genome_location[match(unique(receptor_res_TKO_Rsc_vs_WT_Rsc$chr_name), genome_location$Var1), ]
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

pdf("./Plots/receptors_log2FC_with_genome_location_TKO_Rsc_vs_WT_Rsc.pdf", width = 8, height = 7)
ggplot(data = receptor_res_TKO_Rsc_vs_WT_Rsc, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  # OR red
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_WT_Rsc, grepl("^Olfr", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_WT_Rsc, grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  # Taar OR blue
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_WT_Rsc, grepl("^Taar", gene_name) & !(padj < 0.05 & abs(log2FoldChange) > 0.585)), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(receptor_res_TKO_Rsc_vs_WT_Rsc, grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) > 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") + 
  
  geom_text_repel(data = subset(receptor_res_TKO_Rsc_vs_WT_Rsc, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin = 9.5, ymax = 10, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(-7.5, 10)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(receptor_res_TKO_Rsc_vs_WT_Rsc$genomic_end_relative)), breaks = seq(0, max(receptor_res_TKO_Rsc_vs_WT_Rsc$genomic_end_relative), by = 200), labels = seq(0, max(receptor_res_TKO_Rsc_vs_WT_Rsc$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes") +
  ggtitle("TKO_Rsc_vs_WT_Rsc")
dev.off()

#### plot: volcano plot ####
## check the range of log2FC and padj
res_TKO_Rsc_vs_WT_Rsc <- read.csv("Results/DEseq2_results_TKO_Rsc_vs_WT_Rsc.csv", header = T, row.names = 1) %>% 
  as.data.frame()

##### plot: volcano plot (TKO_Rsc_vs_WT_Rsc padj) #####
range(res_TKO_Rsc_vs_WT_Rsc$log2FoldChange, na.rm = T) 
# [1] -6.104453  5.195857
range(-log10(res_TKO_Rsc_vs_WT_Rsc$padj), na.rm = T)
# 0.000000 9.588024

# check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
res_TKO_Rsc_vs_WT_Rsc[which(res_TKO_Rsc_vs_WT_Rsc$padj == 0), ] 
res_TKO_Rsc_vs_WT_Rsc_functional_OR <- res_TKO_Rsc_vs_WT_Rsc[which(grepl("Olfr", rownames(res_TKO_Rsc_vs_WT_Rsc)) & res_TKO_Rsc_vs_WT_Rsc$gene_type %in% "protein_coding"), ]
res_TKO_Rsc_vs_WT_Rsc_functional_OR_diff <- res_TKO_Rsc_vs_WT_Rsc_functional_OR[which(res_TKO_Rsc_vs_WT_Rsc_functional_OR$padj < 0.05  & abs(res_TKO_Rsc_vs_WT_Rsc_functional_OR$log2FoldChange) > 0.585) ,]
dim(res_TKO_Rsc_vs_WT_Rsc_functional_OR_diff)
# 10 15

res_TKO_Rsc_vs_WT_Rsc_functional_Taar<- res_TKO_Rsc_vs_WT_Rsc[which(grepl("Taar[2-9]", rownames(res_TKO_Rsc_vs_WT_Rsc)) & res_TKO_Rsc_vs_WT_Rsc$gene_type %in% "protein_coding"), ]
res_TKO_Rsc_vs_WT_Rsc_functional_Taar_diff <- res_TKO_Rsc_vs_WT_Rsc_functional_Taar[which(res_TKO_Rsc_vs_WT_Rsc_functional_Taar$padj < 0.05  & abs(res_TKO_Rsc_vs_WT_Rsc_functional_Taar$log2FoldChange) > 0.585) ,]
dim(res_TKO_Rsc_vs_WT_Rsc_functional_Taar_diff)
# 6 14

res_TKO_Rsc_vs_WT_Rsc_rm_functional_OR_diff_functional_Taar_diff <- subset(res_TKO_Rsc_vs_WT_Rsc, !(gene_name %in% c(res_TKO_Rsc_vs_WT_Rsc_functional_OR_diff$gene_name, res_TKO_Rsc_vs_WT_Rsc_functional_Taar_diff$gene_name)))
dim(res_TKO_Rsc_vs_WT_Rsc_rm_functional_OR_diff_functional_Taar_diff)
# 21272    14

res_TKO_Rsc_vs_WT_Rsc_volcano_plot <- res_TKO_Rsc_vs_WT_Rsc %>%
  mutate(group = case_when(
    gene_name %in% res_TKO_Rsc_vs_WT_Rsc_functional_OR_diff$gene_name ~ "DE ORs",
    gene_name %in% res_TKO_Rsc_vs_WT_Rsc_functional_Taar_diff$gene_name ~ "DE Taars",
    padj < 0.05 & abs(log2FoldChange) > 0.585 ~ "Significant DEGs",
    TRUE ~ "Other"
  ))

pdf("./Plots/volcano_plot_TKO_Rsc_vs_WT_Rsc_DEGs_padj.pdf", height = 8, width = 8)
ggplot(res_TKO_Rsc_vs_WT_Rsc_volcano_plot, aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point(data = subset(res_TKO_Rsc_vs_WT_Rsc_volcano_plot, group == c("Significant DEGs", "Other")), 
             aes(color = group), alpha = 0.75, pch = 16, size = 1) +
  geom_point(data = subset(res_TKO_Rsc_vs_WT_Rsc_volcano_plot, group == c("DE Taars", "DE ORs")), 
             aes(color = group), alpha = 0.75, pch = 16, size = 2) +
  geom_text_repel(data = subset(res_TKO_Rsc_vs_WT_Rsc_volcano_plot, group == c("DE Taars", "DE ORs")), 
                  aes(label = gene_name, color = group), size = 3) +
  scale_color_manual(
    values = c(
      "Other" = "grey",
      "Significant DEGs" = "black",
      "DE ORs" = "blue",
      "DE Taars" = "red"
    )
  ) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  xlab("log2 fold change") +
  ylab("-log10 padj") +
  xlim(-10, 10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
  ggtitle("res_TKO_Rsc_vs_WT_Rsc padj")
dev.off()

##### plot: volcano plot (TKO_Rsc_vs_TKO padj) #####
range(res_TKO_Rsc_vs_TKO$log2FoldChange, na.rm = T) 
# [1] -8.576248  7.468644
range(-log10(res_TKO_Rsc_vs_TKO$padj), na.rm = T)
# 0.000 9.234

# check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
res_TKO_Rsc_vs_TKO[which(res_TKO_Rsc_vs_TKO$padj == 0), ] 
res_TKO_Rsc_vs_TKO_functional_OR <- res_TKO_Rsc_vs_TKO[which(grepl("Olfr", rownames(res_TKO_Rsc_vs_TKO)) & res_TKO_Rsc_vs_TKO$gene_type %in% "protein_coding"), ]
res_TKO_Rsc_vs_TKO_functional_Taar<- res_TKO_Rsc_vs_TKO[which(grepl("Taar[2-9]", rownames(res_TKO_Rsc_vs_TKO)) & res_TKO_Rsc_vs_TKO$gene_type %in% "protein_coding"), ]
res_TKO_Rsc_vs_TKO_functional_OR_diff <- res_TKO_Rsc_vs_TKO_functional_OR[which(res_TKO_Rsc_vs_TKO_functional_OR$padj < 0.05  & abs(res_TKO_Rsc_vs_TKO_functional_OR$log2FoldChange) > 0.585) ,]
dim(res_TKO_Rsc_vs_TKO_functional_OR_diff)
# 0 15
res_TKO_Rsc_vs_TKO_functional_Taar_diff <- res_TKO_Rsc_vs_TKO_functional_Taar[which(res_TKO_Rsc_vs_TKO_functional_Taar$padj < 0.05  & abs(res_TKO_Rsc_vs_TKO_functional_Taar$log2FoldChange) > 0.585) ,]
dim(res_TKO_Rsc_vs_TKO_functional_Taar_diff)
# 5 14
res_TKO_Rsc_vs_TKO_rm_functional_OR_diff_functional_Taar_diff <- subset(res_TKO_Rsc_vs_TKO, !(gene_name %in% c(res_TKO_Rsc_vs_TKO_functional_OR_diff$gene_name, res_TKO_Rsc_vs_TKO_functional_Taar_diff$gene_name)))
dim(res_TKO_Rsc_vs_WT_Rsc_rm_functional_OR_diff_functional_Taar_diff)
# 21272    14

res_TKO_Rsc_vs_TKO_volcano_plot <- res_TKO_Rsc_vs_TKO %>%
  mutate(group = case_when(
    gene_name %in% res_TKO_Rsc_vs_TKO_functional_OR_diff$gene_name ~ "DE ORs",
    gene_name %in% res_TKO_Rsc_vs_TKO_functional_Taar_diff$gene_name ~ "DE Taars",
    padj < 0.05 & abs(log2FoldChange) > 0.585 ~ "Significant DEGs",
    TRUE ~ "Other"
  ))

pdf("./Plots/volcano_plot_TKO_Rsc_vs_TKO_DEGs_padj.pdf", height = 8, width = 8)
ggplot(res_TKO_Rsc_vs_TKO_volcano_plot, aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point(data = subset(res_TKO_Rsc_vs_TKO_volcano_plot, group == c("Significant DEGs", "Other")), 
             aes(color = group), alpha = 0.75, pch = 16, size = 1) +
  geom_point(data = subset(res_TKO_Rsc_vs_TKO_volcano_plot, group == c("DE Taars", "DE ORs")), 
             aes(color = group), alpha = 0.75, pch = 16, size = 2) +
  geom_text_repel(data = subset(res_TKO_Rsc_vs_TKO_volcano_plot, group == c("DE Taars", "DE ORs")), 
                  aes(label = gene_name, color = group), size = 3) +
  scale_color_manual(
    values = c(
      "Other" = "grey",
      "Significant DEGs" = "black",
      "DE ORs" = "blue",
      "DE Taars" = "red"
    )
  ) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  xlab("log2 fold change") +
  ylab("-log10 padj") +
  xlim(-10, 10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
  ggtitle("res_TKO_Rsc_vs_TKO padj")
dev.off()
