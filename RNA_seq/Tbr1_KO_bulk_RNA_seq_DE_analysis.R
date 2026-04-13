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
                "circlize","ComplexHeatmap","hutils",
                "broom"
)

lapply(library_Ls,function(x){suppressPackageStartupMessages(library(x,character.only = T))})

## set the workpath
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/2_Tbr1_KO_WT_Het_Homo_bulk_RNA_Seq/202511_Tbr1_KO_WT_Het_Homo_bulk_RNA_seq")
setwd(workpath)

#### featureCounts preparation ####
## Load featureCountsa
featureCounts_with_length <- read.table("featureCounts/featureCounts_geneId.txt", header=T, row.names=1, as.is=T)
# as.is=T means keep all character vectors, do not convert to factors
featureCounts <- featureCounts_with_length[,-grep("Length", colnames(featureCounts_with_length))] 
write.csv(featureCounts,"featureCounts/featureCounts.csv")

##### ensembl_gene_info #####
ensembl_gene_info <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) 
# obtain protein-coding genes and pseudogenes including "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene" and "polymorphic_pseudogene" that OR/TAAR pseudogenes belong. 
# Check https://www.gencodegenes.org/pages/biotypes.html
ensembl_gene_info_protein_coding_and_pseudogene <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) %>% 
  mutate(gene_id_Ensembl_remove_version_meta = str_remove(gene_id_Ensembl, "\\.[0-9]")) %>% 
  filter(gene_type %in% c("protein_coding", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "polymorphic_pseudogene")) %>% 
  filter(grepl("^chr[1-9XY]", chr_name)) 
ensembl_gene_info_protein_coding_and_pseudogene_nodup <- ensembl_gene_info_protein_coding_and_pseudogene %>% 
  arrange(annotation_source) %>% 
  filter(!duplicated(gene_name)) %>% 
  arrange(chr_name, genomic_start, genomic_end)
nrow(ensembl_gene_info_protein_coding_and_pseudogene_nodup)
# 24899

##### featureCounts preparation #####
featureCounts_protein_coding_and_pseudogene_with_length <- featureCounts_with_length[ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl,]
featureCounts_protein_coding_and_pseudogene <- featureCounts_protein_coding_and_pseudogene_with_length %>% 
  select(-"Length")
match_indices <- match(rownames(featureCounts_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl)
featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info <- featureCounts_protein_coding_and_pseudogene_with_length %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices])
rownames(featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info) <- featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info$gene_name
rownames(featureCounts_protein_coding_and_pseudogene) <-  featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info$gene_name
rownames(featureCounts_protein_coding_and_pseudogene_with_length) <-  featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info$gene_name
write.csv(featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info, "featureCounts/featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info.csv")
write.csv(featureCounts_protein_coding_and_pseudogene, "featureCounts/featureCounts_protein_coding_and_pseudogene.csv")
write.csv(featureCounts_protein_coding_and_pseudogene_with_length, "featureCounts/featureCounts_protein_coding_and_pseudogene_with_length.csv")

featureCounts_protein_coding_and_pseudogene_with_length["Tbr1",]
# Length Hetero.1 Hetero.2 Hetero.3 Hetero.4 Homo.1 Homo.2 Homo.3 Homo.4 WT.1 WT.2 WT.3 WT.4
# Tbr1   5486       53       47       39      274     19     29     36     32   55   93   91   81

##### calculate TPM #####
gene_length_kb_protein_coding_and_pseudogene <- featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info$Length / 1000
rpk_protein_coding_and_pseudogene <- featureCounts_protein_coding_and_pseudogene / gene_length_kb_protein_coding_and_pseudogene 
tpm_protein_coding_and_pseudogene <- t(t(rpk_protein_coding_and_pseudogene)/colSums(rpk_protein_coding_and_pseudogene) * 1000000) %>% as.data.frame()
write.csv(tpm_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/tpm_protein_coding_and_pseudogene.csv")

tpm_protein_coding_and_pseudogene["Tbr1",]
# Hetero.1 Hetero.2  Hetero.3 Hetero.4    Homo.1    Homo.2    Homo.3    Homo.4     WT.1     WT.2     WT.3     WT.4
# Tbr1 1.170522  1.04185 0.8987476 5.495681 0.4480702 0.5861865 0.6288393 0.7165054 1.228306 1.859639 2.007506 1.773855

match_indices <- match(rownames(tpm_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)
tpm_protein_coding_and_pseudogene_with_gene_info <- tpm_protein_coding_and_pseudogene %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices])
write.csv(tpm_protein_coding_and_pseudogene_with_gene_info, "DE_analysis_with_replace/Results/tpm_protein_coding_and_pseudogene_with_gene_info.csv")

tpm_protein_coding_and_pseudogene_with_gene_info["Tbr1",]
# Hetero.1 Hetero.2  Hetero.3 Hetero.4    Homo.1    Homo.2    Homo.3    Homo.4     WT.1     WT.2     WT.3     WT.4      gene_type chr_name genomic_start genomic_end gene_name
# Tbr1 1.170522  1.04185 0.8987476 5.495681 0.4480702 0.5861865 0.6288393 0.7165054 1.228306 1.859639 2.007506 1.773855 protein_coding     chr2      61802930    61814114      Tbr1

sum(rowSums(tpm_protein_coding_and_pseudogene) == 0)
# 3246

rownames(tpm_protein_coding_and_pseudogene[which(rowSums(tpm_protein_coding_and_pseudogene) == 0),])
# Genes starting with 'Gm-' or ending with 'Rik' are annotated genes that do not have a canonical name (yet).

quantile(rowSums(tpm_protein_coding_and_pseudogene), probs = seq(0, 1, 0.1))
# 0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100% 
# 0.000000e+00 0.000000e+00 6.925367e-01 6.554404e+00 2.822209e+01 9.138411e+01 1.964891e+02 3.288585e+02 5.410305e+02 9.896759e+02 2.614428e+05

quantile(rowSums(tpm_protein_coding_and_pseudogene[grepl("Olfr", rownames(tpm_protein_coding_and_pseudogene)), 1:9]), probs = seq(0, 1, 0.1))
# 0%         10%         20%         30%         40%         50%         60%         70%         80%         90%        100% 
# 0.0000000   0.2949268   4.6478406   7.9146611  11.4448875  14.6231637  18.1644024  22.9120525  29.6268182  41.4562731 225.0224119 

#### DESeq2 ####
meta <- data.frame(sample = colnames(featureCounts_with_length)[-1],
                   sample_name = str_replace(colnames(featureCounts_with_length)[-1], "\\.", "_"),
                   genotype = str_remove(colnames(featureCounts_with_length)[-1], ".\\d"))

featureCounts <- read.csv("featureCounts/featureCounts_protein_coding_and_pseudogene.csv", header = T, row.names = 1) %>% 
  select(meta$sample) 

colnames(featureCounts) <- meta$sample_name

### load metadata which is the sample categories
meta <- meta[,c("sample_name", "genotype")]
row.names(meta) <- meta$sample_name
dim(meta)
summary(meta)

### Check that sample names match in both files (dat and meta)
# "all" function to check whether all the logic are true
which(row.names(meta) %in% colnames(featureCounts))  # return the position of the former
match(row.names(meta) , colnames(featureCounts))     # return the position of the latter
all(row.names(meta) == colnames(featureCounts))

##### DESeqDataSetFromMatrix #####
### create a DESeqDataSet object
featureCounts[] <- lapply(featureCounts, function(x) as.integer(round(x)))

dds_dat <- DESeqDataSetFromMatrix(countData = featureCounts, colData = meta, design= ~ genotype)

dds_dat$genotype <- relevel(dds_dat$genotype, "WT")  

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
  scale_colour_manual(values =c(WT = "black",
                                Hetero = "orangered",
                                Homo = "springgreen"))+
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
ggsave("DE_analysis_with_replace/Plots/PCA.pdf", p.pca, width = 15, height = 10)

##### DEseq2 #####
dds_DE <- DESeq(dds_dat, test="LRT", reduced=~1, minReplicatesForReplace = 4)
# DESeq2 默认在重复数足够多时，会对极端 Cook’s distance 的样本计数做自动替换，而不是简单删除；默认门槛是该样本所在组至少有 7 个重复。替换值来自该基因在各样本标准化计数上的 trimmed mean。

## Cook’s distance
assays(dds_DE)[["cooks"]]["Tbr1",]
# Hetero_1     Hetero_2     Hetero_3     Hetero_4       Homo_1       Homo_2       Homo_3       Homo_4         WT_1         WT_2         WT_3         WT_4 
# 0.5404101826 0.8068420252 0.9523968792 7.0725293899 0.1095951036 0.0003771794 0.0172165910 0.0561288598 0.1503102117 0.0066624989 0.0845936837 0.0004931904 

##### export the DESeq2 results (all) #####
# set **independentFiltering to FALSE** so that there will be less NA values of padj (see below for padj). 
# Especially for lowly expressed genes such as olfactory receptors. 
# The independent filtering is designed only to filter out low count genes to the extent that they are not enriched with small p-values.
# using the mean of normalized counts as a filter statistic

res_all <- results(dds_DE, independentFiltering = F) %>% 
  as.data.frame()

res_all["Tbr1", ]
# baseMean log2FoldChange     lfcSE     stat       pvalue         padj
# Tbr1 52.28479      -1.475337 0.2414129 38.25719 4.926706e-09 4.826828e-06

## Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(res_all), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_all_protein_coding_and_pseudogene <- res_all %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) %>% 
  na.omit()
write.csv(res_all_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/res_all_protein_coding_and_pseudogene.csv")

nrow(res_all_protein_coding_and_pseudogene[which(res_all_protein_coding_and_pseudogene$padj < 0.05 & abs(res_all_protein_coding_and_pseudogene$log2FoldChange) > 0.585),])
# 138

##### export the DESeq2 results (homo/wt) #####
res_homo_wt <- results(dds_DE, contrast = c("genotype", "Homo", "WT"), test="Wald", independentFiltering = F) %>% 
  as.data.frame() 

match_indices <- match(rownames(res_homo_wt), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_homo_wt_protein_coding_and_pseudogene <- res_homo_wt %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) %>% 
  na.omit()

res_homo_wt_protein_coding_and_pseudogene["Tbr1",]
# baseMean log2FoldChange     lfcSE      stat      pvalue         padj      gene_type chr_name genomic_start genomic_end gene_name
# Tbr1 52.28479      -1.475337 0.2414129 -6.111258 9.88488e-10 2.130587e-06 protein_coding     chr2      61802930    61814114      Tbr1

write.csv(res_homo_wt_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/res_homo_wt_protein_coding_and_pseudogene.csv")

nrow(res_homo_wt_protein_coding_and_pseudogene[which(res_homo_wt_protein_coding_and_pseudogene$padj < 0.05 & abs(res_homo_wt_protein_coding_and_pseudogene$log2FoldChange) > 0.585),])
# 214

##### export the DESeq2 results (het/wt) #####
res_het_wt <- results(dds_DE, contrast = c("genotype", "Hetero", "WT"), test="Wald", independentFiltering = F) %>% 
  as.data.frame() 

## Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(res_het_wt), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_het_wt_protein_coding_and_pseudogene <- res_het_wt %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) %>% 
  na.omit()

write.csv(res_het_wt_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/res_het_wt_protein_coding_and_pseudogene.csv")

nrow(res_het_wt_protein_coding_and_pseudogene[which(res_het_wt_protein_coding_and_pseudogene$padj < 0.05 & abs(res_het_wt_protein_coding_and_pseudogene$log2FoldChange) > 0.585),])
# 15

res_het_wt_protein_coding_and_pseudogene["Tbr1",]
# baseMean log2FoldChange     lfcSE      stat      pvalue      padj      gene_type chr_name genomic_start genomic_end gene_name
# Tbr1 52.28479     -0.6633797 0.2252502 -2.945079 0.003228721 0.3931743 protein_coding     chr2      61802930    61814114      Tbr1

##### export the DESeq2 results (homo/het) #####
res_homo_het <- results(dds_DE, contrast = c("genotype", "Homo", "Hetero"), test="Wald", independentFiltering = F) %>% 
  as.data.frame() 

## Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(res_homo_het), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_homo_het_protein_coding_and_pseudogene <- res_homo_het %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) %>% 
  na.omit()

write.csv(res_homo_het_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/res_homo_het_protein_coding_and_pseudogene.csv")

nrow(res_homo_het_protein_coding_and_pseudogene[which(res_homo_het_protein_coding_and_pseudogene$padj < 0.05 & abs(res_homo_het_protein_coding_and_pseudogene$log2FoldChange) > 0.585),])
# 57

res_homo_het_protein_coding_and_pseudogene["Tbr1",]
# baseMean log2FoldChange     lfcSE    stat      pvalue      padj      gene_type chr_name genomic_start genomic_end gene_name
# Tbr1 52.28479     -0.8119569 0.2495641 -3.2535 0.001139926 0.0790031 protein_coding     chr2      61802930    61814114      Tbr1

##### export normalized data #####
normalized_count_all_protein_coding_and_pseudogene <- DESeq2::counts(dds_DE, normalized = TRUE, replaced = T) %>% as.data.frame()
# replaced = T
# DESeq2 的默认规则是：对某个被判为 outlier 的基因-样本点，用该基因在 所有样本 上的 normalized counts 的 trimmed mean 作为预测中心，再按该样本自己的 size factor 或 normalization factor 调整回去。它不是只看同组样本，也不是简单删掉一个点后求均值。

normalized_count_all_protein_coding_and_pseudogene["Tbr1",]
# Hetero_1 Hetero_2 Hetero_3 Hetero_4  Homo_1   Homo_2   Homo_3   Homo_4    WT_1    WT_2     WT_3     WT_4
# Tbr1  56.1826 47.14707 42.09642 53.72467 21.2969 27.87469 30.77341 33.09999 60.2898 82.5703 92.51478 79.84683

match_indices <- match(rownames(normalized_count_all_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

normalized_count_all_protein_coding_and_pseudogene <- normalized_count_all_protein_coding_and_pseudogene %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) %>% 
  na.omit()

write.csv(normalized_count_all_protein_coding_and_pseudogene, "DE_analysis_with_replace/Results/normalized_count_all_protein_coding_and_pseudogene.csv")

#### plot ####
##### check results #####
res_all_protein_coding_and_pseudogene <- read.csv("DE_analysis_with_replace/Results/res_all_protein_coding_and_pseudogene.csv", header = T, row.names = 1) %>% 
  as.data.frame()
res_homo_wt_protein_coding_and_pseudogene <- read.csv("DE_analysis_with_replace/Results/res_homo_wt_protein_coding_and_pseudogene.csv", header = T, row.names = 1) %>% 
  as.data.frame()

res_homo_wt_protein_coding_and_pseudogene_with_res_all <- res_homo_wt_protein_coding_and_pseudogene %>% 
  mutate(padj_all = res_all_protein_coding_and_pseudogene$padj[match(gene_name, res_all_protein_coding_and_pseudogene$gene_name)])

### check the range of log2FC and padj
range(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange, na.rm = T) 
# [1] -7.766931  4.107565

range(-log10(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj), na.rm = T) 
# [1] 0.00000 22.37475

range(-log10(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all), na.rm = T) 
# [1] 0 Inf

## check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
## see which genes have padj = 0 
res_homo_wt_protein_coding_and_pseudogene_with_res_all[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj == 0), ]
# no 

res_homo_wt_protein_coding_and_pseudogene_with_res_all[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all == 0), ]
# baseMean log2FoldChange     lfcSE       stat    pvalue      padj      gene_type chr_name genomic_start genomic_end gene_name padj_all
# Kdm5d   512.5984      0.3108343 0.9746660  0.3189137 0.7497920 0.9414001 protein_coding     chrY        897788      956786     Kdm5d        0
# Eif2s3y 893.9314     -0.4482201 0.4430784 -1.0116044 0.3117273 0.6662589 protein_coding     chrY       1010543     1028847   Eif2s3y        0
# Uty     298.9154     -1.0135944 0.8951838 -1.1322752 0.2575188 0.6078663 protein_coding     chrY       1096861     1245759       Uty        0
# Ddx3y   926.8799     -0.4744180 0.4972797 -0.9540264 0.3400703 0.6928048 protein_coding     chrY       1260771     1286629     Ddx3y        0

min(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all[res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all != 0])
#  4.210346e-28

### add padj_all_limit limit = 10^(-30)
res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit <- res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all+10^(-30)

range(-log10(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit), na.rm = T) 
# 0 30

write.csv(res_homo_wt_protein_coding_and_pseudogene_with_res_all, file = 'DE_analysis_with_replace/Results/res_homo_wt_protein_coding_and_pseudogene_with_res_all.csv')

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit < 0.05 & 
      abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange) >= 0.585, na.rm = T) 
# 127

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange >= 0.585, na.rm = T) 
# 59

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange <= -0.585, na.rm = T) 
# 68

rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all)[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj < 0.05 & 
                                                                         res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit < 0.05 & 
                                                                         res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange >= 0.585)]
# [1] "Khdc1a"   "Gm28438"  "Gm28661"  "Ccdc150"  "Daw1"     "Fcamr"    "Fcmr"     "Gm2000"   "Ifi213"   "Ifi209"   "Ifi204"   "Ifi211"   "Rsad2"    "Gpr141"   "Agtr1a"   "Gm10044"  "Saysd1"   "Ptpn20"  
# [19] "Ldb3"     "Phf11b"   "Olfm4"    "Tymp"     "Rtp4"     "Stfa3"    "Timmdc1"  "Mx2"      "Mrpl14"   "Gm4951"   "Npas4"    "Ms4a6c"   "Ms4a6b"   "Slc18a2"  "Zbp1"     "Gm14410"  "Zfp931"   "Cpa3"    
# [37] "Riiad1"   "BC021767" "Fcgr1"    "Ifi44"    "Rmdn1"    "Zmynd12"  "Isg15"    "Cfap299"  "Plac8"    "Oasl2"    "Oas2"     "Oas1b"    "Oas1a"    "Gm42421"  "Wdr95"    "Kcne3"    "Ubqlnl"   "Trim30a" 
# [55] "Irf7"     "Ddx60"    "Bst2"     "Jaml"     "Plscr2"  

rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all)[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj < 0.05 & 
                                                                         res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj_all_limit < 0.05 & 
                                                                         res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange <= (-0.585))]
# [1] "Cdk5r2"        "Efhd1"         "Olfr1411"      "Pou2f1"        "Taar2"         "Taar3"         "Taar4"         "Taar5"         "Taar6"         "Taar7a"        "Taar7b"        "Taar7d"       
# [13] "Taar7e"        "Taar7f"        "Taar8b"        "Taar8c"        "Taar9"         "Cdh23"         "Hcn2"          "Diras1"        "Tbc1d30"       "Fgf18"         "Tenm2"         "Olfr393"      
# [25] "Rtn4rl1"       "Tcam1"         "Klhl29"        "Map3k9"        "Scgn"          "Dcdc2a"        "Fam169a"       "Rgs7bp"        "Gprin2"        "Gucy1b2"       "Siah3"         "Rtp1"         
# [37] "Dscam"         "Umodl1"        "Kctd16"        "Gnal"          "Tmem151a"      "Olfr1432"      "Gpr158"        "Tbr1"          "Olfr1196"      "Lrrc4c"        "Pak7"          "Rbm12"        
# [49] "Soga1"         "Kcns1"         "Kcnc4"         "Gabbr2"        "Zpld2"         "Fbxo2"         "Hipk2"         "Tmem178b"      "Olfr683"       "Olfr513"       "Olfr532"       "Cend1"        
# [61] "Muc6"          "Slc6a2"        "Jph3"          "Ldlr"          "Pik3cb"        "Obp1a"         "Frmpd3"        "A730046J19Rik"

res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR <- res_homo_wt_protein_coding_and_pseudogene_with_res_all[grepl(pattern = "^Olfr", rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all)), ]
res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR <- res_homo_wt_protein_coding_and_pseudogene_with_res_all[grepl(pattern = "^Taar[2-9]", rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all)), ] 

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 & 
      abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange) >= 0.585, na.rm = T) 
# 7

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange >= 0.585, na.rm = T) 
# 0

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange <= -0.585, na.rm = T) 
# 7

## upregulated OR
rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR)[which((res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
                                                                             res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 &
                                                                             res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange >= 0.585))]

## downregulated OR
rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR)[which((res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
                                                                             res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 &
                                                                             res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange <= -0.585))]
# [1] "Olfr1411" "Olfr393"  "Olfr1432" "Olfr1196" "Olfr683"  "Olfr513"  "Olfr532"

## differentially expreseed OR
Diff_OR <- rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR)[which((res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
                                                                                        res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 & 
                                                                                        res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange <= -0.585))]
# "Olfr1411" "Olfr393"  "Olfr1432" "Olfr1196" "Olfr683"  "Olfr513"  "Olfr532" 

OR_info <- read.csv("/Users/yaleikong/bioinformatics_analysis/database/olfactory_receptor_information/Olfactory_receptor_info_arranged_by_QianLi/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv",
                    row.names = 1, header = T)

Diff_OR %in% rownames(OR_info)[which(OR_info$OR_Class == "Class_II" & OR_info$gene_type == "protein_coding")]
# [1]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
Diff_OR %in% rownames(OR_info)[which(OR_info$OR_Class == "Class_I" & OR_info$gene_type == "protein_coding")]
# [1] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
# 1个Class I OR 6个Class II OR 均为functional receptor 

sum(OR_info$OR_Class == "Class_II")
# 1250

sum(OR_info$OR_Class == "Class_II" & OR_info$gene_type == "protein_coding" )
# 1006

sum(OR_info$OR_Class == "Class_I")
# 158

sum(OR_info$OR_Class == "Class_I" & OR_info$gene_type == "protein_coding" )
# 131

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj_all_limit < 0.05 & 
      abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$log2FoldChange) >= 0.585, na.rm = T) 
# 13

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$log2FoldChange >= 0.585, na.rm = T) 
# 0

sum(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj_all_limit < 0.05 & 
      res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$log2FoldChange <= -0.585, na.rm = T) 
# 13

rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR)[which((res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj < 0.05 & 
                                                                               res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj_all_limit < 0.05 & 
                                                                               abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$log2FoldChange) >= 0.585))]
# [1] "Taar2"  "Taar3"  "Taar4"  "Taar5"  "Taar6"  "Taar7a" "Taar7b" "Taar7d" "Taar7e" "Taar7f" "Taar8b" "Taar8c" "Taar9" 

##### plot 1: volcano plot (Homo vs. WT) #####
## OR_diff
res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff <- res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj < 0.05 & 
                                                                                                                                               res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$padj_all_limit < 0.05 & 
                                                                                                                                               abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$log2FoldChange) >= 0.585 & 
                                                                                                                                               res_homo_wt_protein_coding_and_pseudogene_with_res_all_OR$gene_type %in% "protein_coding"),]
dim(res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff)
# [1]  7 13

## TAAR_diff
res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff <- res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR[which(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj < 0.05 & 
                                                                                                                                                   res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$padj_all_limit < 0.05 & 
                                                                                                                                                   abs(res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$log2FoldChange) >= 0.585 & 
                                                                                                                                                   res_homo_wt_protein_coding_and_pseudogene_with_res_all_TAAR$gene_type %in% "protein_coding"),]
dim(res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff)
# [1] 13 13

## rm_functional_OR_diff_functional_TAAR_diff_Gucy1b2
res_homo_wt_protein_coding_and_pseudogene_with_res_all_rm_functional_OR_diff_functional_TAAR_diff_Gucy1b2_Tbr1 <- subset(res_homo_wt_protein_coding_and_pseudogene_with_res_all, !(gene_name %in% c(res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff$gene_name, 
                                                                                                                                                                                               res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff$gene_name, 
                                                                                                                                                                                               "Gucy1b2", "Tbr1")))
dim(res_homo_wt_protein_coding_and_pseudogene_with_res_all_rm_functional_OR_diff_functional_TAAR_diff_Gucy1b2)
# [1] 21533    13

pdf("DE_analysis_with_replace/Plots/volcano_plot_Homo_vs_WT_DEGs.pdf", height = 4, width = 4)
ggplot(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_rm_functional_OR_diff_functional_TAAR_diff_Gucy1b2_Tbr1, 
       aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point(aes(colour = padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), alpha = 0.75, pch = 16, size = 1) + 
  scale_color_manual(values = c("grey", "black")) +
  # The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff, alpha = 0.75, pch = 16, size = 2, color = "blue") + 
  # geom_text(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff, aes(label = gene_name)) + 
  
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff, alpha = 0.75, pch = 16, size = 2, color = "red") + 
  # geom_text(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff, aes(label = gene_name)) + 
  
  ## Gucy1b2
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Gucy1b2",], alpha = 0.75, pch = 16, size = 1, color = "green") + 
  # geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Gucy1b2",], aes(label = gene_name), color = "green", size = 2) + 
  
  ## Tbr1
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1",], alpha = 0.75, pch = 16, size = 1, color = "black") + 
  # geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1",], aes(label = gene_name), color = "black", size = 2) + 
  
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  xlab("log2 fold change") + 
  ylab("-log10 padj") + 
  xlim(floor(min(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange)), ceiling(max(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,  ceiling(max(-log10(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj)))))+
  ggtitle("OR: blue, TAAR: red")
dev.off()

pdf("DE_analysis_with_replace/Plots/volcano_plot_Homo_vs_WT_DEGs_with_DE_receptor_name.pdf", height = 4, width = 4)
ggplot(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_rm_functional_OR_diff_functional_TAAR_diff_Gucy1b2_Tbr1, 
       aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point(aes(colour = padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), alpha = 0.75, pch = 16, size = 0.5) + 
  scale_color_manual(values = c("grey", "black")) +
  # The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  
  ## OR
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff, alpha = 0.75, pch = 16, size = 1, color = "blue") + 
  geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_OR_diff, aes(label = gene_name), color = "blue", size = 2, max.overlaps = 20) + 
  
  ## Taar
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff, alpha = 0.75, pch = 16, size = 1, color = "red") + 
  geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all_functional_TAAR_diff, aes(label = gene_name), color = "red", size = 2) + 
  
  ## Gucy1b2
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Gucy1b2",], alpha = 0.75, pch = 16, size = 1, color = "green") + 
  geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Gucy1b2",], aes(label = gene_name), color = "green", size = 2) + 
  
  ## Tbr1
  geom_point(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1",], alpha = 0.75, pch = 16, size = 1, color = "black") + 
  geom_text_repel(data = res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1",], aes(label = gene_name), color = "black", size = 2) + 
  
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  xlab("log2 fold change") + 
  ylab("-log10 padj") + 
  xlim(floor(min(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange)), ceiling(max(res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,  ceiling(max(-log10(res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj)))))+
  ggtitle("OR: blue, TAAR: red")
dev.off()

##### plot 2: Taar normalized expression counts bar plot #####
res_homo_wt_protein_coding_and_pseudogene_with_res_all <- read.csv(file = 'DE_analysis_with_replace/Results/res_homo_wt_protein_coding_and_pseudogene_with_res_all.csv',
                                                                   header = T, row.names = 1)
res_het_wt_protein_coding_and_pseudogene <- read.csv(file = 'DE_analysis_with_replace/Results/res_het_wt_protein_coding_and_pseudogene.csv',
                                                     header = T, row.names = 1)
res_homo_het_protein_coding_and_pseudogene <- read.csv(file = 'DE_analysis_with_replace/Results/res_homo_het_protein_coding_and_pseudogene.csv',
                                                       header = T, row.names = 1)

normalized_count_all_protein_coding_and_pseudogene <- read.csv("DE_analysis_with_replace/Results/normalized_count_all_protein_coding_and_pseudogene.csv", header = T, row.names = 1)

Taar_normalized_expression <- normalized_count_all_protein_coding_and_pseudogene[which(grepl(pattern = "^Taar[2-9]", rownames(normalized_count_all_protein_coding_and_pseudogene)) & 
                                                                                         normalized_count_all_protein_coding_and_pseudogene$gene_type %in% "protein_coding"), 
                                                                                 grepl("Het|WT|Homo", colnames(normalized_count_all_protein_coding_and_pseudogene))]

Taar_normalized_expression_for_plot <- Taar_normalized_expression %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(Taar_normalized_expression), 
                names_to = 'sample_name', 
                values_to = 'normalized_counts') %>% 
  left_join(res_homo_wt_protein_coding_and_pseudogene_with_res_all[,c("gene_name","padj_all_limit")], by = "gene_name")  %>%
  mutate(padj_homo_vs_wt = res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj[match(gene_name, res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_name)],
         padj_homo_vs_het = res_homo_het_protein_coding_and_pseudogene$padj[match(gene_name, res_homo_het_protein_coding_and_pseudogene$gene_name)],
         padj_het_vs_wt = res_het_wt_protein_coding_and_pseudogene$padj[match(gene_name, res_het_wt_protein_coding_and_pseudogene$gene_name)]) %>% 
  mutate(genotype = str_remove(sample_name, "_\\d")) %>% 
  mutate(genotype = factor(genotype, level = c("WT", "Hetero", "Homo"))) 

pdf("DE_analysis_with_replace/Plots/bar_plot_taar_normalized_counts.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_expression_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_expression_for_plot,  sample_name == "WT_1"), 
            aes(label = sprintf("%.3f", padj_all_limit), 
                y = ceiling(max(Taar_normalized_expression_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_expression_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图内填充的颜色
  
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

Taar_normalized_expression_compared_to_control_for_plot <- as.data.frame(Taar_normalized_expression/rowMeans(Taar_normalized_expression[, grep("WT",colnames(Taar_normalized_expression))]))  %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(Taar_normalized_expression), 
                names_to = 'sample_name', 
                values_to = 'normalized_counts') %>% 
  left_join(res_homo_wt_protein_coding_and_pseudogene_with_res_all[,c("gene_name","padj_all_limit")], by = "gene_name")  %>%
  mutate(padj_homo_vs_wt = res_homo_wt_protein_coding_and_pseudogene_with_res_all$padj[match(gene_name, res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_name)],
         padj_homo_vs_het = res_homo_het_protein_coding_and_pseudogene$padj[match(gene_name, res_homo_het_protein_coding_and_pseudogene$gene_name)],
         padj_het_vs_wt = res_het_wt_protein_coding_and_pseudogene$padj[match(gene_name, res_het_wt_protein_coding_and_pseudogene$gene_name)]) %>% 
  mutate(genotype = str_remove(sample_name, "_\\d")) %>% 
  mutate(genotype = factor(genotype, level = c("WT", "Hetero", "Homo"))) %>% 
  mutate(
    star_wt_het = case_when(
      padj_het_vs_wt < 0.001 ~ "***",
      padj_het_vs_wt < 0.01  ~ "**",
      padj_het_vs_wt < 0.05  ~ "*",
      TRUE                   ~ "n.s."
    ),
    star_wt_homo = case_when(
      padj_homo_vs_wt < 0.001 ~ "***",
      padj_homo_vs_wt < 0.01  ~ "**",
      padj_homo_vs_wt < 0.05  ~ "*",
      TRUE                    ~ "n.s."
    ),
    star_het_homo = case_when(
      padj_homo_vs_het < 0.001 ~ "***",
      padj_homo_vs_het < 0.01  ~ "**",
      padj_homo_vs_het < 0.05  ~ "*",
      TRUE                     ~ "n.s."
    ),
    star_all = case_when(
      padj_all_limit < 0.001 ~ "***",
      padj_all_limit < 0.01  ~ "**",
      padj_all_limit < 0.05  ~ "*",
      TRUE                   ~ "n.s."
    )
  )


pdf("DE_analysis_with_replace/Plots/bar_plot_taar_normalized_counts_compared_to_WT.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_expression_compared_to_control_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_expression_compared_to_control_for_plot,  sample_name == "WT_1"), 
            aes(label = sprintf("%.3f", padj_all_limit), 
                y = ceiling(max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图内填充的颜色
  
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

pdf("DE_analysis_with_replace/Plots/bar_plot_taar_normalized_counts_compared_to_WT_pairwise_test.pdf", width = 12, height = 7)
ggplot(data = Taar_normalized_expression_compared_to_control_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + 
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.4, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, color = "black", fill = "white") + #散点
  
  geom_text(data = subset(Taar_normalized_expression_compared_to_control_for_plot,  sample_name == "WT_1"), 
            aes(label = star_all, 
                y = max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts)*1.5), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_expression_compared_to_control_for_plot,  sample_name == "WT_1"), 
            aes(label = star_wt_homo, 
                y = max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts)*1.4), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_expression_compared_to_control_for_plot,  sample_name == "WT_1"), 
            aes(label = star_wt_het, 
                y = max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts)*1.3), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  geom_text(data = subset(Taar_normalized_expression_compared_to_control_for_plot,  sample_name == "WT_1"), 
            aes(label = star_het_homo, 
                y = max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts)*1.2), 
            position = position_dodge(width = 0.9), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Taar_normalized_expression_compared_to_control_for_plot$normalized_counts))*1.1))+ 
  scale_fill_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("gray60", "blue4", "green4")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts") +
  ggtitle("all\nhomo_vs_wt\nhet_vs_wt\nhomo_vs_het")
dev.off()

##### plot 3: log2FC and genome location (Taar gene cluster) ####
taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all <- res_homo_wt_protein_coding_and_pseudogene_with_res_all[grepl(pattern = "^Taar|Slc18b1|Stx7|Moxd1", row.names(res_homo_wt_protein_coding_and_pseudogene_with_res_all)) & 
                                                                                                                        grepl("protein_coding", res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_type), ]

pdf("DE_analysis_without_cutoff/Plots/taar_log2FC_with_genome_location.pdf", width = 7, height = 7)
ggplot(data = taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all, aes(genomic_start, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all, (padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585) == TRUE), size = 2, pch = 16, na.rm = T, color = "blue")+ 
  geom_point(data = subset(taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all, (padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585) == F), size = 2, pch = 1, na.rm = T, color = "blue")+ 
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))-1, ceiling(max(taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))+2))+
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
  
  geom_rect(aes(xmin = genomic_start, xmax = genomic_end, ymin = 1, ymax = 1.5), colour = "black", fill = "black")+
  
  geom_text(data = data.frame(label = c("Slc18b1", "Taar", "Stx7", "Moxd1"), 
                              x = c(23815000, 23900000, 24168000, 24260000), 
                              y = rep(2, 4)), 
            mapping = aes(x = x, y = y, label = label))+
  geom_text(data = data.frame(label = gsub("Taar", "", taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_name[grepl("Taar",taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_name)]), 
                              x = ((taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_start+taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end)/2)[grepl("Taar",taar_res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_name)], 
                              y = rep(2, 14)), 
            mapping = aes(x = x, y = y, label = label), 
            size = 2
  )

dev.off()

##### plot 4: log2FC with genome location (OR and Taar gene cluster)####
olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all <- res_homo_wt_protein_coding_and_pseudogene_with_res_all[grepl(pattern = "Olfr|Taar", rownames(res_homo_wt_protein_coding_and_pseudogene_with_res_all)) & grepl("protein_coding", res_homo_wt_protein_coding_and_pseudogene_with_res_all$gene_type), ] %>% 
  arrange(factor(chr_name, levels = paste0("chr", c(1:19, "X", "Y"))))

olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all <- olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all %>% 
  mutate(genomic_start_relative = 1:nrow(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all))

write.csv(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, file = "DE_analysis_with_replace/Results/olfactory_receptors_res_homo_wt.csv")
# Generate genome location data
genome_location <- as.data.frame(table(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$chr_name))
genome_location <- genome_location[match(unique(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$chr_name), genome_location$Var1), ]
genome_location$start <- 0
genome_location <- genome_location %>%
  mutate(end = cumsum(Freq), # Create a cumulative frequency column
         start = if_else(row_number() == 1, start, start + lag(end))) 
genome_location$Var1 <- factor(genome_location$Var1, levels = paste0("chr", c(1:19, "X", "Y")))

col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)
col_for_chr <- col_for_chr_universal

# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y"))

col_for_chr <- col_for_chr[match(chr_all, unique(genome_location$Var1))]
col_for_chr <- col_for_chr[!is.na(col_for_chr)]

pdf("DE_analysis_without_cutoff/Plots/Olfactory_receptors_log2FC_with_genome_location.pdf", width = 8, height = 7)
ggplot(data = olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Olfr", gene_name)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Olfr", gene_name) & padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name)) , size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name) & padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") +
  
  geom_text_repel(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin =  ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))*0.95 , ymax =  ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))*0.98, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(floor(min(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))-1, ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange)))) +
  scale_x_continuous(expand = c(0, 0), 
                     limits = c(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative)), breaks = seq(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative), by = 200), 
                     labels = seq(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes")
dev.off()

pdf("DE_analysis_without_cutoff/Plots/Olfactory_receptors_log2FC_with_genome_location_without_gene_name.pdf", width = 8, height = 7)
ggplot(data = olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Olfr", gene_name)), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Olfr", gene_name) & padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name)) , size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name) & padj < 0.05 & padj_all < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = TRUE, color = "blue") +
  
  # geom_text_repel(data = subset(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin =  ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))*0.95 , ymax =  ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))*0.98, fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange))-1, ceiling(max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$log2FoldChange)))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative)), breaks = seq(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative), by = 200), labels = seq(0, max(olfactory_receptor_res_homo_wt_protein_coding_and_pseudogene_with_res_all$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes")
dev.off()

##### plot: Tbr1 normalized counts ####
normalized_count_all_protein_coding_and_pseudogene <- read.csv("DE_analysis_with_replace/Results/normalized_count_all_protein_coding_and_pseudogene.csv",
                                                               header = T, row.names = 1)
Tbr1_normalized_count <- normalized_count_all_protein_coding_and_pseudogene["Tbr1",grepl("Het|WT|Homo", colnames(normalized_count_all_protein_coding_and_pseudogene))]
res_homo_wt_protein_coding_and_pseudogene_with_res_all <- read.csv(file = 'DE_analysis_with_replace/Results/res_homo_wt_protein_coding_and_pseudogene_with_res_all.csv',
                                                                   header = T, row.names = 1)
res_het_wt_protein_coding_and_pseudogene <- read.csv(file = 'DE_analysis_with_replace/Results/res_het_wt_protein_coding_and_pseudogene.csv',
                                                     header = T, row.names = 1)
res_homo_het_protein_coding_and_pseudogene <- read.csv(file = 'DE_analysis_with_replace/Results/res_homo_het_protein_coding_and_pseudogene.csv',
                                                       header = T, row.names = 1)

res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1",]


## 1. 提取 Tbr1 normalized counts
tbr1_df <- data.frame(
  sample = colnames(Tbr1_normalized_count),
  normalized_count = as.numeric(Tbr1_normalized_count[1, ]),
  stringsAsFactors = FALSE
)

## 2. 样本分组
tbr1_df$group <- sub("_[0-9]+$", "", tbr1_df$sample)
tbr1_df$group <- factor(tbr1_df$group, levels = c("WT", "Hetero", "Homo"))

## 3. 计算每组均值和标准差（barplot 用）
summary_df <- tbr1_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(normalized_count),
    sd = sd(normalized_count),
    .groups = "drop"
  )

## 4. 提取 Tbr1 的 pairwise padj
padj_het_wt  <- res_het_wt_protein_coding_and_pseudogene["Tbr1", "padj"]
padj_homo_wt <- res_homo_wt_protein_coding_and_pseudogene_with_res_all["Tbr1", "padj"]
padj_homo_het <- res_homo_het_protein_coding_and_pseudogene["Tbr1", "padj"]

format_padj <- function(p) {
  if (is.na(p)) {
    "padj = NA"
  } else if (p < 0.001) {
    paste0("padj = ", format(p, scientific = TRUE, digits = 2))
  } else {
    paste0("padj = ", signif(p, 3))
  }
}

label_het_wt   <- format_padj(padj_het_wt)
label_homo_wt  <- format_padj(padj_homo_wt)
label_homo_het <- format_padj(padj_homo_het)

## 5. 设置显著性标注高度
y_max <- max(tbr1_df$normalized_count, na.rm = TRUE)

ann_df <- data.frame(
  x1 = c(1, 1, 2),
  x2 = c(2, 3, 3),
  y  = c(y_max * 1.10, y_max * 1.22, y_max * 1.34),
  label = c(label_het_wt, label_homo_wt, label_homo_het)
)

## 6. 作图
p <- ggplot(summary_df, aes(x = group, y = mean)) +
  geom_col(width = 0.6, fill = "grey80", color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  geom_point(
    data = tbr1_df,
    aes(x = group, y = normalized_count),
    inherit.aes = FALSE,
    size = 2.5,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_segment(
    data = ann_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = ann_df,
    aes(x = x1, xend = x1, y = y, yend = y - y_max * 0.02),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = ann_df,
    aes(x = x2, xend = x2, y = y, yend = y - y_max * 0.02),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ann_df,
    aes(x = (x1 + x2) / 2, y = y + y_max * 0.03, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    x = NULL,
    y = "Normalized counts",
    title = "Tbr1 normalized counts"
  ) +
  coord_cartesian(ylim = c(0, y_max * 1.48)) +
  theme_classic(base_size = 14)

p
