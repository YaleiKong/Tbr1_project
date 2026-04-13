Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list=ls())
getwd()  

## set seed for repeat
set.seed(24)

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
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/3_Tbr1_KO_Cre_ires_GFP_FACS_bulk_RNA_Seq/20230709_Tbr1_KO_GFP_FACS_Bulk_RNA_Seq_in_Smart_Seq_final_results/20240926_DE_analysis_results_featureCounts_multioverlap_final_results")
setwd(workpath)

## create a result folder
folder_1 <- "Results/"
if (!file.exists(folder_1)) { dir.create(folder_1) }

folder_2 <- "Plots/"
if (!file.exists(folder_2)) { dir.create(folder_2) }

#### featureCounts preparation ####
##### featureCounts ##### 
featureCounts <-  read.table("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/20230709_Tbr1_ko_Cre_ires_EGFP_FACS/20230709_Tbr1_KO_GFP_FACS_Bulk_RNA_Seq_in_Smart_Seq_final_results/featureCounts/featureCounts_with_header_fractionOfOverlapping_geneId.txt", header = T, row.names = 1)
meta <- read.table("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/20230709_Tbr1_ko_Cre_ires_EGFP_FACS/20230709_Tbr1_KO_GFP_FACS_Bulk_RNA_Seq_in_Smart_Seq_final_results/20230709_meta.txt", header = T)
colnames(featureCounts) <- ifelse(colnames(featureCounts) %in% meta$name, 
                                           meta$sample_name_new[match(colnames(featureCounts), meta$name)],
                                           colnames(featureCounts))
featureCounts_cre_egfp <- read.table("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/20230709_Tbr1_ko_Cre_ires_EGFP_FACS/20230709_Tbr1_KO_GFP_FACS_Bulk_RNA_Seq_in_Smart_Seq_final_results/featureCounts/featureCounts_Cre_Egfp_with_header.txt", header = T, row.names = 1)
featureCounts_cre_egfp <- featureCounts_cre_egfp[,1:13]
colnames(featureCounts_cre_egfp) <- colnames(featureCounts)
featureCounts <- rbind(featureCounts, featureCounts_cre_egfp)                                          

##### ensembl_gene_info ####
ensembl_gene_info <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) 
# obtain protein-coding genes and pseudogenes including "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene" and "polymorphic_pseudogene" that OR/TAAR pseudogenes belong. Check https://www.gencodegenes.org/pages/biotypes.html
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

featureCounts_protein_coding_and_pseudogene <- featureCounts[ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl,]
match_indices <- match(rownames(featureCounts_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl)
featureCounts_protein_coding_and_pseudogene_with_gene_info <- featureCounts_protein_coding_and_pseudogene %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices])
rownames(featureCounts_protein_coding_and_pseudogene_with_gene_info) <- featureCounts_protein_coding_and_pseudogene_with_gene_info$gene_name
write.csv(featureCounts_protein_coding_and_pseudogene_with_gene_info, "featureCounts/featureCounts_protein_coding_and_pseudogene_with_gene_info.csv")

#### calculate TPM ####
gene_length_kb <-  featureCounts_protein_coding_and_pseudogene_with_gene_info$Length / 1000
rpk_protein_coding_and_pseudogene_filter <- featureCounts_protein_coding_and_pseudogene_filter / gene_length_kb 
tpm_protein_coding_and_pseudogene_filter <- t(t(rpk_protein_coding_and_pseudogene_filter)/colSums(rpk_protein_coding_and_pseudogene_filter) * 1000000) %>% as.data.frame()
write.csv(tpm_protein_coding_and_pseudogene_filter, "Results/tpm_protein_coding_and_pseudogene_filter.csv")
# Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(tpm_protein_coding_and_pseudogene_filter), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

tpm_raw_count_protein_coding_and_pseudogene <- tpm_protein_coding_and_pseudogene_filter %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) 
write.csv(tpm_raw_count_protein_coding_and_pseudogene, "Results/tpm_raw_count_protein_coding_and_pseudogene.csv")

#### DESeq2 ####
featureCounts <- read.csv("featureCounts/featureCounts_protein_coding_and_pseudogene_filter.csv", header = T, row.names = 1) 

###load metadata which is the sample categories
meta <- meta %>% 
  mutate(sample_name = sample_name_new,
         genotype = ifelse(grepl("Het", sample_name), "Het", "Homo"),
         celltype = ifelse(grepl("GFP_P", sample_name), "GFP_P", "GFP_N"),
         sample_type = paste0(genotype, "_", celltype)) %>% 
  select(-c("name", "sample_name_new"))
row.names(meta) <- meta$sample_name

### create a DESeqDataSet object
featureCounts <- featureCounts %>% 
  select(meta$sample_name)
### Check that sample names match in both files (dat and meta)
# "all" function to check whether all the logic are true
all(row.names(meta) == colnames(featureCounts))

featureCounts[] <- lapply(featureCounts, function(x) as.integer(round(x)))
dds_dat <- DESeqDataSetFromMatrix(countData = featureCounts, colData = meta, design = ~ sample_type)
# design= ~ batch + condition: The design indicates how to model the samples, here, that we want to measure the effect of the condition, controlling for batch differences.
dds_dat$sample_type <- relevel(dds_dat$sample_type, "Het_GFP_N")

##### PCA #####
vsd <- vst(dds_dat, blind = FALSE)
DESeq2::plotPCA(vsd, intgroup = "sample_type")

pca_data <- DESeq2::plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)
# ggplot2
percentVar <- round(100 * attr(pca_data, "percentVar"))
plot_pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(Het_GFP_N = "black",
                                Het_GFP_P = "orangered",
                                Homo_GFP_N = "springgreen",
                                Homo_GFP_P = "steelblue"))+
  geom_text_repel(aes(label = pca_data$name), color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 15),
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title= element_blank(),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm'),
        panel.grid = element_blank())
ggsave("Plots/DESeq2_PCA.pdf", plot_pca, width = 15, height = 10)

##### DEseq2 #####
dds_DE <- DESeq(dds_dat)

##### export the DESeq2 results (Homo GFP_P/GFP_N) #####
res_Homo_GFP_P_N <- results(dds_DE, contrast = c("sample_type", "Homo_GFP_P", "Homo_GFP_N"), test = "Wald", independentFiltering = F) %>% 
  as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(res_Homo_GFP_P_N), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_Homo_GFP_P_N_protein_coding_and_pseudogene <- res_Homo_GFP_P_N %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(res_Homo_GFP_P_N_protein_coding_and_pseudogene, "Results/res_Homo_GFP_P_N_protein_coding_and_pseudogene.csv")

dim(res_Homo_GFP_P_N_protein_coding_and_pseudogene[res_Homo_GFP_P_N_protein_coding_and_pseudogene$padj < 0.05 & abs(res_Homo_GFP_P_N_protein_coding_and_pseudogene$log2FoldChange) > 0.585,])
# [1] 2021   11
##### export the DESeq2 results (Het GFP_P/GFP_N) #####
res_Het_GFP_P_N <- results(dds_DE, contrast = c("sample_type", "Het_GFP_P", "Het_GFP_N"), test = "Wald", independentFiltering = F) %>% 
  as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(res_Het_GFP_P_N), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_Het_GFP_P_N_protein_coding_and_pseudogene <- res_Het_GFP_P_N %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(res_Het_GFP_P_N_protein_coding_and_pseudogene, "Results/res_Het_GFP_P_N_protein_coding_and_pseudogene.csv")

dim(res_Het_GFP_P_N_protein_coding_and_pseudogene[res_Het_GFP_P_N_protein_coding_and_pseudogene$padj < 0.05 & abs(res_Het_GFP_P_N_protein_coding_and_pseudogene$log2FoldChange) > 0.585,])
# [1] 2920   11

##### export the DESeq2 results (GFP_P Homo/Het) #####
res_GFP_P_Homo_Het <- results(dds_DE, contrast = c("sample_type", "Homo_GFP_P", "Het_GFP_P"), test = "Wald", independentFiltering = F) %>% 
  as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(res_GFP_P_Homo_Het), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_GFP_P_Homo_Het_protein_coding_and_pseudogene <- res_GFP_P_Homo_Het %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()

write.csv(res_GFP_P_Homo_Het_protein_coding_and_pseudogene, "Results/res_GFP_P_Homo_Het_protein_coding_and_pseudogene.csv")

dim(res_GFP_P_Homo_Het_protein_coding_and_pseudogene[res_GFP_P_Homo_Het_protein_coding_and_pseudogene$padj < 0.05 & abs(res_GFP_P_Homo_Het_protein_coding_and_pseudogene$log2FoldChange) > 0.585,])
# [1] 1795   11

##### export the DESeq2 results (GFP_N Homo/Het) #####
res_GFP_N_Homo_Het <- results(dds_DE, contrast = c("sample_type", "Homo_GFP_N", "Het_GFP_N"), test = "Wald", independentFiltering = F) %>% 
  as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(res_GFP_N_Homo_Het), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_GFP_N_Homo_Het_protein_coding_and_pseudogene <- res_GFP_N_Homo_Het %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()

write.csv(res_GFP_N_Homo_Het_protein_coding_and_pseudogene, "Results/res_GFP_N_Homo_Het_protein_coding_and_pseudogene.csv")

dim(res_GFP_N_Homo_Het_protein_coding_and_pseudogene[res_GFP_N_Homo_Het_protein_coding_and_pseudogene$padj < 0.05 & abs(res_GFP_N_Homo_Het_protein_coding_and_pseudogene$log2FoldChange) > 0.585,])
# [1] 11 11
##### export the normalized counts #####
normalized_count <- DESeq2::counts(dds_DE, normalized = TRUE) %>% as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(normalized_count), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

normalized_count_protein_coding_and_pseudogene <- normalized_count %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(normalized_count_protein_coding_and_pseudogene, "Results/DEseq2_normalized_count_protein_coding_and_pseudogene.csv")

#### plots ####
normalized_count_protein_coding_and_pseudogene <- read.csv("Results/DEseq2_normalized_count_protein_coding_and_pseudogene.csv", header = T, row.names = 1)
res_Het_GFP_P_N_protein_coding_and_pseudogene <- read.csv("Results/res_Het_GFP_P_N_protein_coding_and_pseudogene.csv", header = T, row.names = 1)
res_Homo_GFP_P_N_protein_coding_and_pseudogene <- read.csv("Results/res_Homo_GFP_P_N_protein_coding_and_pseudogene.csv", header = T, row.names = 1)
res_GFP_P_Homo_Het_protein_coding_and_pseudogene <- read.csv("Results/res_GFP_P_Homo_Het_protein_coding_and_pseudogene.csv", header = T, row.names = 1)
res_GFP_N_Homo_Het_protein_coding_and_pseudogene <- read.csv("Results/res_GFP_N_Homo_Het_protein_coding_and_pseudogene.csv", header = T, row.names = 1)

##### plot 1: heatmap (Taar + diff OR)  #####
TAAR_normalized_count <- normalized_count_protein_coding_and_pseudogene %>%
  filter(grepl("Taar", gene_name) & gene_type == "protein_coding") %>%
  mutate(Het_GFP_P = rowMeans(.[,grepl("Het.*GFP_P", colnames(.))])) %>% 
  arrange(desc(Het_GFP_P)) %>% 
  select(-Het_GFP_P)

# res_GFP_P_Homo_vs_Het upregulate
res_GFP_P_Homo_Het_OR_up <- rownames(res_GFP_P_Homo_Het_protein_coding_and_pseudogene)[which(grepl("Olfr", rownames(res_GFP_P_Homo_Het_protein_coding_and_pseudogene)) &
                                                                                               grepl("protein_coding", res_GFP_P_Homo_Het_protein_coding_and_pseudogene$gene_type) &
                                                                                               res_GFP_P_Homo_Het_protein_coding_and_pseudogene$padj < 0.05 &
                                                                                               res_GFP_P_Homo_Het_protein_coding_and_pseudogene$log2FoldChange > 0.585)]

res_GFP_P_Homo_Het_up_OR_normalized_count <- normalized_count_protein_coding_and_pseudogene[res_GFP_P_Homo_Het_OR_up,] %>%
  mutate(Homo_GFP_P = rowMeans(.[,grepl("Homo.*GFP_P", colnames(.))])) %>% 
  arrange(desc(Homo_GFP_P)) %>% 
  select(-Homo_GFP_P)

# res_GFP_P_Homo_vs_Het downregulate
res_GFP_P_Homo_Het_OR_down <- rownames(res_GFP_P_Homo_Het_protein_coding_and_pseudogene)[which(grepl("Olfr", rownames(res_GFP_P_Homo_Het_protein_coding_and_pseudogene)) &
                                                                                                 grepl("protein_coding", res_GFP_P_Homo_Het_protein_coding_and_pseudogene$gene_type) &
                                                                                                 res_GFP_P_Homo_Het_protein_coding_and_pseudogene$padj < 0.05 &
                                                                                                 res_GFP_P_Homo_Het_protein_coding_and_pseudogene$log2FoldChange < -0.585)]
res_GFP_P_Homo_Het_down_OR_normalized_count <- normalized_count_protein_coding_and_pseudogene[res_GFP_P_Homo_Het_OR_down,] %>%
  mutate(Het_GFP_P = rowMeans(.[,grepl("Het.*GFP_P", colnames(.))])) %>% 
  arrange(desc(Het_GFP_P)) %>% 
  select(-Het_GFP_P)

colnames_order <- c(paste0("Tbr1_KO_Het_", 1:3, "_GFP_N"), paste0("Tbr1_KO_Homo_", 1:3, "_GFP_N"), paste0("Tbr1_KO_Het_", 1:3, "_GFP_P"), paste0("Tbr1_KO_Homo_", 1:3, "_GFP_P"))
TAAR_and_diff_OR_normalized_count_for_heatmap <- rbind(TAAR_normalized_count, res_GFP_P_Homo_Het_up_OR_normalized_count, res_GFP_P_Homo_Het_down_OR_normalized_count) 
write.csv(TAAR_and_diff_OR_normalized_count_for_heatmap, file = "Results/TAAR_and_diff_OR_normalized_count_for_heatmap.csv")

TAAR_and_diff_OR_normalized_count_for_heatmap <- TAAR_and_diff_OR_normalized_count_for_heatmap[,colnames_order] %>% 
  filter(rowSums(.)!=0)

annotation_col <- data.frame(
  sample_type = factor(c(rep("Het_GFP_N", 3), rep("Homo_GFP_N", 3), rep("Het_GFP_P", 3), rep("Homo_GFP_P", 3)))
)
rownames(annotation_col) <- colnames(TAAR_and_diff_OR_normalized_count_for_heatmap) 

TAAR <- intersect(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap), rownames(TAAR_normalized_count))
length(TAAR)
# 14
upregulated_OR <-intersect(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap), rownames(res_GFP_P_Homo_Het_up_OR_normalized_count))
length(upregulated_OR)
# 221
downregulated_OR <- intersect(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap), rownames(res_GFP_P_Homo_Het_down_OR_normalized_count))
length(downregulated_OR)
# 30

annotation_row <- data.frame(
  gene_group = factor(c(rep("TAAR", length(TAAR)), rep("upregulated_OR", length(upregulated_OR)), rep("downregulated_OR", length(downregulated_OR))))
)
rownames(annotation_row) <- rownames(TAAR_and_diff_OR_normalized_count_for_heatmap)

annotation_colors <- list(
  sample_type = c(Het_GFP_N = "gray60", Homo_GFP_N = "gray40", Het_GFP_P = "springgreen3", Homo_GFP_P = "springgreen4"),
  gene_group = c(TAAR = "deeppink3", upregulated_OR = "steelblue3", downregulated_OR = "royalblue4")
)

pdf("Plots/TAAR_and_diff_OR_normalized_count_scale_heatmap.pdf", width = 8, height = 6)
pheatmap(TAAR_and_diff_OR_normalized_count_for_heatmap,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F,
         show_colnames = F,
         scale = "row", #对行标准化,这一步很重要
         border= F, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "scaled normalized counts",
         heatmap_legend_param = list(title = NULL),
         fontsize_col = 8,
         angle_col = "0")
dev.off()

pdf("Plots/TAAR_and_diff_OR_normalized_count_heatmap.pdf", width = 8, height = 6)
pheatmap(TAAR_and_diff_OR_normalized_count_for_heatmap,
         top_annotation = HeatmapAnnotation(df = annotation_col, col = annotation_colors),
         right_annotation = rowAnnotation(df = annotation_row, col = annotation_colors),
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F,
         show_colnames = F,
         breaks = seq(0, 1000,length.out = 101),
         border= F, 
         color = colorRampPalette(c("gray100", "red3"))(100), 
         main = "normalized counts",
         heatmap_legend_param = list(title = NULL),
         fontsize_col = 8,
         angle_col = "0")
dev.off()

##### plot 2: heatmap (marker) ######
apoptosis_markers <- c("Bcl2l11","Bbc3", "Pmaip1","Hrk", "Bmf", "Bax", "Bak1", "Casp3", "Casp7", "Casp9")
progenitor_markers <- c("Neurog1","Neurod1","Neurod2","Top2a","Mki67")
imOSN_markers <- c("Gap43","Gnas","Gng8","Hdac2","Stmn4")
mOSN_markers <- c("Omp","Gnal","Gng13","Adcy3","Cnga4","Cnga2","Cngb1","Slc17a6","Stoml3")

apoptosis_markers_normalized_count <- normalized_count_protein_coding_and_pseudogene[apoptosis_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
progenitor_markers_normalized_count <- normalized_count_protein_coding_and_pseudogene[progenitor_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
imOSN_markers_normalized_count <- normalized_count_protein_coding_and_pseudogene[imOSN_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
mOSN_markers_normalized_count <- normalized_count_protein_coding_and_pseudogene[mOSN_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
markers_normalized_count <- rbind(apoptosis_markers_normalized_count, progenitor_markers_normalized_count, imOSN_markers_normalized_count, mOSN_markers_normalized_count)
write.csv(markers_normalized_count, file = "Results/markers_normalized_count.csv")
annotation_col <- data.frame(
  sample_type = c(rep("Het_GFP_N", 3), rep("Homo_GFP_N", 3), rep("Het_GFP_P", 3), rep("Homo_GFP_P", 3))
)
rownames(annotation_col) <- colnames(TAAR_and_diff_OR_normalized_count_for_heatmap) 

annotation_row <- data.frame(
  gene_group = c(rep("apoptosis", length(apoptosis_markers)), rep("progenitor", length(progenitor_markers)), 
                 rep("imOSN", length(imOSN_markers)), rep("mOSN", length(mOSN_markers))))
rownames(annotation_row) <- rownames(markers_normalized_count)

annotation_colors <- list(
  sample_type = c(Het_GFP_N = "gray60", Homo_GFP_N = "gray40", Het_GFP_P = "springgreen3", Homo_GFP_P = "springgreen4"),
  gene_group = c(apoptosis = "lightskyblue4", progenitor = "deepskyblue", imOSN = "royalblue1", mOSN = "darkblue")
)

pdf("Plots/markers_normalized_count_heatmap.pdf", width = 6, height = 6)
pheatmap(markers_normalized_count,
         top_annotation = HeatmapAnnotation(df = annotation_col, col = annotation_colors),
         right_annotation = rowAnnotation(df = annotation_row, col = annotation_colors),
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = T,
         show_colnames = F,
         border= F, 
         breaks = seq(0, 10000,length.out = 101),
         color = colorRampPalette(c("white", "red"))(100), 
         main = "normalized counts",
         heatmap_legend_param = list(title = NULL),
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = "0")
dev.off()


