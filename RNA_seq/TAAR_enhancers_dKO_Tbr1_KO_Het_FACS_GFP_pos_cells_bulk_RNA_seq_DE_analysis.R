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
                "circlize","ComplexHeatmap","hutils",
                "clusterProfiler", "org.Mm.eg.db", "ReactomePA", "GOSemSim", "DOSE", "meshes", "topGO", "Rgraphviz"
                
)
# reshape2 
# makes it easy to transform data between wide and long formats.
# ggbreak
# An implementation of scale functions for setting axis breaks of a 'gg' plot.
# gtools
# to use mixedsort, check: http://veda.cs.uiuc.edu/TCGA_classify/RWR/DRaWR/suppressPackageStartupMessages(library/gtools/html/mixedsort.html
# ComplexHeatmap
# to add legend in Circos plot
# hutils
# to use Switch, which is the vectorized form of base::switch. like ifelse is vevtorized form of if else. Check: https://www.delftstack.com/howto/r/use-a-vectorized-if-function-with-multiple-conditions-in-r/
# http://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GOSemSim")
# BiocManager::install("DOSE")
# BiocManager::install("meshes")
# BiocManager::install("ReactomePA")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")

setwd("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/4_Taar_enhancer1_and_2_DKO_Tbr1_KO_Cre_ires_GFP_FACS_bulk_RNA_Seq/20230925_Taar_enhancer1_and_2_dko_Tbr1_ko_FACS_in_Smart_Seq/DE_Analysis_20241014")
getwd() 

## create a result folder
folder_1 <- "Results/"
if (!file.exists(folder_1)) { dir.create(folder_1) }

folder_2 <- "Plots/"
if (!file.exists(folder_2)) { dir.create(folder_2) }

#### featureCounts preparation ####
featureCounts_with_gene_length <- read.table(file = "featureCounts/featureCounts_with_header_fractionOfOverlapping_geneId.txt", header = T, row.names = 1, as.is=T)
featureCounts_egfp_with_gene_length <- read.table(file = "featureCounts/featureCounts_cre_with_header.txt", header = T, row.names = 1, as.is=T)
featureCounts_cre_with_gene_length <- read.table(file = "featureCounts/featureCounts_EGFP_with_header.txt", header = T, row.names = 1, as.is=T)

sample_info <- data.frame(
  sample_number = paste0("sample_", 1:9),
  sample_name = c("Tbr1_Het_GFP_dko_homo_1", "Tbr1_Het_GFP_dko_het_1", "Tbr1_Het_GFP_dko_het_2", 
                  "Tbr1_Het_GFP_dko_homo_2", "Tbr1_Het_GFP_dko_homo_3", "Tbr1_Het_GFP_dko_het_3", 
                  "Tbr1_Het_GFP_dko_homo_4", "Tbr1_Het_GFP_dko_homo_5"),
  stringsAsFactors = FALSE # Prevents automatic conversion of strings to factors
)
sample_info$genotype <- str_remove(sample_info$sample_name, "_\\d$")
rownames(sample_info) <- sample_info$sample_name

colnames(featureCounts_with_gene_length)
colnames(featureCounts_with_gene_length) <- c("gene_length", sample_info$sample_name)
write.csv(featureCounts_with_gene_length, file = "featureCounts/featureCounts_with_gene_length_with_sample_names.csv")

colnames(featureCounts_egfp_with_gene_length) <- c("gene_length", sample_info$sample_name)
colnames(featureCounts_cre_with_gene_length) <- c("gene_length", sample_info$sample_name)
featureCounts_egfp_cre_with_gene_length <- rbind(featureCounts_egfp_with_gene_length, featureCounts_cre_with_gene_length)
write.csv(featureCounts_egfp_cre_with_gene_length, file = "featureCounts/featureCounts_egfp_cre_with_gene_length_with_sample_names.csv")

### ensembl_gene preparation
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

featureCounts_with_gene_length_protein_coding_and_pseudogene <- featureCounts_with_gene_length[ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl,]
match_indices <- match(rownames(featureCounts_with_gene_length_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_id_Ensembl)
featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info <- featureCounts_with_gene_length_protein_coding_and_pseudogene %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices])
rownames(featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info) <- featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info$gene_name
write.csv(featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info, "featureCounts/featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info.csv")

#### calculate TPM ####
gene_length_kb <-  featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info$gene_length / 1000
rpk_protein_coding_and_pseudogene <- featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info[,grepl("Tbr",colnames(featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info))] / gene_length_kb 
tpm_protein_coding_and_pseudogene <- t(t(rpk_protein_coding_and_pseudogene)/colSums(rpk_protein_coding_and_pseudogene) * 1000000) %>% as.data.frame()

# Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(tpm_protein_coding_and_pseudogene), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

tpm_protein_coding_and_pseudogene <- tpm_protein_coding_and_pseudogene %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices],
                gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices]) 
write.csv(tpm_protein_coding_and_pseudogene, "Results/tpm_protein_coding_and_pseudogene.csv")

#### calculate Cre/GFP TPM ####
Cre_EGFP_gene_length_kb <- featureCounts_egfp_cre_with_gene_length$gene_length / 1000
rpk_Cre_EGFP <- featureCounts_egfp_cre_with_gene_length[,grepl("Tbr",colnames(featureCounts_egfp_cre_with_gene_length))] / Cre_EGFP_gene_length_kb 
tpm_Cre_EGFP <- t(t(rpk_Cre_EGFP)/colSums(rpk_protein_coding_and_pseudogene) * 1000000) %>% 
  as.data.frame()
write.csv(tpm_Cre_EGFP, "Results/tpm_Cre_EGFP.csv")
tpm_Cre_EGFP[c("Cre","EGFP"),]

#### DESeq2####
meta <- sample_info[,c("sample_name", "genotype")]
featureCounts <- read.csv("featureCounts/featureCounts_with_gene_length_protein_coding_and_pseudogene_with_gene_info.csv", header = T, row.names = 1) %>% 
  select(meta$sample_name)
all(row.names(meta) == colnames(featureCounts))
featureCounts[] <- lapply(featureCounts, function(x) as.integer(round(x)))
dds_dat <- DESeqDataSetFromMatrix(countData = featureCounts, colData = meta, design = ~ genotype)
# design= ~ batch + condition: The design indicates how to model the samples, here, that we want to measure the effect of the condition, controlling for batch differences.
dds_dat$genotype <- relevel(dds_dat$genotype, "Tbr1_Het_GFP_dko_het")

##### PCA #####
vsd <- vst(dds_dat, blind = FALSE)
pca_data <- DESeq2::plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
# ggplot2
percentVar <- round(100 * attr(pca_data, "percentVar"))
plot_pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(Tbr1_Het_GFP_dko_het = "springgreen",
                                Tbr1_Het_GFP_dko_homo = "steelblue"))+
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

##### export the DESeq2 results #####
res_homo_vs_het <- results(dds_DE, contrast = c("genotype", "Tbr1_Het_GFP_dko_homo", "Tbr1_Het_GFP_dko_het"), test="Wald", independentFiltering = FALSE) %>% 
  as.data.frame() 
# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(res_homo_vs_het), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_homo_vs_het <- res_homo_vs_het %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(res_homo_vs_het, "Results/res_homo_vs_het.csv")

dim(res_homo_vs_het[res_homo_vs_het$padj < 0.05 & abs(res_homo_vs_het$log2FoldChange) > 0.585,])
# [1] 531  11

##### export the normalized counts #####
normalized_count <- DESeq2::counts(dds_DE, normalized = TRUE) %>% as.data.frame()

# Precute matches to avoid multiple calls to `match`
match_indices <- match(rownames(normalized_count), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

normalized_count <- normalized_count %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(normalized_count, "Results/normalized_count.csv")

#### plot####
##### plot 1: heatmap (functional taar and diff OR)  #####
normalized_count <- read.csv("Results/normalized_count.csv", header = T, row.names = 1)
res_homo_vs_het <- read.csv("Results/res_homo_vs_het.csv", header = T, row.names = 1)

functional_TAAR_normalized_count <- normalized_count %>%
  filter(grepl("Taar[2-9]", gene_name) & gene_type %in% "protein_coding") %>%
  mutate(het = rowMeans(.[,grepl("het", colnames(.))])) %>% 
  arrange(desc(het)) %>% 
  select(-het)

up_functional_OR <- res_homo_vs_het$gene_name[which(grepl("Olfr", res_homo_vs_het$gene_name) & 
                                                      res_homo_vs_het$gene_type %in% "protein_coding" & 
                                                      res_homo_vs_het$padj < 0.05 & 
                                                      res_homo_vs_het$log2FoldChange > 0.585)]
length(up_functional_OR)
# [1] 41
up_functional_OR_normalized_count <- normalized_count[up_functional_OR,] %>%
  mutate(homo = rowMeans(.[,grepl("homo", colnames(.))])) %>% 
  arrange(desc(homo)) %>% 
  select(-homo)

down_functional_OR <- res_homo_vs_het$gene_name[which(grepl("Olfr", res_homo_vs_het$gene_name) & 
                                                        res_homo_vs_het$gene_type %in% "protein_coding" & 
                                                        res_homo_vs_het$padj < 0.05 & 
                                                        res_homo_vs_het$log2FoldChange < -0.585)]
length(down_functional_OR)
# [1] 3
down_functional_OR_normalized_count <- normalized_count[down_functional_OR,] %>%
  mutate(het = rowMeans(.[,grepl("het", colnames(.))])) %>% 
  arrange(desc(het)) %>% 
  select(-het)

TAAR_and_diff_OR_normalized_count_for_heatmap <- rbind(functional_TAAR_normalized_count, up_functional_OR_normalized_count, down_functional_OR_normalized_count) 
TAAR_and_diff_OR_normalized_count_for_heatmap <- TAAR_and_diff_OR_normalized_count_for_heatmap[,meta$sample_name] %>% 
  filter(rowSums(.)!=0) %>% 
  select(c(colnames(.)[grepl("het", colnames(.))], colnames(.)[grepl("homo", colnames(.))])) %>% 
  as.matrix()
write.csv(TAAR_and_diff_OR_normalized_count_for_heatmap, "Results/TAAR_and_diff_OR_normalized_count_for_heatmap.csv")

TAAR <- rownames(TAAR_and_diff_OR_normalized_count_for_heatmap)[which(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap) %in% rownames(functional_TAAR_normalized_count))]
length(TAAR)
# 13
upregulated_OR <- rownames(TAAR_and_diff_OR_normalized_count_for_heatmap)[which(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap) %in% rownames(up_functional_OR_normalized_count))]
length(upregulated_OR)
# 41
downregulated_OR <- rownames(TAAR_and_diff_OR_normalized_count_for_heatmap)[which(rownames(TAAR_and_diff_OR_normalized_count_for_heatmap) %in% rownames(down_functional_OR_normalized_count))]
length(downregulated_OR)
# 3

annotation_col <- data.frame(
  sample_type = factor(c(rep("Tbr1_Het_GFP_dko_het", 3), rep("Tbr1_Het_GFP_dko_homo", 5)))
)
rownames(annotation_col) <- colnames(TAAR_and_diff_OR_normalized_count_for_heatmap) 

annotation_row <- data.frame(
  gene_group = c(rep("TAAR", length(TAAR)), rep("upregulated_OR", length(upregulated_OR)), rep("downregulated_OR", length(downregulated_OR))))

rownames(annotation_row) <- rownames(TAAR_and_diff_OR_normalized_count_for_heatmap)

annotation_colors <- list(
  sample_type = c(Tbr1_Het_GFP_dko_het = "springgreen3", Tbr1_Het_GFP_dko_homo = "springgreen4"),
  gene_group = c(TAAR = "deeppink3", upregulated_OR = "steelblue3", downregulated_OR = "royalblue4")
)

pdf("Plots/TAAR_and_diff_OR_normalized_count_heatmap.pdf", width = 8, height = 6)
pheatmap(TAAR_and_diff_OR_normalized_count_for_heatmap,
         top_annotation = HeatmapAnnotation(df = annotation_col, col = annotation_colors),
         right_annotation = rowAnnotation(df = annotation_row, col = annotation_colors),
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F,
         show_colnames = F,
         breaks = seq(0, 1000, length.out = 101),
         border= F, 
         color = colorRampPalette(c("gray100", "red3"))(100), 
         main = "normalized counts",
         heatmap_legend_param = list(title = NULL),
         fontsize_col = 8,
         angle_col = "0")
dev.off()

##### plot 2: heatmap (markers) #####
colnames_order <- c(paste0("Tbr1_Het_GFP_dko_het_", c(1:2,4)), paste0("Tbr1_Het_GFP_dko_homo_", 1:5))
apoptosis_markers <- c("Bcl2l11","Bbc3", "Pmaip1","Hrk", "Bmf", "Bax", "Bak1", "Casp3", "Casp7", "Casp9")
progenitor_markers <- c("Neurog1","Neurod1","Neurod2","Top2a","Mki67")
imOSN_markers <- c("Gap43","Gnas","Gng8","Hdac2","Stmn4")
mOSN_markers <- c("Omp","Gnal","Gng13","Adcy3","Cnga4","Cnga2","Cngb1","Slc17a6","Stoml3")

apoptosis_markers_normalized_count <- normalized_count[apoptosis_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
progenitor_markers_normalized_count <- normalized_count[progenitor_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
imOSN_markers_normalized_count <- normalized_count[imOSN_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
mOSN_markers_normalized_count <- normalized_count[mOSN_markers,] %>% 
  select(colnames_order) %>% 
  arrange(desc(rowSums(.)))
markers_normalized_count <- rbind(apoptosis_markers_normalized_count, progenitor_markers_normalized_count, imOSN_markers_normalized_count, mOSN_markers_normalized_count) %>% 
  as.matrix()
write.csv(markers_normalized_count, file = "Results/markers_normalized_count.csv")

annotation_col <- data.frame(
  sample_type = factor(c(rep("Tbr1_Het_GFP_dko_het", 3), rep("Tbr1_Het_GFP_dko_homo", 5)))
)
rownames(annotation_col) <- colnames(markers_normalized_count) 

annotation_row <- data.frame(
  gene_group = c(rep("apoptosis", length(apoptosis_markers)), rep("progenitor", length(progenitor_markers)), 
                 rep("imOSN", length(imOSN_markers)), rep("mOSN", length(mOSN_markers))))
markers <- rownames(markers_normalized_count)
rownames(annotation_row) <- markers

annotation_colors <- list(
  sample_type = c(Tbr1_Het_GFP_dko_het = "springgreen3", Tbr1_Het_GFP_dko_homo = "springgreen4"),
  gene_group = c(apoptosis = "lightskyblue4", progenitor = "deepskyblue", imOSN = "royalblue1", mOSN = "darkblue")
)

pdf("Plots/markers_normalized_count_heatmap.pdf", width = 8, height = 6)
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

