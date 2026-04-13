Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list=ls())
getwd()  

## set seed for repeat
set.seed(42)

## library 
library_Ls <- c(
  "DESeq2", "LSD", "gProfileR",
  "stringr", "readxl",
  "dplyr", "tidyr", "tidyverse",
  "ggplot2", "pheatmap", "ggrepel", "gg.gap", "RColorBrewer",
  "grid", "gridExtra", "gplots", "reshape2", "ggpubr", "ggbreak", "gtools", "textshape",
  "circlize", "ComplexHeatmap", "hutils",
  "broom"
)

invisible(lapply(
  library_Ls,
  function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
))

## set the workpath
workpath <- c("/Users/yaleikong/bioinformatics_analysis/analysis_results/2_Tbr1_project/Tbr1_bulk_RNA_Seq/5_Tbr1_flox_Goofy_or_Omp_Cre_bulk_RNA_Seq/20240628_Tbr1_flox_Goofy_or_Omp_Cre_bulk_RNA_Seq/DE_Analysis_20241014")
setwd(workpath)

## create a result folder
folder_1 <- "Results/"
if (!file.exists(folder_1)) { dir.create(folder_1) }

folder_2 <- "Plots/"
if (!file.exists(folder_2)) { dir.create(folder_2) }

#### featureCounts preparation ####
## Load featureCountsa
featureCounts_with_length <- read.table("featureCounts/featureCounts_geneId.txt", header=T, row.names=1, as.is=T)
# as.is=T means keep all character vectors, do not convert to factors

Tbr1_exon_featureCounts <- read.csv("featureCounts/Tbr1_exon_featureCounts_with_length_with_sample_names.csv", header = T, row.names = 1)

###### sample info #####
meta <- read.csv("metadata.csv", header = T)
rownames(meta) <- meta$sample_name

## change the colnames of featureCounts
# 提取原始列名
original_colnames <- colnames(featureCounts_with_length)

# 创建一个映射表，将 sample 列名映射到 sample_name
colname_mapping <- setNames(meta$sample_name, meta$sample)

# 使用 ifelse 和 %in% 语句替换列名
new_colnames <- ifelse(original_colnames %in% meta$sample, colname_mapping[original_colnames], original_colnames)

# 赋值新的列名给数据框
colnames(featureCounts_with_length) <- new_colnames

write.csv(featureCounts_with_length,"featureCounts/featureCounts_with_length_with_sample_names.csv")

featureCounts <- featureCounts_with_length[,-grep("Length", colnames(featureCounts_with_length))] %>% 
  relocate(meta$sample_name)
write.csv(featureCounts,"featureCounts/featureCounts_with_sample_names.csv")

## change the colnames of featureCounts
# 提取原始列名
original_colnames <- colnames(Tbr1_exon_featureCounts)

# 创建一个映射表，将 sample 列名映射到 sample_name
colname_mapping <- setNames(meta$sample_name, meta$sample)

# 使用 ifelse 和 %in% 语句替换列名
new_colnames <- ifelse(original_colnames %in% meta$sample, colname_mapping[original_colnames], original_colnames)

colnames(Tbr1_exon_featureCounts) <- new_colnames
write.csv(Tbr1_exon_featureCounts,"featureCounts/Tbr1_exon_featureCounts_with_length_with_sample_names.csv")

Tbr1_exon_featureCounts <- Tbr1_exon_featureCounts %>% 
  as.data.frame() %>% 
  dplyr::select(-"Length") %>% 
  relocate(meta$sample_name)
write.csv(Tbr1_exon_featureCounts,"featureCounts/Tbr1_exon_featureCounts_with_sample_names.csv")

##### ensembl_gene_info ####
ensembl_gene_info <- read.table(file = "/Users/yaleikong/bioinformatics_analysis/database/reference/Mus_musculus/GRCm38/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, as.is = T, fill = T) 
# olfacotry receptor genes "protein_coding", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "polymorphic_pseudogene"
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

##### featureCounts_protein_coding_and_pseudogene #####
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

#### calculate TPM #####
gene_length_kb_protein_coding_and_pseudogene <- featureCounts_protein_coding_and_pseudogene_with_length_with_gene_info$Length / 1000
rpk_protein_coding_and_pseudogene <- featureCounts_protein_coding_and_pseudogene / gene_length_kb_protein_coding_and_pseudogene 
tpm_protein_coding_and_pseudogene <- t(t(rpk_protein_coding_and_pseudogene)/colSums(rpk_protein_coding_and_pseudogene) * 1000000) %>% as.data.frame()
write.csv(tpm_protein_coding_and_pseudogene, "Results/tpm_protein_coding_and_pseudogene.csv")

sum(rowSums(tpm_protein_coding_and_pseudogene) == 0)
# 3278
rownames(tpm_protein_coding_and_pseudogene[which(rowSums(tpm_protein_coding_and_pseudogene) == 0),])
# Genes starting with 'Gm-' or ending with 'Rik' are annotated genes that do not have a canonical name (yet).

quantile(rowSums(tpm_protein_coding_and_pseudogene), probs = seq(0, 1, 0.1))
# 0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100% 
# 0.000000e+00 0.000000e+00 9.884595e-01 8.524789e+00 4.273914e+01 1.124010e+02 2.265028e+02 3.807136e+02 6.278021e+02 1.188381e+03 2.834775e+05 

quantile(rowSums(tpm_protein_coding_and_pseudogene[grepl("Olfr", rownames(tpm_protein_coding_and_pseudogene)), 1:9]), probs = seq(0, 1, 0.1))
# 0%         10%         20%         30%         40%         50%         60%         70%         80%         90%        100% 
# 0.0000000   0.4178173   7.7781023  15.7470331  22.3599191  30.0923124  38.7258682  47.3298750  61.3056941  90.0788381 743.7469033 

tpm_protein_coding_and_pseudogene_Omp <- tpm_protein_coding_and_pseudogene[,grepl("Omp", colnames(tpm_protein_coding_and_pseudogene))]
write.csv(tpm_protein_coding_and_pseudogene_Omp,file="./Results/tpm_protein_coding_and_pseudogene_Omp.csv")

tpm_protein_coding_and_pseudogene_Goofy <- tpm_protein_coding_and_pseudogene[,!grepl("Omp", colnames(tpm_protein_coding_and_pseudogene))]
write.csv(tpm_protein_coding_and_pseudogene_Goofy,file="./Results/tpm_protein_coding_and_pseudogene_Goofy.csv")

tpm_protein_coding_and_pseudogene_Omp <- read.csv("./Results/tpm_protein_coding_and_pseudogene_Omp.csv", row.names = 1, header = T)

tpm_protein_coding_and_pseudogene_Omp["Tbr1",]
# Tbr1_flox_Omp_Cre_3 Omp_Cre_3 Tbr1_flox_Omp_Cre_4 Omp_Cre_4 Omp_Cre_1 Tbr1_flox_Omp_Cre_1 Omp_Cre_2 Tbr1_flox_Omp_Cre_2
# Tbr1           0.4389565  2.196067           0.8792267  1.123445   1.11743            1.199724 0.8750885           0.8393681

tpm_protein_coding_and_pseudogene_Goofy <- read.csv("./Results/tpm_protein_coding_and_pseudogene_Goofy.csv", row.names = 1, header = T)

tpm_protein_coding_and_pseudogene_Goofy["Tbr1",]
# Tbr1_flox_Goofy_Cre_4 Tbr1_flox_4 Tbr1_flox_Goofy_Cre_1 Tbr1_flox_Goofy_Cre_2 Tbr1_flox_1 Tbr1_flox_Goofy_Cre_3 Tbr1_flox_2 Tbr1_flox_3
# Tbr1              0.562805    1.823644             0.9235383             0.5415748   0.9388611              1.103276   0.9019063    1.655877

##### exon (Goofy) #####
Tbr1_exon_featureCounts_with_length <- read.csv("featureCounts/Tbr1_exon_featureCounts_with_length_with_sample_names.csv", header = T, row.names = 1)
Tbr1_exon_featureCounts <- read.csv("featureCounts/Tbr1_exon_featureCounts_with_sample_names.csv", header = T, row.names = 1)
Tbr1_exon_length_kb <- Tbr1_exon_featureCounts_with_length$Length / 1000
Tbr1_exon_rpk <- Tbr1_exon_featureCounts / Tbr1_exon_length_kb 
Tbr1_exon_tpm <- t(t(Tbr1_exon_rpk)/colSums(featureCounts_protein_coding_and_pseudogene) * 1000000) %>% as.data.frame()
rownames(Tbr1_exon_tpm) <- paste0("exon", 1:6)

Tbr1_exon_tpm_plot_Goofy <- Tbr1_exon_tpm[,!grepl("Omp",colnames(Tbr1_exon_tpm))]  %>% 
  tibble::rownames_to_column(var = 'exon_name') %>% 
  pivot_longer( cols = colnames(Tbr1_exon_tpm)[!grepl("Omp",colnames(Tbr1_exon_tpm))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  mutate(genotype = ifelse(grepl("Goofy", sampletype), "Tbr1_flox_Goofy_Cre", "Tbr1_flox")) %>% 
  mutate(genotype = factor(genotype, level = c("Tbr1_flox", "Tbr1_flox_Goofy_Cre"))) 

t_test_results <- Tbr1_exon_tpm_plot_Goofy %>%
  group_by(exon_name) %>%
  summarise(t_test = list(tidy(t.test(normalized_counts ~ genotype))), .groups = 'drop') %>%
  unnest(t_test)
# exon_name  estimate estimate1 estimate2 statistic p.value parameter conf.low conf.high method                  alternative
# <chr>         <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>     <dbl> <chr>                   <chr>      
#   1 exon1     -0.0706       0.269     0.339  -0.542    0.615       4.26 -0.424       0.283 Welch Two Sample t-test two.sided  
# 2 exon2      0.173        0.173     0       1.57     0.215       3    -0.178       0.524 Welch Two Sample t-test two.sided  
# 3 exon3      0.104        0.307     0.203   0.554    0.601       5.68 -0.362       0.570 Welch Two Sample t-test two.sided  
# 4 exon4      0.000561     0.480     0.479   0.00160  0.999       4.73 -0.918       0.919 Welch Two Sample t-test two.sided  
# 5 exon5     -0.0544       1.26      1.31   -0.0576   0.956       4.73 -2.53        2.42  Welch Two Sample t-test two.sided  
# 6 exon6      0.366        0.682     0.316   2.94     0.0471      3.69  0.00787     0.723 Welch Two Sample t-test two.sided  

Tbr1_exon_tpm_plot_Goofy <- Tbr1_exon_tpm_plot_Goofy %>% 
  left_join(t_test_results[,c("exon_name","p.value")], by = "exon_name") 

pdf("./Plots/bar_plot_Tbr1_exon_normalized_counts_Goofy.pdf", width = 12, height = 7)
ggplot(data = Tbr1_exon_tpm_plot_Goofy, aes(exon_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.4, position = position_dodge(0.75)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(x) mean(x) + sd(x) / sqrt(length(x)), 
               fun.max = function(x) mean(x) - sd(x) / sqrt(length(x)), 
               width = 0.2, position = position_dodge(0.75)) + #errorbar
  
  geom_text(data = subset(Tbr1_exon_tpm_plot_Goofy, sampletype == "Tbr1_flox_1"), 
            aes(label = sprintf("%.3f", p.value), 
                y = ceiling(max(Tbr1_exon_tpm_plot_Goofy$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  geom_point(position = position_dodge(0.75), size = 1, pch = 21, color = "black", fill = "black") + #散点
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Tbr1_exon_tpm_plot_Goofy$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts")

dev.off()

##### exon (omp) #####
Tbr1_exon_featureCounts_with_length <- read.csv("featureCounts/Tbr1_exon_featureCounts_with_length_with_sample_names.csv", header = T, row.names = 1)
Tbr1_exon_featureCounts <- read.csv("featureCounts/Tbr1_exon_featureCounts_with_sample_names.csv", header = T, row.names = 1)
Tbr1_exon_length_kb <- Tbr1_exon_featureCounts_with_length$Length / 1000
Tbr1_exon_rpk <- Tbr1_exon_featureCounts / Tbr1_exon_length_kb 
Tbr1_exon_tpm <- t(t(Tbr1_exon_rpk)/colSums(featureCounts_protein_coding_and_pseudogene) * 1000000) %>% as.data.frame()
rownames(Tbr1_exon_tpm) <- paste0("exon", 1:6)

Tbr1_exon_tpm_plot_Omp <- Tbr1_exon_tpm[,grep("Omp",colnames(Tbr1_exon_tpm))]  %>% 
  tibble::rownames_to_column(var = 'exon_name') %>% 
  pivot_longer( cols = colnames(Tbr1_exon_tpm)[grep("Omp",colnames(Tbr1_exon_tpm))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  mutate(genotype = ifelse(grepl("^Omp", sampletype), "Omp_Cre", "Tbr1_flox_Omp_Cre")) %>% 
  mutate(genotype = factor(genotype, level = c("Omp_Cre", "Tbr1_flox_Omp_Cre"))) 

t_test_results <- Tbr1_exon_tpm_plot_Omp %>%
  group_by(exon_name) %>%
  summarise(t_test = list(tidy(t.test(normalized_counts ~ genotype))), .groups = 'drop') %>%
  unnest(t_test)
# # A tibble: 6 × 11
# exon_name estimate estimate1 estimate2 statistic p.value parameter conf.low conf.high method                  alternative
# <chr>        <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>     <dbl> <chr>                   <chr>      
#   1 exon1       0.0561     0.240     0.184     0.354  0.737       5.48  -0.341      0.454 Welch Two Sample t-test two.sided  
# 2 exon2       0.111      0.265     0.155     0.541  0.615       4.40  -0.438      0.660 Welch Two Sample t-test two.sided  
# 3 exon3       0.592      0.733     0.141     1.12   0.340       3.25  -1.03       2.21  Welch Two Sample t-test two.sided  
# 4 exon4       0.787      0.972     0.185     4.05   0.0180      3.71   0.230      1.34  Welch Two Sample t-test two.sided  
# 5 exon5       0.426      1.33      0.901     0.479  0.649       5.92  -1.76       2.61  Welch Two Sample t-test two.sided  
# 6 exon6       0.372      0.667     0.295     2.70   0.0609      3.57  -0.0289     0.774 Welch Two Sample t-test two.sided 

Tbr1_exon_tpm_plot_Omp <- Tbr1_exon_tpm_plot_Omp %>% 
  left_join(t_test_results[,c("exon_name","p.value")], by = "exon_name") 

pdf("./Plots/bar_plot_Tbr1_exon_normalized_counts_Omp.pdf", width = 12, height = 7)
ggplot(data = Tbr1_exon_tpm_plot_Omp, aes(exon_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.4, position = position_dodge(0.75)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(x) mean(x) + sd(x) / sqrt(length(x)), 
               fun.max = function(x) mean(x) - sd(x) / sqrt(length(x)), 
               width = 0.2, position = position_dodge(0.75)) + #errorbar
  
  geom_text(data = subset(Tbr1_exon_tpm_plot_Omp, sampletype == "Omp_Cre_1"), 
            aes(label = sprintf("%.3f", p.value), 
                y = ceiling(max(Tbr1_exon_tpm_plot_Omp$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  geom_point(position = position_dodge(0.75), size = 1, pch = 21, color = "black", fill = "black") + #散点
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(Tbr1_exon_tpm_plot_Omp$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 3)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts")
dev.off()


##### DESeq2 (Goofy) #####
### create a DESeqDataSet object
meta_Goofy <- meta[!grepl("Omp",rownames(meta)),]
featureCounts_protein_coding_and_pseudogene_Goofy <- featureCounts_protein_coding_and_pseudogene[,!grepl("Omp", colnames(featureCounts_protein_coding_and_pseudogene))] %>% 
  select(rownames(meta_Goofy))

featureCounts_protein_coding_and_pseudogene_Goofy[] <- lapply(featureCounts_protein_coding_and_pseudogene_Goofy, function(x) as.integer(round(x)))
dds_dat_Goofy <- DESeqDataSetFromMatrix(countData = featureCounts_protein_coding_and_pseudogene_Goofy, colData = meta_Goofy, design= ~ genotype)
dds_dat_Goofy$genotype <- relevel(dds_dat_Goofy$genotype, "Tbr1_flox")  

###### PCA ######
vsd <- vst(dds_dat_Goofy, blind = FALSE)
pca_data <- DESeq2::plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
# ggplot2
percentVar <- round(100 * attr(pca_data, "percentVar"))
p.pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(Tbr1_flox_Goofy_Cre = "orangered",
                                Tbr1_flox = "springgreen"))+
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
ggsave("Plots/Goofy_PCA.pdf", p.pca, width = 15, height = 10)

###### DEseq2 ######
dds_DE_Goofy <- DESeq(dds_dat_Goofy)
assays(dds_DE_Goofy)
# List of length 4
# names(4): counts mu H cooks

assays(dds_DE_Goofy)[["cooks"]]["Tbr1",]
# Tbr1_flox_Goofy_Cre_1 Tbr1_flox_Goofy_Cre_2 Tbr1_flox_Goofy_Cre_3 Tbr1_flox_Goofy_Cre_4           Tbr1_flox_1           Tbr1_flox_2           Tbr1_flox_3           Tbr1_flox_4 
# 0.02555787            0.12855590            0.19668948            0.06305629            0.09647254            0.13020586            0.26498148            0.02827729 

###### export the DESeq2 results ######
res_Goofy_cko <- results(dds_DE_Goofy, contrast = c("genotype", "Tbr1_flox_Goofy_Cre", "Tbr1_flox"), test="Wald", independentFiltering = F)
# set **independentFiltering to FALSE** so that there will be less NA values of padj (see below for padj). 
# Especially for lowly expressed genes such as olfactory receptors. 
# The independent filtering is designed only to filter out low count genes to the extent that they are not enriched with small p-values.

# PrecGoofyute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(res_Goofy_cko), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_Goofy_cko <- res_Goofy_cko %>%
  as.data.frame() %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(res_Goofy_cko, "Results/res_Goofy_cko.csv")

all.equal(
  res_Goofy_cko$padj,
  p.adjust(res_Goofy_cko$pvalue, method = "BH")
)
# TRUE
# 不同的 p 值得到相同的 padj，不代表出错；这是因为 BH 校正不是逐个独立变换，而是排序后再做单调化，尾部很多值会被压成同一个数。

quantile(res_Goofy_cko$pvalue)
# 0%          25%          50%          75%         100% 
# 1.778877e-15 4.710522e-01 6.679972e-01 8.342980e-01 9.999921e-01 

quantile(res_Goofy_cko$padj)
# 0%          25%          50%          75%         100% 
# 3.732617e-11 9.999921e-01 9.999921e-01 9.999921e-01 9.999921e-01 

rownames(res_Goofy_cko[res_Goofy_cko$padj < 0.05 & abs(res_Goofy_cko$log2FoldChange) >= 0.585,])
# [1] "Taar2"   "Taar4"   "Taar5"   "Taar6"   "Taar7a"  "Taar7b"  "Taar7d"  "Taar7e"  "Taar7f"  "Taar9"   "Wfdc18"  "Gucy1b2"

dim(res_Goofy_cko[res_Goofy_cko$padj < 0.05 & abs(res_Goofy_cko$log2FoldChange) >= 0.585,])
# [1] 12 11

###### export the normalized counts ######
normalized_count_Goofy_cko <- DESeq2::counts(dds_DE_Goofy, normalized = TRUE) %>% as.data.frame()

# PrecGoofyute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(normalized_count_Goofy_cko), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

normalized_count_Goofy_cko <- normalized_count_Goofy_cko %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(normalized_count_Goofy_cko, "Results/normalized_count_Goofy_cko.csv")

#### plot (Goofy) ####
##### plot 1: volcano plot #####
### check the range of log2FC and padj
res_Goofy_cko <- read.csv("./Results/res_Goofy_cko.csv", header = T, row.names = 1) %>% as.data.frame()

range(res_Goofy_cko$log2FoldChange, na.rm = T) 
# [1] -8.019236  5.187309
range(-log10(res_Goofy_cko$padj), na.rm = T) 
# [1] 3.433988e-06 1.042799e+01

# check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
# see which genes have padj = 0 
res_Goofy_cko[which(res_Goofy_cko$padj == 0), ] # no 

res_Goofy_cko_diff <- res_Goofy_cko[res_Goofy_cko$padj < 0.05 & 
                                      abs(res_Goofy_cko$log2FoldChange) >= 0.585, ]
dim(res_Goofy_cko_diff)
# 12 11
res_Goofy_cko_diff_functional_OR <- res_Goofy_cko[grepl(pattern = "^Olfr", rownames(res_Goofy_cko)) & 
                                                    res_Goofy_cko$gene_type %in% "protein_coding" & 
                                                    res_Goofy_cko$padj < 0.05 & 
                                                    abs(res_Goofy_cko$log2FoldChange) >= 0.585, ]
dim(res_Goofy_cko_diff_functional_OR)
# 0 11
res_Goofy_cko_diff_functional_TAAR <- res_Goofy_cko[grepl(pattern = "^Taar[2-9]", rownames(res_Goofy_cko)) & 
                                                      res_Goofy_cko$gene_type %in% "protein_coding"  & 
                                                      res_Goofy_cko$padj < 0.05 & 
                                                      abs(res_Goofy_cko$log2FoldChange) >= 0.585, ] 
dim(res_Goofy_cko_diff_functional_TAAR)
# 10 11
res_Goofy_cko_other_genes <- res_Goofy_cko %>% 
  filter(!rownames(.) %in% c(rownames(res_Goofy_cko_diff_functional_OR), c(res_Goofy_cko_diff_functional_TAAR)))

pdf("./Plots/volcano_plot_DEGs_Goofy_cre.pdf", height = 8, width = 8)
ggplot(data = res_Goofy_cko_other_genes, aes(x = log2FoldChange, y = -log10(padj), 
                                             label = gene_name, 
                                             colour = padj < 0.05 & abs(log2FoldChange) >= 0.585)) +
  geom_point(alpha = 0.75, pch = 16, size = 1) + 
  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_color_manual(values = c("grey", "black")) +
  
  geom_point(data = res_Goofy_cko_diff_functional_TAAR, color = "blue", size = 2, pch = 16)+
  geom_text_repel(data = res_Goofy_cko_diff_functional_TAAR, size = 3, colour = 'blue') + 
  
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  
  xlab("log2 fold change") + 
  ylab("-log10 padj") + 
  
  xlim(floor(min(res_Goofy_cko$log2FoldChange)), ceiling(max(res_Goofy_cko$log2FoldChange))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(-log10(res_Goofy_cko$padj))))) + 
  
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
# remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
dev.off()

##### plot 2: Taar normalized expression counts bar plot #####
normalized_count_Goofy_cko <- read.csv("./Results/normalized_count_Goofy_cko.csv", header = T, row.names = 1)
res_Goofy_cko <- read.csv("./Results/res_Goofy_cko.csv", header = T, row.names = 1) %>% as.data.frame

functional_TAAR_normalized_count_Goofy <- normalized_count_Goofy_cko[grepl(pattern = "^Taar[2-9]", rownames(normalized_count_Goofy_cko)) & normalized_count_Goofy_cko$gene_type %in% "protein_coding", ]

functional_TAAR_normalized_count_Goofy_for_plot <- as.data.frame(functional_TAAR_normalized_count_Goofy[, grepl("Tbr1",colnames(functional_TAAR_normalized_count_Goofy))]) %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(functional_TAAR_normalized_count_Goofy)[grepl("Tbr1",colnames(functional_TAAR_normalized_count_Goofy))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  left_join(res_Goofy_cko[,c("gene_name","padj")], by = "gene_name")  %>%
  mutate(genotype = ifelse(grepl("Goofy", sampletype), "Tbr1_flox_Goofy_Cre", "Tbr1_flox")) %>% 
  mutate(genotype = factor(genotype, level = c("Tbr1_flox", "Tbr1_flox_Goofy_Cre"))) 

pdf("./Plots/bar_plot_taar_normalized_counts_Goofy.pdf", width = 12, height = 7)
ggplot(data = functional_TAAR_normalized_count_Goofy_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.2, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "black") + #散点
  
  geom_text(data = subset(functional_TAAR_normalized_count_Goofy_for_plot,  sampletype == "Tbr1_flox_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(functional_TAAR_normalized_count_Goofy_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(functional_TAAR_normalized_count_Goofy_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts")

dev.off()

functional_TAAR_normalized_count_compared_to_control_Goofy_for_plot <- as.data.frame(functional_TAAR_normalized_count_Goofy[, grepl("Tbr1",colnames(functional_TAAR_normalized_count_Goofy))]/rowMeans(functional_TAAR_normalized_count_Goofy[, grepl("Tbr1_flox_\\d",colnames(functional_TAAR_normalized_count_Goofy))])) %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(functional_TAAR_normalized_count_Goofy)[grepl("Tbr1",colnames(functional_TAAR_normalized_count_Goofy))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  left_join(res_Goofy_cko[,c("gene_name","padj")], by = "gene_name")  %>%
  mutate(genotype = ifelse(grepl("Goofy", sampletype), "Tbr1_flox_Goofy_Cre", "Tbr1_flox")) %>% 
  mutate(genotype = factor(genotype, level = c("Tbr1_flox", "Tbr1_flox_Goofy_Cre"))) 

pdf("./Plots/bar_plot_taar_normalized_counts_compared_to_control_Goofy.pdf", width = 12, height = 7)
ggplot(data = functional_TAAR_normalized_count_compared_to_control_Goofy_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.2, position = position_dodge(0.9))+ #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "black")+ #散点
  
  geom_text(data = subset(functional_TAAR_normalized_count_compared_to_control_Goofy_for_plot,  sampletype == "Tbr1_flox_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(functional_TAAR_normalized_count_compared_to_control_Goofy_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(functional_TAAR_normalized_count_compared_to_control_Goofy_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts compared to control")
dev.off()

##### plot 5: log2FC and genome location (Taar gene cluster) ####
res_Goofy_cko <- read.csv("./Results/res_Goofy_cko.csv", header = T, row.names = 1)
taar_res_Goofy_cko <- res_Goofy_cko[grepl(pattern = "^Taar[2-9]|Slc18b1|Stx7|Moxd1", row.names(res_Goofy_cko)) & res_Goofy_cko$gene_type %in% "protein_coding", ]
taar_res_Goofy_cko$gene_name <- rownames(taar_res_Goofy_cko)

pdf("./Plots/taar_log2FC_with_genome_location_Goofy.pdf", width = 7, height = 7)
ggplot(data = taar_res_Goofy_cko, aes(genomic_start, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(taar_res_Goofy_cko, padj < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = T, color = "blue")+ #散点
  geom_point(data = subset(taar_res_Goofy_cko, (padj < 0.05 & abs(log2FoldChange) >= 0.585) == F), size = 2, pch = 1, na.rm = T, color = "blue")+ #散点
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(taar_res_Goofy_cko$log2FoldChange)), 2.5))+
  # floor(min(taar_res_Goofy_cko$log2FoldChange)) -8
  scale_x_continuous(expand = c(0, 0), limits = c(23790000, 24310000), breaks = seq(23800000, 24300000, by = 100000) , labels = seq(23800000, 24300000, by = 100000)/1000000)+
  
  geom_text_repel(size = 4, colour = 'blue')+
  
  theme_classic()+ ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(), 
        panel.grid = element_blank())+ #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
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
  geom_text(data = data.frame(label = gsub("Taar", "", taar_res_Goofy_cko$gene_name[grepl("Taar",taar_res_Goofy_cko$gene_name)]), 
                              x = ((taar_res_Goofy_cko$genomic_start+taar_res_Goofy_cko$genomic_end)/2)[grepl("Taar",taar_res_Goofy_cko$gene_name)], 
                              y = rep(2, 14)), 
            mapping = aes(x = x, y = y, label = label), 
            size = 2)

dev.off()

##### plot 6: log2FC with genome location (OR and Taar)####
# Read data
res_Goofy_cko <- read.csv("./Results/res_Goofy_cko.csv", header = TRUE, row.names = 1)
chr_all <- paste0("chr", c(1:19, "X", "Y"))
functional_receptor_res_Goofy_cko <- res_Goofy_cko[grepl(pattern = "Olfr|Taar[2-9]", rownames(res_Goofy_cko)) & res_Goofy_cko$gene_type %in% "protein_coding", ] %>% 
  mutate(chr_name = factor(chr_name, levels = chr_all)) %>% 
  arrange(chr_name) %>% 
  mutate(genomic_start_relative = 1:nrow(.),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(.))
write.csv(functional_receptor_res_Goofy_cko, file = "Results/olfactory_receptor_res_Goofy_cko.csv")

# Generate genome location data
genome_location <- as.data.frame(table(functional_receptor_res_Goofy_cko$chr_name))
genome_location <- genome_location[match(unique(functional_receptor_res_Goofy_cko$chr_name), genome_location$Var1), ]
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

pdf("./Plots/functional_receptors_log2FC_with_genome_location_Goofy.pdf", width = 8, height = 7)
ggplot(data = functional_receptor_res_Goofy_cko, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(functional_receptor_res_Goofy_cko, grepl("^Olfr", gene_name) == TRUE), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(functional_receptor_res_Goofy_cko, (grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) >= 0.585) == TRUE), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(functional_receptor_res_Goofy_cko, grepl("^Taar", gene_name) == TRUE), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(functional_receptor_res_Goofy_cko, (grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) >= 0.585) == TRUE), size = 2, pch = 16, na.rm = TRUE, color = "blue") +
  
  geom_text_repel(data = subset(functional_receptor_res_Goofy_cko, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin =  ceiling(max(functional_receptor_res_Goofy_cko$log2FoldChange))*0.95, ymax =  ceiling(max(functional_receptor_res_Goofy_cko$log2FoldChange)), fill = Var1), inherit.aes = FALSE) +
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(functional_receptor_res_Goofy_cko$log2FoldChange)), ceiling(max(functional_receptor_res_Goofy_cko$log2FoldChange)))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(functional_receptor_res_Goofy_cko$genomic_end_relative)), breaks = seq(0, max(functional_receptor_res_Goofy_cko$genomic_end_relative), by = 200), labels = seq(0, max(functional_receptor_res_Goofy_cko$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes")
dev.off()

#### DESeq2 (Omp) #####
### create a DESeqDataSet object
meta_Omp <- meta[grepl("Omp",rownames(meta)),]
featureCounts_protein_coding_and_pseudogene_Omp <- featureCounts_protein_coding_and_pseudogene[,grepl("Omp", colnames(featureCounts_protein_coding_and_pseudogene))] %>% 
  select(rownames(meta_Omp))
all(row.names(meta_Omp) == colnames(featureCounts_protein_coding_and_pseudogene_Omp))

featureCounts_protein_coding_and_pseudogene_Omp[] <- lapply(featureCounts_protein_coding_and_pseudogene_Omp, function(x) as.integer(round(x)))
dds_dat_Omp <- DESeqDataSetFromMatrix(countData = featureCounts_protein_coding_and_pseudogene_Omp, colData = meta_Omp, design= ~ genotype)
dds_dat_Omp$genotype <- relevel(dds_dat_Omp$genotype, "Omp_Cre")  

###### PCA ######
vsd <- vst(dds_dat_Omp, blind = FALSE)
pca_data <- DESeq2::plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
# ggplot2
percentVar <- round(100 * attr(pca_data, "percentVar"))
p.pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3.5) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(Tbr1_flox_Omp_Cre = "orangered",
                                Omp_Cre = "springgreen"))+
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
ggsave("Plots/Omp_PCA.pdf", p.pca, width = 15, height = 10)

###### DEseq2 ######
dds_DE_Omp <- DESeq(dds_dat_Omp)

###### export the DESeq2 results ######
res_Omp_cko <- results(dds_DE_Omp, contrast = c("genotype", "Tbr1_flox_Omp_Cre", "Omp_Cre"), test="Wald", independentFiltering = F)
# set **independentFiltering to FALSE** so that there will be less NA values of padj (see below for padj). 
# Especially for lowly expressed genes such as olfactory receptors. 
# The independent filtering is designed only to filter out low count genes to the extent that they are not enriched with small p-values.

# Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(res_Omp_cko), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

res_Omp_cko <- res_Omp_cko %>%
  as.data.frame() %>% 
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(res_Omp_cko, "Results/res_Omp_cko.csv")

all.equal(
  res_Omp_cko$padj,
  p.adjust(res_Omp_cko$pvalue, method = "BH")
)
# TRUE

quantile(res_Omp_cko$pvalue)
# 0%          25%          50%          75%         100% 
# 1.143159e-05 1.943557e-01 3.868584e-01 6.743715e-01 9.998926e-01 

quantile(res_Omp_cko$padj)
# 0%       25%       50%       75%      100% 
# 0.2381429 0.7460807 0.7736721 0.8990971 0.9998926

dim(res_Omp_cko[res_Omp_cko$padj < 0.05 & abs(res_Omp_cko$log2FoldChange) >= 0.585,])
# [1]  0 11

###### export the normalized counts ######
normalized_count_Omp_cko <- DESeq2::counts(dds_DE_Omp, normalized = TRUE) %>% as.data.frame()

# Precompute matches to avoid multiple calls to `match()`
match_indices <- match(rownames(normalized_count_Omp_cko), ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name)

normalized_count_Omp_cko <- normalized_count_Omp_cko %>%
  dplyr::mutate(gene_type = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_type[match_indices],
                chr_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$chr_name[match_indices],
                genomic_start = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_start[match_indices],
                genomic_end = ensembl_gene_info_protein_coding_and_pseudogene_nodup$genomic_end[match_indices],
                gene_name = ensembl_gene_info_protein_coding_and_pseudogene_nodup$gene_name[match_indices]) %>% 
  na.omit()
write.csv(normalized_count_Omp_cko, "Results/normalized_count_Omp_cko.csv")

#### plot (omp) ####
##### plot 1: volcano plot #####
### check the range of log2FC and padj
res_Omp_cko <- read.csv("./Results/res_Omp_cko.csv", header = T, row.names = 1) %>% as.data.frame
range(res_Omp_cko$log2FoldChange, na.rm = T) 
# [1] -4.601634  4.338417
range(-log10(res_Omp_cko$padj), na.rm = T) 
# [1] 4.663885e-05 6.231624e-01

# check are there genes of padj = 0, which causes issues when plotting. So change those padj to a small number.
# see which genes have padj = 0 
res_Omp_cko[which(res_Omp_cko$padj == 0), ] # no 

res_Omp_cko_diff <- res_Omp_cko[res_Omp_cko$padj < 0.05 & 
                                  abs(res_Omp_cko$log2FoldChange) >= 0.585, ]
dim(res_Omp_cko_diff)
# 0 11
res_Omp_cko_diff_functional_OR <- res_Omp_cko[grepl(pattern = "^Olfr", rownames(res_Omp_cko)) & 
                                                res_Omp_cko$gene_type %in% "protein_coding" & 
                                                res_Omp_cko$padj < 0.05 & 
                                                abs(res_Omp_cko$log2FoldChange) >= 0.585, ]
dim(res_Omp_cko_diff_functional_OR)
# 0 11
res_Omp_cko_diff_functional_TAAR <- res_Omp_cko[grepl(pattern = "^Taar[2-9]", rownames(res_Omp_cko)) & 
                                                  res_Omp_cko$gene_type %in% "protein_coding"  & 
                                                  res_Omp_cko$padj < 0.05 & 
                                                  abs(res_Omp_cko$log2FoldChange) >= 0.585, ] 
dim(res_Omp_cko_diff_functional_TAAR)
# 0 11
res_Omp_cko_other_genes <- res_Omp_cko %>% 
  filter(!rownames(.) %in% c(rownames(res_Omp_cko_diff_functional_OR), c(res_Omp_cko_diff_functional_TAAR)))

pdf("./Plots/volcano_plot_DEGs_Omp_cre.pdf", height = 8, width = 8)
ggplot(data = res_Omp_cko_other_genes, aes(x = log2FoldChange, y = -log10(padj), 
                                           label = gene_name, 
                                           colour = padj < 0.05 & abs(log2FoldChange) >= 0.585)) +
  geom_point(alpha = 0.75, pch = 16, size = 1) + 
  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_color_manual(values = c("grey", "black")) +
  
  geom_point(data = res_Omp_cko_diff_functional_TAAR, color = "blue", size = 2, pch = 16)+
  geom_text_repel(data = res_Omp_cko_diff_functional_TAAR, size = 3, colour = 'blue') + 
  
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  
  xlab("log2 fold change") + 
  ylab("-log10 padj") + 
  
  xlim(floor(min(res_Omp_cko$log2FoldChange)), ceiling(max(res_Omp_cko$log2FoldChange))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(-log10(res_Omp_cko$padj))))) + 
  
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.75), vjust = 1), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
# remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
dev.off()

##### plot 2: Taar normalized expression counts bar plot #####
normalized_count_Omp_cko <- read.csv("./Results/normalized_count_Omp_cko.csv", header = T, row.names = 1)
res_Omp_cko <- read.csv("./Results/res_Omp_cko.csv", header = T, row.names = 1) %>% as.data.frame

functional_TAAR_normalized_count_Omp <- normalized_count_Omp_cko[grepl(pattern = "^Taar[2-9]", rownames(normalized_count_Omp_cko)) & normalized_count_Omp_cko$gene_type %in% "protein_coding", ]

functional_TAAR_normalized_count_Omp_for_plot <- as.data.frame(functional_TAAR_normalized_count_Omp[, grepl("Cre",colnames(functional_TAAR_normalized_count_Omp))]) %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(functional_TAAR_normalized_count_Omp)[grepl("Cre",colnames(functional_TAAR_normalized_count_Omp))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  left_join(res_Omp_cko[,c("gene_name","padj")], by = "gene_name")  %>%
  mutate(genotype = ifelse(grepl("Tbr1", sampletype), "Tbr1_flox_Omp_Cre", "Omp_Cre")) %>% 
  mutate(genotype = factor(genotype, level = c("Omp_Cre", "Tbr1_flox_Omp_Cre"))) 

pdf("./Plots/bar_plot_taar_normalized_counts_Omp.pdf", width = 12, height = 7)
ggplot(data = functional_TAAR_normalized_count_Omp_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.2, position = position_dodge(0.9)) + #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "black") + #散点
  
  geom_text(data = subset(functional_TAAR_normalized_count_Omp_for_plot,  sampletype == "Omp_Cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(functional_TAAR_normalized_count_Omp_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(functional_TAAR_normalized_count_Omp_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts")

dev.off()

functional_TAAR_normalized_count_compared_to_control_Omp_for_plot <- as.data.frame(functional_TAAR_normalized_count_Omp[, grepl("Cre",colnames(functional_TAAR_normalized_count_Omp))]/rowMeans(functional_TAAR_normalized_count_Omp[, grepl("^Omp",colnames(functional_TAAR_normalized_count_Omp))])) %>% 
  tibble::rownames_to_column(var = 'gene_name') %>% 
  pivot_longer( cols = colnames(functional_TAAR_normalized_count_Omp)[grepl("Cre",colnames(functional_TAAR_normalized_count_Omp))], 
                names_to = 'sampletype', 
                values_to = 'normalized_counts') %>% 
  left_join(res_Omp_cko[,c("gene_name","padj")], by = "gene_name")  %>%
  mutate(genotype = ifelse(grepl("Tbr1", sampletype), "Tbr1_flox_Omp_Cre", "Omp_Cre")) %>% 
  mutate(genotype = factor(genotype, level = c("Omp_Cre", "Tbr1_flox_Omp_Cre"))) 

pdf("./Plots/bar_plot_taar_normalized_counts_compared_to_control_Omp.pdf", width = 12, height = 7)
ggplot(data = functional_TAAR_normalized_count_compared_to_control_Omp_for_plot, aes(gene_name, normalized_counts, group = genotype)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill=genotype, color=genotype), width = 0.8, position = position_dodge(0.9)) + # bar
  stat_summary(geom = "errorbar", 
               fun.min = function(group) mean(group) + sd(group) / sqrt(length(group)), 
               fun.max = function(group) mean(group) - sd(group) / sqrt(length(group)), 
               width = 0.2, position = position_dodge(0.9))+ #errorbar
  
  geom_point(position = position_dodge(0.9), size = 1, pch = 21, color = "black", fill = "black")+ #散点
  
  geom_text(data = subset(functional_TAAR_normalized_count_compared_to_control_Omp_for_plot,  sampletype == "Omp_Cre_1"), 
            aes(label = sprintf("%.3f", padj), 
                y = ceiling(max(functional_TAAR_normalized_count_compared_to_control_Omp_for_plot$normalized_counts))*0.97), 
            position = position_dodge(width = 0.75), 
            size =  2) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(functional_TAAR_normalized_count_compared_to_control_Omp_for_plot$normalized_counts))))+ 
  scale_fill_manual(values = rep("white", 2)) + #改变柱状图轮廓的颜色
  scale_color_manual(values = c("springgreen3", "red1")) + #改变柱状图内填充的颜色
  
  theme_classic() +  ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(0.75),
                                 angle = 90),
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.75)),
        legend.title = element_blank()) +
  labs(x = "", y = "normalized counts compared to control")
dev.off()

##### plot 5: log2FC and genome location (Taar gene cluster) ####
res_Omp_cko <- read.csv("./Results/res_Omp_cko.csv", header = T, row.names = 1)
taar_res_Omp_cko <- res_Omp_cko[grepl(pattern = "^Taar[2-9]|Slc18b1|Stx7|Moxd1", row.names(res_Omp_cko)) & res_Omp_cko$gene_type %in% "protein_coding", ]
taar_res_Omp_cko$gene_name <- rownames(taar_res_Omp_cko)

pdf("./Plots/taar_log2FC_with_genome_location_Omp.pdf", width = 7, height = 7)
ggplot(data = taar_res_Omp_cko, aes(genomic_start, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(taar_res_Omp_cko, padj < 0.05 & abs(log2FoldChange) >= 0.585), size = 2, pch = 16, na.rm = T, color = "blue")+ #散点
  geom_point(data = subset(taar_res_Omp_cko, (padj < 0.05 & abs(log2FoldChange) >= 0.585) == F), size = 2, pch = 1, na.rm = T, color = "blue")+ #散点
  
  scale_y_continuous(expand = c(0, 0), limits = c(-8, 2.5))+
  # floor(min(taar_res_Goofy_cko$log2FoldChange)) -8
  scale_x_continuous(expand = c(0, 0), limits = c(23790000, 24310000), breaks = seq(23800000, 24300000, by = 100000) , labels = seq(23800000, 24300000, by = 100000)/1000000)+
  
  geom_text_repel(size = 4, colour = 'blue')+
  
  theme_classic()+ ##可以去除顶部和右边的边框
  theme(legend.title = element_blank(), 
        panel.grid = element_blank())+ #底色改为白色无网格线
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
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
  geom_text(data = data.frame(label = gsub("Taar", "", taar_res_Omp_cko$gene_name[grepl("Taar",taar_res_Omp_cko$gene_name)]), 
                              x = ((taar_res_Omp_cko$genomic_start+taar_res_Omp_cko$genomic_end)/2)[grepl("Taar",taar_res_Omp_cko$gene_name)], 
                              y = rep(2, 14)), 
            mapping = aes(x = x, y = y, label = label), 
            size = 2)

dev.off()

##### plot 6: log2FC with genome location (OR and Taar)####
# Read data
res_Omp_cko <- read.csv("./Results/res_Omp_cko.csv", header = TRUE, row.names = 1)
chr_all <- paste0("chr", c(1:19, "X", "Y"))

functional_receptor_res_Omp_cko <- res_Omp_cko[grepl(pattern = "Olfr|Taar[2-9]", rownames(res_Omp_cko)) & res_Omp_cko$gene_type %in% "protein_coding", ] %>% 
  mutate(chr_name = factor(chr_name, levels = chr_all)) %>% 
  arrange(chr_name) %>% 
  mutate(genomic_start_relative = 1:nrow(.),
         genomic_end_relative = genomic_start_relative + 1,
         gene_name = rownames(.))
write.csv(functional_receptor_res_Omp_cko, file = "Results/olfactory_receptor_res_Omp_cko.csv")

# Generate genome location data
genome_location <- as.data.frame(table(functional_receptor_res_Omp_cko$chr_name))
genome_location <- genome_location[match(unique(functional_receptor_res_Omp_cko$chr_name), genome_location$Var1), ]
genome_location$start <- 0
genome_location <- genome_location %>%
  mutate(end = cumsum(Freq), # Create a cumulative frequency column
         start = if_else(row_number() == 1, start, start + lag(end))) 
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)
col_for_chr <- col_for_chr_universal

# obtain the colors for chr used in step1.1
col_for_chr <- col_for_chr[match(chr_all, unique(genome_location$Var1))]
col_for_chr <- col_for_chr[!is.na(col_for_chr)]

pdf("./Plots/functional_receptors_log2FC_with_genome_location_Omp.pdf", width = 8, height = 7)
ggplot(data = functional_receptor_res_Omp_cko, aes(genomic_start_relative, log2FoldChange, label = gene_name)) + 
  geom_point(data = subset(functional_receptor_res_Omp_cko, grepl("^Olfr", gene_name) == TRUE), size = 2, pch = 1, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(functional_receptor_res_Omp_cko, (grepl("^Olfr", gene_name) & padj < 0.05 & abs(log2FoldChange) >= 0.585) == TRUE), size = 2, pch = 16, na.rm = TRUE, color = "red") + 
  geom_point(data = subset(functional_receptor_res_Omp_cko, grepl("^Taar", gene_name) == TRUE), size = 2, pch = 1, na.rm = TRUE, color = "blue") + 
  geom_point(data = subset(functional_receptor_res_Omp_cko, (grepl("^Taar", gene_name) & padj < 0.05 & abs(log2FoldChange) >= 0.585) == TRUE), size = 2, pch = 16, na.rm = TRUE, color = "blue") +
  
  geom_text_repel(data = subset(functional_receptor_res_Omp_cko, grepl("^Taar", gene_name) == TRUE), size = 4, colour = 'blue', nudge_x = -0.2, nudge_y = -0.2,max.overlaps = 20) +
  
  geom_rect(data = genome_location, aes(xmin = start, xmax = end, ymin =  4*0.95, ymax = 4, fill = Var1), inherit.aes = FALSE) +
  # ceiling(max(functional_receptor_res_Goofy_cko$log2FoldChange)) 4
  scale_fill_manual(values = col_for_chr) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(-8, 4)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(functional_receptor_res_Omp_cko$genomic_end_relative)), breaks = seq(0, max(functional_receptor_res_Omp_cko$genomic_end_relative), by = 200), labels = seq(0, max(functional_receptor_res_Omp_cko$genomic_end_relative), by = 200)) +
  
  theme_classic() +
  theme(plot.title = element_text(size = rel(0.75)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1)), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = rel(0.8)), 
        legend.title = element_blank()) +
  
  labs(x = "Relative chromosome coordinates", y = "log2 Gene expression changes")
dev.off()



