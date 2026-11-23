################################################################################################################
# sbatch_fastp_bowtie2_samtools_picard_genrich.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

inputFolder="raw_data"
outputFolder="processed_data"

for input in $inputFolder/*R1.fastq.gz
do
    cell=${input##*/}
    cell=${cell%%_R1.fastq.gz}
    echo $cell

    fileR1=$input
    fileR2=${input/_R1/_R2}
    echo $fileR1
    echo $fileR2   

    mkdir $outputFolder/$cell 
    
    sbatch ./fastp_bowtie2_samtools_picard_genrich.slurm $cell $outputFolder $fileR1 $fileR2
done
################################################################################################################
# fastp_bowtie2_samtools_picard_genrich.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=trim_align_genrich
#SBATCH --partition=cpu 
#SBATCH -n 10 
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3

#########
# fastp
#########
source activate fastp

fastp --in1 $3 --in2 $4 --out1 $2/$1/"${1}_R1.trimmed.fq.gz" --out2 $2/$1/"${1}_R2.trimmed.fq.gz" \
 -g -q 20 -u 40 -l 30 -g -h $2/$1/"${1}.fastp.html" &> $2/$1/"${1}.fastp.log"

conda deactivate

#########
# QC
#########
source activate trimmomatic
fastqc $3 $4
fastqc  $2/$1/"${1}_R1.trimmed.fq.gz"  $2/$1/"${1}_R2.trimmed.fq.gz" 
conda deactivate

#########
# bowtie2
#########
source activate bowtie2
bowtie2 --end-to-end --very-sensitive \
 -x /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/bowtie2_index_Gencode_GRCm38_primary_assembly_genome/genome_GRCm38_primary_assembly_Gencode \
-1 $2/$1/"${1}_R1.trimmed.fq.gz"  -2 $2/$1/"${1}_R2.trimmed.fq.gz" \
-X 2000 -p 10 -S $2/$1/"${1}.fastp.bowtie2.sam" &> $2/$1/"${1}.bowtie2.log"
conda deactivate


##############################
# samtools (exclude chrM)
##############################
module load samtools

# find chr([1-9]|1[0-9]|X|Y) and =; add @ SQ lines present in the header
awk 'BEGIN {OFS = FS = "\t"} /^@/ || (/chr([1-9]|1[0-9]|X|Y)/ && /=/)' $2/$1/"${1}.fastp.bowtie2.sam" > $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"

rm $2/$1/"${1}.fastp.bowtie2.sam"

samtools view -bS -F 4 -q 20 $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"  -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam"

samtools sort $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam" -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"

samtools index $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"

################################
#### picard to remove duplicates
################################
source activate picard

#### remove duplicates
picard MarkDuplicates \
      I=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"  \
      O=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam" \
      REMOVE_DUPLICATES=true \
      ASSUME_SORT_ORDER=coordinate \
      M=$2/$1/"${1}_remove_dup_metrics.txt"

conda deactivate

################################
#### genrich 
################################
samtools sort -n $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam" -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.bam"
samtools index $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.bam"

source activate genrich

Genrich  -t $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.bam"  -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.genrich.qvaluelessThan0.05.narrowpeak" \
 -j  -y  -r  -e chrM  -p 0.05 -q 0.05 \
 -v -f $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.genrich.qvaluelessThan0.05.log"

conda deactivate

################################################################################################################
# sbatch_bam_to_bw_BPM.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=bamCoverage 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

for input in processed_data/*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam}
    echo $cell

    sbatch ./bam_to_bw.slurm $cell $input 
done

################################################################################################################
# bam_to_bw_BPM.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=bw 
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools 
samtools sort $2 -o processed_data/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bam"
samtools index processed_data/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bam"

module load miniconda3
source activate deeptools
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
bamCoverage -b processed_data/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bam" \
--normalizeUsing BPM \
--effectiveGenomeSize 2652783500 \
-o processed_data/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw" \
--numberOfProcessors 10 \
--binSize 10 \
--smoothLength 30 \
--ignoreForNormalization chrX \
--extendReads \
--centerReads

conda deactivate

################################################################################################################
# enhancer_heatmap_profile.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=heatmap_profile
#SBATCH --partition=cpu
#SBATCH -n 20
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate deeptools

# 63 GI enhancers
computeMatrix reference-point -S processed_data/Lsl_Tbr1_ires_tdTomato_Het_1/Lsl_Tbr1_ires_tdTomato_Het_1.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Lsl_Tbr1_ires_tdTomato_Homo_3/Lsl_Tbr1_ires_tdTomato_Homo_3.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
-R /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers.bed \
--referencePoint center -b 2000 -a 2000 -bs 10 \
-out heatmap_profile/OR_63_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz

plotHeatmap -m heatmap_profile/OR_63_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz \
-out heatmap_profile/OR_63_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.pdf \
--refPointLabel "center" --samplesLabel Lsl_Tbr1_ires_tdTomato_Het_1 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2 Lsl_Tbr1_ires_tdTomato_Homo_3 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4 \
--regionsLabel "GI enhancers" \
--colorMap RdYlBu_r --heatmapWidth 3 --heatmapHeight 10

# J element
computeMatrix reference-point -S processed_data/Lsl_Tbr1_ires_tdTomato_Het_1/Lsl_Tbr1_ires_tdTomato_Het_1.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Lsl_Tbr1_ires_tdTomato_Homo_3/Lsl_Tbr1_ires_tdTomato_Homo_3.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
-R /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/J_element.bed \
--referencePoint center -b 2000 -a 2000 -bs 10 \
-out heatmap_profile/J_element.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz

plotHeatmap -m heatmap_profile/J_element.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz \
-out heatmap_profile/J_element.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.pdf \
--refPointLabel "center" --samplesLabel Lsl_Tbr1_ires_tdTomato_Het_1 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2 Lsl_Tbr1_ires_tdTomato_Homo_3 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4 \
--regionsLabel "J element" \
--colorMap RdYlBu_r --heatmapWidth 3 --heatmapHeight 5

# Taar enhancers
computeMatrix reference-point -S processed_data/Lsl_Tbr1_ires_tdTomato_Het_1/Lsl_Tbr1_ires_tdTomato_Het_1.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Lsl_Tbr1_ires_tdTomato_Homo_3/Lsl_Tbr1_ires_tdTomato_Homo_3.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
processed_data/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4/Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.bw \
-R /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/Taar_enhancers.bed \
--referencePoint center -b 2000 -a 2000 -bs 10 \
-out heatmap_profile/Taar_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz

plotHeatmap -m heatmap_profile/Taar_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.matrix.gz \
-out heatmap_profile/Taar_enhancers.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.sorted.pdf \
--refPointLabel "center" --samplesLabel Lsl_Tbr1_ires_tdTomato_Het_1 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Het_2 Lsl_Tbr1_ires_tdTomato_Homo_3 Goofy_Cre_Lsl_Tbr1_ires_tdTomato_Homo_4 \
--regionsLabel "Taar enhancers" \
--colorMap RdYlBu_r --heatmapWidth 3 --heatmapHeight 6

conda deactivate
