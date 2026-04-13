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
#SBATCH -n 5 
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3

#########
# fastp
#########
source activate fastp

fastp --in1 $3 --in2 $4 --out1 $2/$1/"${1}_R1.trimmed.fq.gz" --out2 $2/$1/"${1}_R2.trimmed.fq.gz" \
 -g -q 20 -u 40 -l 30 -g -h $2/$1/"${1}.fastp.html" &> $2/$1/"${1}.fastp.log"

# removes the polyG tails that arise from lack of signal in NextSeq/NovaSeq technologies. 
# -g: It is enabled for Nextseq/Novaseq data by default, and you can specify -g to enable it for any data, or specify -G to disable it.  
# -q: Quality threshold per base required. Default: 15, which means that a Phred quality score of at least 15 is required
# -u: Percent of bases allowed to be below the quality threshold to keep the read (0~100). Default 40 means 40% bases can fail the quality threshold. If more bases fail, the read is removed.
# -l 50: this specifies that if a read is shorter than 50 basepairs after all filters, it should be removed.

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

# --end-to-end       entire read must align; no clipping (on)
#    OR
# --local            local alignment; ends might be soft clipped (off); trim some read characters to maximize the alignment score

# Presets:                 Same as:
#   For --end-to-end:
#    --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
#    --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
#    --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
#    --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# 
#   For --local:
#    --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#    --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#    --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
#    --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

# Paired-end:
      # -I/--minins <int>  minimum fragment length (0)
      # -X/--maxins <int>  maximum fragment length (500)
      # --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
      # --no-mixed         suppress unpaired alignments for paired reads
      # --no-discordant    suppress discordant alignments for paired reads
      # --dovetail         concordant when mates extend past each other
      # --no-contain       not concordant when one mate alignment contains other
      # --no-overlap       not concordant when mates overlap at all
# Performance:
      # -p/--threads <int> number of alignment threads to launch (1)

##############################
# samtools (exclude chrM)
##############################
module load samtools
# single core

# find chr([1-9]|1[0-9]|X|Y) and =; add @ SQ lines present in the header
awk 'BEGIN {OFS = FS = "\t"} /^@/ || (/chr([1-9]|1[0-9]|X|Y)/ && /=/)' $2/$1/"${1}.fastp.bowtie2.sam" > $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"
# BEGIN {OFS = FS = "\t"}: This initializes awk by setting the Output Field Separator (OFS) and the Field Separator (FS) to a tab character. This is done because SAM files are tab-delimited.
# /^@/: This matches any line starting with "@", which in SAM files are header lines. These lines are included in the output.
# /chr([1-9]|1[0-9]|X|Y)/: This part uses a regular expression to match lines where the chromosome field is 'chr1' to 'chr19', 'chrX', or 'chrY'. 
# && /=/: This checks for lines containing the '=' character, which in the context of a SAM file typically represents a match to the reference chromosome. 

rm $2/$1/"${1}.fastp.bowtie2.sam"

samtools view -bS -F 4 -q 20 $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"  -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam"
# -S: In older versions of SAMtools, this flag indicated the input is in SAM format. However, in more recent versions of SAMtools (1.0 and later), the -S option has been deprecated because the tool automatically detects the input format (SAM or BAM).
# -b: This option tells samtools view to output the result in BAM format. BAM is a binary format that is more compressed than SAM, which is plain text.
# -F 4: exclude segment unmapped
# -q 20: “Filter on read mapping quality (phred scale)": >=20

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

#### check Insert Sizes
# picard CollectInsertSizeMetrics \
#       I=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam" \
#       O=$2/$1/"${1}_insert_size_metrics.txt" \
#       H=$2/$1/"${1}_insert_size_histogram.pdf" \
#       M=0.5
# Note: If processing a small file, set the minimum percentage option (M) to 0.5, otherwise an error may occur.

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
# -j ATAC mode 
# -y Keep unpaired alignments
# -r remove pcr duplicates
# -e ChrM remove mitochondrial reads
# -v verbose mode
# -f log
# "If a -q threshold is not specified, q-values are not calculated (reported as -1). "
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
