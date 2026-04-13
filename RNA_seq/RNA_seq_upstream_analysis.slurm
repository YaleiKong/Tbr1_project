################################################################################################################
# sbatch_trim_align_featurecount.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch_trim_align
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

inputFolder="raw_data"
outputFolder="processed_trimmomatic_hisat2_featurecount"

for input in $inputFolder/WT*_R1.fastq.gz
do
    cell=${input##*/}
    cell=${cell%%_R1.fastq.gz}
    echo $cell

    fileR1=$input
    fileR2=${input/R1/R2}
    echo $fileR1
    echo $fileR2

    # mkdir $outputFolder/$cell

    sbatch ./trim_align_featurecount.slurm $cell $outputFolder $fileR1 $fileR2
done
################################################################################################################
# sbatch_trim_align_featurecount.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=trim_align_featurecount
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate trimmomatic

#####
fastqc $3 $4 -o $2/$1

#####
trimmomatic PE -phred33 -threads 10 \
$3 $4 \
$2/$1/"${1}_R1P.fq.gz" $2/$1/"${1}_R1U.fq.gz" \
$2/$1/"${1}_R2P.fq.gz" $2/$1/"${1}_R2U.fq.gz" \
ILLUMINACLIP:/lustre/home/acct-medlqian/medlqian-loop3/.conda/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50

#####
fastqc $2/$1/"${1}_R1P.fq.gz" $2/$1/"${1}_R1U.fq.gz" $2/$1/"${1}_R2P.fq.gz" $2/$1/"${1}_R2U.fq.gz" -o $2/$1
conda deactivate

#####
source activate hisat2 # version 2.2.1
module load samtools

hisat2 -p 10 --rna-strandness RF \
-x /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/hisat2_index_Gencode_GRCm38_primary_assembly_genome/genome_GRCm38_primary_assembly_Gencode \
-1 $2/$1/"${1}_R1P.fq.gz"  -2 $2/$1/"${1}_R2P.fq.gz" -U $2/$1/"${1}_R1U.fq.gz",$2/$1/"${1}_R2U.fq.gz" \
-S $2/$1/"${1}.trimmomatic.hisat2.sam"

#####
# find chr([1-9]|1[0-9]|X|Y) and =; add @ SQ lines present in the header
awk 'BEGIN {OFS = FS = "\t"} /^@/ || (/chr([1-9]|1[0-9]|X|Y)/ && /=/)' $2/$1/"${1}.trimmomatic.hisat2.sam" > $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.sam"
# BEGIN {OFS = FS = "\t"}: This initializes awk by setting the Output Field Separator (OFS) and the Field Separator (FS) to a tab character. This is done because SAM files are tab-delimited.
# /^@/: This matches any line starting with "@", which in SAM files are header lines. These lines are included in the output.
# /chr([1-9]|1[0-9]|X|Y)/: This part uses a regular expression to match lines where the chromosome field is 'chr1' to 'chr19', 'chrX', or 'chrY'. 
# && /=/: This checks for lines containing the '=' character, which in the context of a SAM file typically represents a match to the reference chromosome. 

samtools view -bS -F 4 -q 20 $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.sam" -o $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.bam"

samtools sort $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.bam" -o $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.bam" 

samtools index $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.bam"
 
conda deactivate

##########
source activate picard
picard MarkDuplicates \
      I=$2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.bam"   \
      O=$2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam" \
      REMOVE_DUPLICATES=true \
      ASSUME_SORT_ORDER=coordinate \
      M=$2/$1/"${1}_remove_dup_metrics.txt"

conda deactivate

########
source activate subread
featureCounts -T 10 -p -C -s 2 \
-O --fraction \
-a /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/gencode.vM25.Comprehensive_gene_annotation_ALL.chr_patch_hapl_scaff.annotation.gtf \
-t exon -g gene_id \
-o $2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.featureCounts.multiOverlap.txt" \
$2/$1/"${1}.trimmomatic.hisat2.excluede_chrM.mapped.Q20.sorted.bam" 

conda deactivate

################################################################################################################
# sbatch_trim_align_featurecount.slurm
################################################################################################################
