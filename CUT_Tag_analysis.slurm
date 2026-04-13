################################################################################################################
# do_trim_align.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=do_trim_align
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

inputFolder="raw_data"
outputFolder="processed_trimmomatic_bowtie2"

for input in $inputFolder/*_R1_001.fastq.gz
do
    cell=${input##*/}
    cell=${cell%%_S*}
    echo $cell

    file1R1=$input
    file1R2=${input/_R1/_R2}
    echo $file1R1
    echo $file1R2

    mkdir $outputFolder/$cell

    sbatch ./trim_align.slurm $cell $outputFolder $file1R1 $file1R2
done
################################################################################################################
# trim_align.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=trim_align
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#############################
#### fastqc 
#############################
module load fastqc/0.11.7-gcc-9.2.0
fastqc $3 $4

#############################
#### trimmomatic 
#############################
module load miniconda3/4.8.2
source activate trimmomatic
module load jdk/12.0.2_10-gcc-9.2.0

java -jar /lustre/home/acct-medlqian/medlqian-loop3/.conda/envs/trimmomatic/share/trimmomatic-0.39-2/trimmomatic.jar PE -phred33 \
$3 $4 $2/$1/"${1}_R1P.fq.gz" $2/$1/"${1}_R1U.fq.gz" \
$2/$1/"${1}_R2P.fq.gz" $2/$1/"${1}_R2U.fq.gz" \
ILLUMINACLIP:/lustre/home/acct-medlqian/medlqian-loop3/.conda/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50

conda deactivate

#############################
#### align to mouse genome 
#############################
module load miniconda3
source activate bowtie2

module load samtools

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 -p 10 \
-x /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/bowtie2_index_Gencode_GRCm38_primary_assembly_genome/genome_GRCm38_primary_assembly_Gencode \
-1 $2/$1/"${1}_R1P.fq.gz"  -2 $2/$1/"${1}_R2P.fq.gz" -U $2/$1/"${1}_R1U.fq.gz",$2/$1/"${1}_R2U.fq.gz" \
-S $2/$1/"${1}_trimmomatic_bowtie2.sam"  &>  $2/$1/"${1}_trimmomatic_bowtie2.txt" 

samtools view -bS -F 4 -q 20 $2/$1/"${1}_trimmomatic_bowtie2.sam" -o $2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.bam"

samtools sort $2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.bam" -o $2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.bam" 

samtools index $2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.bam" 

rm $2/$1/"${1}_trimmomatic_bowtie2.sam" 


################################
#### picard to remove duplicates
################################
source activate picard

#### remove duplicates
picard MarkDuplicates \
      I=$2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.bam"   \
      O=$2/$1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.bam" \
      REMOVE_DUPLICATES=true \
      ASSUME_SORT_ORDER=coordinate \
      M=$2/$1/"${1}_remove_dup_metrics.txt"

################################################################################################################
# sbatch_bam_to_bw_BPM.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=bamCoverage 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err


for input in */*_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.bam}
    echo $cell

    sbatch ./bam_to_bw_BPM.slurm $cell $input 
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
samtools sort $1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.bam" -o $1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.bam"
samtools index $1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.bam"

module load miniconda3
source activate deeptools
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
bamCoverage -b $1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.bam" \
--normalizeUsing BPM \
--effectiveGenomeSize 2652783500 \
-o $1/"${1}_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.BPM.bw" \
--numberOfProcessors 10 \
--binSize 10 \
--smoothLength 30 \
--ignoreForNormalization chrX \
--extendReads \
--centerReads

conda deactivate

################################################################################################################
# calculate_enhancer_depth_rpm.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=calculate_enhancer_depth_rpm.slurm
#SBATCH --partition=cpu 
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools
bedFile=/lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed
for input in */*_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.bam
do
    cell=${input##*/}
    cell=${cell%%_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.bam}
    
    samtools depth -aa -b $bedFile \
    $input -o "${cell}/${cell}.enhancer.samtools.depth.txt"
    
    # Process each region in the BED file
    while read -r chr start end name; do
        sum=0
        # Process each line in the depth file
        while read -r depthChr depthPos depthVal; do
            # Check if the chromosome matches and the position is within the range
            if [[ $depthChr == $chr && $depthPos -ge $start && $depthPos -le $end ]]; then
                # Add the depth to the sum
                sum=$(($sum + $depthVal))
            fi
        done < "${cell}/${cell}.enhancer.samtools.depth.txt"
        # Print the result for the current region
        echo -e "$chr\t$start\t$end\t$name\t$sum"
    done < "$bedFile" > "${cell}/${cell}.enhancer.samtools.depth.sum.txt"

    awk '{print $1, $2, $3, $4, $5/($3-$2)}'  "${cell}/${cell}.enhancer.samtools.depth.sum.txt" > "${cell}/${cell}.enhancer.samtools.depth.sum.dividedByLength.txt"

    ### normalize rpm
    input_library_size=$(samtools view -c -F 4 $input) 
    echo $input_library_size
    scale_factor=$(echo "scale=6; 1000000/$input_library_size" | bc)
    echo $scale_factor
    awk -v scale_factor="$scale_factor" 'BEGIN {OFS="\t"} {print $0, $5 * scale_factor}' "${cell}/${cell}.enhancer.samtools.depth.sum.dividedByLength.txt" > "${cell}/${cell}.enhancer.samtools.depth.sum.dividedByLength.RPM.txt"
done

################################################################################################################
# OR_gene_body_BPM_heatmap_profile.slurm
################################################################################################################
#!/bin/bash
#SBATCH --job-name=peak_heatmap
#SBATCH --partition=cpu
#SBATCH -n 20
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate deeptools

grep "Olfr" /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/gencode.vM25.Comprehensive_gene_annotation_ALL.chr_patch_hapl_scaff.annotation.gtf >> heatmap_profile/OR.gtf
grep "protein_coding" heatmap_profile/OR.gtf >> heatmap_profile/functional_OR.gtf

computeMatrix scale-regions -S processed_trimmomatic_bowtie2/En1-Tg_IgG/En1-Tg_IgG_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.BPM.bw processed_trimmomatic_bowtie2/En1-Tg_TBR1/En1-Tg_TBR1_trimmomatic_bowtie2.mapped.Q20.sorted.remove_duplicates.sort.BPM.bw \
--regionsFileName heatmap_profile/functional_OR.gtf \
--regionBodyLength 1000 \
--startLabel TSS --endLabel TES \
--beforeRegionStartLength 1000 --afterRegionStartLength 1000 \
--missingDataAsZero --skipZeros \
-o heatmap_profile/OR_TSS_TES_BPM_heatmap_profile.matrix.gz

plotHeatmap -m heatmap_profile/OR_TSS_TES_BPM_heatmap_profile.matrix.gz --colorMap RdYlBu_r -out heatmap_profile/OR_TSS_TES_BPM_heatmap_profile_adjust.pdf --samplesLabel En1-Tg_IgG En1-Tg_TBR1 --yMax 1 \
--heatmapHeight 10 --heatmapWidth 3 

conda deactivate
