#!/usr/bin
# Last modified Nov 11 2025
# 12 CPU, 40G
# conda activate polyfun (rnaseqc)

#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=7-0:0:0


module load gcc/12.3
module load star/2.7.11a
module load samtools/1.18 
module load picard/3.1.0

thread=12
sample=$1

##Step 1. STAR: Build index (this command). Run once for a given genome + annotation.
## hg38.fa from UCSC genome rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
## gencode.v26.annotation.gtf from gencode
## This should be read length – 1 (for 151 bp paired-end reads, use 150). 
## Check read length from fastq: zcat xx.fastq.gz | head -n 2 | tail -n 1 | wc -c ; 151 bp read length
## cd /data/sbcs/GuoLab/backup/liq17/ref/STAR_index
# salloc -N 1 -n 1 -c 8 --mem=80G -t 8:00:00 # 1h
# fa_file='/data/sbcs/GuoLab/backup/liq17/ref/UCSC/hg38.fa'
# gtf_file='/data/sbcs/GuoLab/backup/liq17/ref/gencode/gencode.v26.annotation.gtf'
# STAR --runMode genomeGenerate \
     # --genomeDir ./STAR_150bpindex_hg38/ \
     # --genomeFastaFiles ${fa_file} \
     # --sjdbGTFfile ${gtf_file}  \
     # --sjdbOverhang 150 \
     # --runThreadN 8


############################### setting ###############################################
WORK_DIR='/data/sbcs/GuoLab/backup/RNA-seq/AfricanAmerican/processed/Expr'
echo "........................"
echo "Run analysis in $WORK_DIR"
echo "........................"

input="/data/sbcs/GuoLab/backup/RNA-seq/AfricanAmerican/batch_202409/fastq_files/usftp21.novogene.com/01.RawData/"

FQ1=$input/$sample/${sample}_1.fq.gz
FQ2=$input/$sample/${sample}_2.fq.gz

# STAR index path
index='/data/sbcs/GuoLab/backup/liq17/ref/STAR_index/STAR_150bpindex_hg38' 

# Create output directories
mkdir -p $WORK_DIR/$sample
output=$WORK_DIR/$sample
prefix="$sample."

################################ run #################################################
START=$(date +%s)

STAR --runMode alignReads \
     --runThreadN ${thread} \
     --genomeDir ${index} \
     --twopassMode Basic \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFilterType BySJout \
     --outFilterScoreMinOverLread 0.33 \
     --outFilterMatchNminOverLread 0.33 \
     --limitSjdbInsertNsj 1200000 \
     --readFilesIn $FQ1 $FQ2 \
     --readFilesCommand zcat \
     --outFileNamePrefix ${output}/${prefix} \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs None \
     --alignSoftClipAtReferenceEnds Yes \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMtype BAM Unsorted \
     --outSAMunmapped Within \
     --genomeLoad NoSharedMemory \
     --outSAMattributes NH HI AS nM NM ch \
     --outSAMattrRGline ID:rg1 SM:sm1

# Sort and index BAM
samtools sort --threads ${thread} ${output}/${prefix}Aligned.out.bam -o ${output}/${prefix}Aligned.sortedByCoord.out.bam
samtools index ${output}/${prefix}Aligned.sortedByCoord.out.bam

# Cleanup
rm ${output}/${prefix}Aligned.out.bam

# Save splice junction info
# First line: keeps the first-pass splice junctions around.
# Next three lines: compress splice junction + count outputs from STAR so you don’t waste space.
cp ${output}/${prefix}_STARpass1/SJ.out.tab ${output}/${sample}SJ.pass1.out.tab
gzip ${output}/${sample}SJ.pass1.out.tab
gzip ${output}/${prefix}SJ.out.tab
gzip ${output}/${prefix}ReadsPerGene.out.tab

# Make duplicated alignments
java -Xmx16g -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${output}/${prefix}Aligned.sortedByCoord.out.bam O=${output}/${prefix}Aligned.sortedByCoord.out.md.bam M=${output}/${prefix}Aligned.sortedByCoord.out.md_metrics.txt
samtools index ${output}/${prefix}Aligned.sortedByCoord.out.md.bam

/data/sbcs/GuoLab/backup/liq17/biosoft/rnaseqc.v2.4.2.linux /data/sbcs/GuoLab/backup/liq17/ref/gencode/gencode.v26.genes.gtf ${output}/${prefix}Aligned.sortedByCoord.out.md.bam /data/sbcs/GuoLab/backup/RNA-seq/13713-XG/processed/rnaseqc_out -s ${sample_id} --mapping-quality=255 --base-mismatch=6

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for ${sample_id} STAR and RNA-SeQC"