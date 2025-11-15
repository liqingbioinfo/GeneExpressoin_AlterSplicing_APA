#!/usr/bin
module purge
module load GCC/6.4.0-2.28
module load STAR/2.5.4b
module load SAMtools/1.6

thread=8

batch=$1
sample=$2

############################### setting ###############################################
WORK_DIR=`pwd|perl -F"/" -lane '$p=join"/",@F[0..($#F)];print $p'`
echo "........................"
echo "run analysis in $WORK_DIR"
echo "........................"

if [[ $batch == "batch1" ]]; then
    input=/nobackup/sbcs/chenz27/Adenoma_RNA-seq/RNA-seq221214294725
else
    input=/nobackup/sbcs/chenz27/Adenoma_RNA-seq/RNA-seq221214294725-batch2
fi

gz1=`grep "$sample" $input/fq.gz.list |sed -n '1p'`
gz2=`grep "$sample" $input/fq.gz.list |sed -n '2p'`

FQ1=$input/$gz1
FQ2=$input/$gz2

# index=/nobackup/sbcs/chenz27/New_Public_data/STAR-index/hg19/150bp  # index genome using STAR 2.7.11
index=/nobackup/sbcs/chenz27/New_Public_data/STAR-index/hg19_star_2.5.4/150bp # undex genome using STAR 2.5.4 

if [[ ! -d $batch/$sample ]]; then
    mkdir $WORK_DIR/$batch
    mkdir $WORK_DIR/$batch/$sample
fi

output=$WORK_DIR/$batch/$sample
prefix="$sample."

################################ run #################################################
START=$(date +%s)

STAR  --runMode alignReads \
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

samtools sort --threads ${thread} ${output}/${prefix}Aligned.out.bam  -o ${output}/${prefix}Aligned.sortedByCoord.out.bam
samtools index ${output}/${prefix}Aligned.sortedByCoord.out.bam

rm ${output}/${prefix}Aligned.out.bam

cp ${output}/${prefix}_STARpass1/SJ.out.tab ${output}/${sample_id}.SJ.pass1.out.tab
gzip ${output}/${sample_id}.SJ.pass1.out.tab
gzip ${output}/${prefix}SJ.out.tab
gzip ${output}/${prefix}ReadsPerGene.out.tab

# mark duplicate
module load picard/2.18.27

java -Xmx16g -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${output}/${prefix}Aligned.sortedByCoord.out.bam O=${output}/${prefix}Aligned.sortedByCoord.out.md.bam M=${output}/${prefix}Aligned.sortedByCoord.out.md_metrics.txt
samtools index ${output}/${prefix}Aligned.sortedByCoord.out.md.bam

#rm ${output}/${prefix}Aligned.sortedByCoord.out.bam
#rm ${output}/${prefix}Aligned.sortedByCoord.out.bam.bai

# quantification
/nobackup/sbcs/scratch/chenzs/software/rnaseqc.v2.4.2.linux ${index}/gencode.v19.annotation.patched_contigs.genes.gtf ${output}/${prefix}Aligned.sortedByCoord.out.md.bam ${output} -s ${sample_id} --mapping-quality=255 --base-mismatch=6

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for ${sample_id} STAR and RNA-SeQC"
