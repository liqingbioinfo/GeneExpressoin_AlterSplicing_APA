data=/home/ec2-user/project1/chenzs/UVA_CRC_RNA-seq/rnaseq

module load star/2.7.11a

##To generate STAR index
/nobackup/sbcs/chenz27/soft/STAR-2.7.11a/source/STAR --runMode genomeGenerate \
     --genomeDir ./ \
     --genomeFastaFiles ./Homo_sapiens_assembly19.fasta \
     --sjdbGTFfile ./gencode.v19.annotation.patched_contigs.gtf  \
     --sjdbOverhang 149 \
     --runThreadN 8

##Use STAR to map RNA seq from fastq to BAM
thread=16
sample_id=$1
fq_folder=${data}
index_folder=/home/ec2-user/project1/chenzs/public_data/STAR-index/hg38/100bp
mkdir /home/ec2-user/project1/chenzs/crctwas/UVA_CRC_RNA-seq/processed/${sample_id}
output=/home/ec2-user/project1/chenzs/crctwas/UVA_CRC_RNA-seq/processed/${sample_id}
prefix="${sample_id}."

START=$(date +%s)

STAR --runMode alignReads \
     --runThreadN ${thread} \
     --genomeDir ${index_folder} \
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
     --readFilesIn ${fq_folder}/${sample_id}_R1.fq.gz ${fq_folder}/${sample_id}_R2.fq.gz \
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
cp ${output}/${prefix}_STARpass1/SJ.out.tab ${output}/${sample_id}.SJ.pass1.out.tab
gzip ${output}/${sample_id}.SJ.pass1.out.tab
gzip ${output}/${prefix}SJ.out.tab
gzip ${output}/${prefix}ReadsPerGene.out.tab

picard MarkDuplicates I=${output}/${prefix}Aligned.sortedByCoord.out.bam O=${output}/${prefix}Aligned.sortedByCoord.out.md.bam M=${output}/${prefix}Aligned.sortedByCoord.out.md_metrics.txt

/home/ec2-user/project1/chenzs/soft/rnaseqc.v2.4.2.linux ${index_folder}/gencode.v26.annotation.gene.gtf ${output}/${prefix}Aligned.sortedByCoord.out.md.bam ${output} -s ${sample_id} --mapping-quality 255 --base-mismatch 6

END=$(date +%s)
DIFF=$(( $END - $START  ))
echo "It took $DIFF seconds for ${sample_id} to align and count"

