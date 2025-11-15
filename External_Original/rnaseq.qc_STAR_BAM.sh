#!/usr/bin
#RNA-seq quality control and metrics, typically on STAR-aligned BAM files. 
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

# index=/nobackup/sbcs/chenz27/New_Public_data/STAR-index/hg19/150bp  # index genome using STAR 2.7.11
index=/nobackup/sbcs/chenz27/New_Public_data/STAR-index/hg19_star_2.5.4/150bp # undex genome using STAR 2.5.4

output=$WORK_DIR/$batch/$sample

# quantification
/nobackup/sbcs/scratch/chenzs/software/rnaseqc.v2.4.2.linux ${index}/gencode.v19.annotation.patched_contigs.genes.gtf ${output}/${sample}.Aligned.sortedByCoord.out.md.bam ${output} -s ${sample} --mapping-quality=255 --base-mismatch=6

