
 # reference: https://davidaknowles.github.io/leafcutter/articles/Usage.html
   
 # bam files
 ls /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/batch1/*/*.Aligned.sortedByCoord.out.md.bam |awk -F"/" '{print $7"\t"$8"\t"$0}' > batch1_aligned_bam.list 
 ls /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/batch2/*/*.Aligned.sortedByCoord.out.md.bam |awk -F"/" '{print $7"\t"$8"\t"$0}' > batch2_aligned_bam.list 

 # Converting bams to juncs
 parallel --colsep "\t" -q echo regtools junctions extract -a 8 -m 50 -M 500000 -s 0 {3} \| gzip -c \> {1}_junc/{2}.regtools_junc.txt.gz  :::: batch1_batch2_aligned_bam.list|parallel -j 8 bash -c
 ls /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/splice/batch1_junc/*.junc |perl -F"\t" -lane 'if(/(L01-Ts-\d+A)/){print "batch1\t$1\t$_"}' > batch1_junc.list
 ls /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/splice/batch2_junc/*.junc |perl -F"\t" -lane 'if(/(L01-Ts-\d+A)/){print "batch2\t$1\t$_"}' > batch2_junc.list

 # juncs files for 190 included samples
 pip3 install qtl
 compNcol.hash.pl -c 1,2 -q ../expression/batch1_batch2_qc.txt.used.samples -d 1,2 -db batch1_junc.list -e |grep -v nocommon |cut -f 5  > batch1_batch2_qc.txt.used.samples.junc.list
 compNcol.hash.pl -c 1,2 -q ../expression/batch1_batch2_qc.txt.used.samples -d 1,2 -db batch2_junc.list -e |grep -v nocommon |cut -f 5  >> batch1_batch2_qc.txt.used.samples.junc.list
 
 # run in /nobackup/sbcs/chenz27/anaconda3
 # Generating intron excision ratios with LeafCutter

 python /nobackup/sbcs/chenz27/soft/leafcutter/clustering/leafcutter_cluster_regtools.py -j batch1_batch2_qc.txt.used.samples.junc.list -m 30 -l 500000 -o CRC_adenoma_190samples_Splice

 python3 cluster_prepare_fastqtl.py batch1_batch2_qc.txt.used.samples.junc.list  gencode.v19.GRCh37.genes.exons.txt.gz gencode.v19.annotation.patched_contigs.genes.gtf CRC_adenoma_190samples_Splice batch1_batch2_qc.txt.used.samples.list  

