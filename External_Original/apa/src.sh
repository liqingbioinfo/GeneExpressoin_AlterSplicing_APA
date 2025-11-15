 ##
 ## call APA using DaPars
 ## reference: https://github.com/3UTR/DaPars2
 ## 

 # bam files
 cat ../splice/batch1_aligned_bam.list ../splice/batch2_aligned_bam.list  > batch1_batch2_aligned_bam.list

 ## 
 ## get 3UTR annotation for reference genome from UCSC tables
 ##
 
 python /nobackup/sbcs/chenz27/soft/DaPars2-2.1/src/DaPars_Extract_Anno.py -b Gencode_19_GeneAnnotation.bed -s Gencode_19_IDmapping.bed -o Gencode_19_3UTR_annotation.bed


 ##
 ## convert bam files to wig
 ## 
 
 perl run_apa.pl /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/expression/batch1_batch2_aligned_bam.list.used 


 ##
 ## run DaPars
 ##
 
 sed 's/\.Aligned\.sortedByCoord\.out\.md\.bam//g' ../expression/batch1_batch2_aligned_bam.list.used |cat | parallel --colsep "\t" -q echo grep \'"^Mapped Reads\b"\' {3}.metrics.tsv \|awk \''{print"wig/{2}.wig""\t"$3}'\' |bash  > mapping_wig_location_with_depth.txt

 ls /nobackup/sbcs/chenz27/Adenoma_RNA-seq/processed/apa/wig/*.wig |cat | tr '\n' ','|sed 's/,$/\n/'  >> asian_CRC_Adenoma_APA_configure.txt 

 cat chrList.txt |parallel -q echo python /nobackup/sbcs/chenz27/soft/DaPars2-2.1/src/Dapars2_Multi_Sample.py asian_Adenoma_APA_configure.txt {} |parallel -j 8 bash -c 

 cat output_*/Adenoma_APA_Asian_result_temp.*.txt > CRC_adenoma_190samples_APA.txt
