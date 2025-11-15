#parallel -q echo sed -n \''1p;2p;3p;9p;10p;15p;24p;52p'\' {} \|paste -s -d \''\t'\' \| perl -pne \''s/\|/\t/g'\' \| sed \''s/ //g'\' \> {}.qc ::: $(cat metrics.list) |parallel -j 8 bash -c
cat batch1/*/*.metrics.tsv.qc > batch1_qc.txt
cat batch2/*/*.metrics.tsv.qc > batch2_qc.txt

# umr: uniquely mapped reads
perl -F"\t" -lane 'if(($F[3] >= 0.2) & ($F[11] <= 0.4) & ($F[13] <=0.3) & ($F[15] > 10000000)){print $_}' batch1_qc.txt  > batch1_qc.txt.pass
perl -F"\t" -lane 'if(($F[3] >= 0.2) & ($F[11] <= 0.4) & ($F[13] <=0.3) & ($F[15] > 10000000)){print $_}' batch2_qc.txt  > batch2_qc.txt.pass

perl -F"\t" -lane 'if(($F[3] >= 0.2) & ($F[11] <= 0.4) & ($F[13] <=0.3) & ($F[15] > 10000000)){}else{print $_}' batch1_qc.txt > batch1_qc.txt.unpass
perl -F"\t" -lane 'if(($F[3] >= 0.2) & ($F[11] <= 0.4) & ($F[13] <=0.3) & ($F[15] > 10000000)){}else{print $_}' batch2_qc.txt > batch2_qc.txt.unpass
