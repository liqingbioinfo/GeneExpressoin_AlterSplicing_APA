import pandas as pd
import qtl.annotation

annot = qtl.annotation.Annotation('/scratch/sbcs/chenzs/CRC_TWAS/RNA/analyze-11302020/STAR-index/100bp/gencode.v19.annotation.patched_contigs.genes.gtf')
exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                        for g in annot.genes for e in g.transcripts[0].exons],
                       columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
exon_df.to_csv('gencode.v19.GRCh37.genes.exons.txt.gz', sep='\t', index=False)
