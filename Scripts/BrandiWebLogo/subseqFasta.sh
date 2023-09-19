
cat cyp51*_aln.fasta > cyp51.fasta
mafft.bat cyp51.fasta > cyp51_aligned.fasta
run_pasta.py -d protein -i cyp51_aligned.fasta

sed -i '' 's/ //g' cyp51_final_align.fasta
cat cyp51_final_align.fasta | grep '>' | sed 's/^.*>//g' > cyp51_final_align.IDs
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '107..125' -o subseq/cyp51_final_align.SRS1.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '206..213' -o subseq/cyp51_final_align.SRS2.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '228..238' -o subseq/cyp51_final_align.SRS3.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '282..302' -o subseq/cyp51_final_align.SRS4.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '293..298' -o subseq/cyp51_final_align.AGXDTT.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '356..359' -o subseq/cyp51_final_align.EXXR.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '362..371' -o subseq/cyp51_final_align.SRS5.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '412..414' -o subseq/cyp51_final_align.PER.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '447..456' -o subseq/cyp51_final_align.FXXGXXXCXG.fasta
faChooseSubseq.pl -i cyp51_final_align.fasta -n cyp51_final_align.IDs -s '491..499' -o subseq/cyp51_final_align.SRS6.fasta





