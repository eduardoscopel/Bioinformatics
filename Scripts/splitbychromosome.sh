ar=(sacCer3chr01 sacCer3chr02 sacCer3chr03 sacCer3chr04 sacCer3chr05 sacCer3chr06 sacCer3chr07 sacCer3chr08 sacCer3chr09 sacCer3chr10 sacCer3chr11 sacCer3chr12 sacCer3chr13 sacCer3chr14 sacCer3chr15 sacCer3chr16)
for file in *.bam
do
for chr in "${ar[@]}"
do
samtools view $file $chr -b > ploidyNGS/input/${file%%.bam}${chr:7:11}.bam
done
done

ar=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR)
for file in *.bam
do
for chr in "${ar[@]}"
do
samtools view $file $chr -b > ploidyNGS/input/${file%%.bam}$chr.bam
done
done


