
for file in *.fa;
do
  strain=${file%%.fa};
  sed -i -e "s/sacCer3chrI\>/${strain}chr1/g;
 s/sacCer3chrII\>/${strain}chr2/g;
 s/sacCer3chrIII\>/${strain}chr3/g;
 s/sacCer3chrIV\>/${strain}chr4/g;
 s/sacCer3chrV\>/${strain}chr5/g;
 s/sacCer3chrVI\>/${strain}chr6/g;
 s/sacCer3chrVII\>/${strain}chr7/g;
 s/sacCer3chrVIII\>/${strain}chr8/g;
 s/sacCer3chrIX\>/${strain}chr9/g;
 s/sacCer3chrX\>/${strain}chr10/g;
 s/sacCer3chrXI\>/${strain}chr11/g;
 s/sacCer3chrXII\>/${strain}chr12/g;
 s/sacCer3chrXIII\>/${strain}chr13/g;
 s/sacCer3chrXIV\>/${strain}chr14/g;
 s/sacCer3chrXV\>/${strain}chr15/g;
 s/sacCer3chrXVI\>/${strain}chr16/g" $file;
done

for file in *.fa;
do
  sed -i -e "s/chrI\>/chr1/g;
 s/chrII\>/chr2/g;
 s/chrIII\>/chr3/g;
 s/chrIV\>/chr4/g;
 s/chrV\>/chr5/g;
 s/chrVI\>/chr6/g;
 s/chrVII\>/chr7/g;
 s/chrVIII\>/chr8/g;
 s/chrIX\>/chr9/g;
 s/chrX\>/chr10/g;
 s/chrXI\>/chr11/g;
 s/chrXII\>/chr12/g;
 s/chrXIII\>/chr13/g;
 s/chrXIV\>/chr14/g;
 s/chrXV\>/chr15/g;
 s/chrXVI\>/chr16/g" $file;
done



for file in *.fa;
 do
 python /home/es47540/scripts/splitfa.py $file;
done

cat *chr1.fasta > ../chr_fa/chr1.fa
cat *chr2.fasta > ../chr_fa/chr2.fa
cat *chr3.fasta > ../chr_fa/chr3.fa
cat *chr4.fasta > ../chr_fa/chr4.fa
cat *chr5.fasta > ../chr_fa/chr5.fa
cat *chr6.fasta > ../chr_fa/chr6.fa
cat *chr7.fasta > ../chr_fa/chr7.fa
cat *chr8.fasta > ../chr_fa/chr8.fa
cat *chr9.fasta > ../chr_fa/chr9.fa
cat *chr10.fasta > ../chr_fa/chr10.fa
cat *chr11.fasta > ../chr_fa/chr11.fa
cat *chr12.fasta > ../chr_fa/chr12.fa
cat *chr13.fasta > ../chr_fa/chr13.fa
cat *chr14.fasta > ../chr_fa/chr14.fa
cat *chr15.fasta > ../chr_fa/chr15.fa
cat *chr16.fasta > ../chr_fa/chr16.fa

for file in chr*.fa;
do
  fastaLC2n.pl -i $file -N -A -o temp_${file%%.fa}.fa # remove lowercase, ambiguous, and Ns
done


for file in chr*.fa;
do
  fastaLC2n.pl -i $file -N -o amb_${file%%.fa}.fa # remove lowercase, and Ns (keep ambiguous)
done



for file in temp*;
do
  alcat.pl -i $file -o ../alignment/${file:5:13};
done

alcat.pl -i 'chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa' -o wg_align.fa



for file in *.fa
do
  strain=${file%%.fa}
  sed -i "s/sacCer3chrXVI:554957-557073/$strain/g" $file
  fastaLC2n.pl -i $file -N -o amb_${file%%.fa}.fa
done

cat amb_* > ../../CEN-tree/621.CEN16.fa
