#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N chgchrnumber
#PBS -l nodes=1:ppn=4
#PBS -l walltime=6:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae


cd /scratch/es47540/rawdata/aneuploidy/

for file in *.fa
do
  strain=${file%%.fa}
  perl -pi -e 's/sacCer3chrI\>/${strain}chr1/g;
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
   s/sacCer3chrXVI\>/${strain}chr16/g' $file

done


for file in *.fa; do   strain=${file%%.fa};
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
 s/sacCer3chrXVI\>/${strain}chr16/g" $file;  done
