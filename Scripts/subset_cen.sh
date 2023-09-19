#!/bin/bash
#SBATCH --partition=bensasson_p
#SBATCH --job-name=subset_cen
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --mail-user=eduardoscopel@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/fasta


for file in *.fa
do
  strain=${file%%.fa}
  faChooseSubseq.pl -i $file -n sacCer3chrI -s "150465..152582" -o CEN1/${strain}.CEN1.fa
  faChooseSubseq.pl -i $file -n sacCer3chrII -s "237207..239323" -o CEN2/${strain}.CEN2.fa
  faChooseSubseq.pl -i $file -n sacCer3chrIII -s "113385..115501" -o CEN3/${strain}.CEN3.fa
  faChooseSubseq.pl -i $file -n sacCer3chrIV -s "448711..450821" -o CEN4/${strain}.CEN4.fa
  faChooseSubseq.pl -i $file -n sacCer3chrV -s "150987..153104" -o CEN5/${strain}.CEN5.fa
  faChooseSubseq.pl -i $file -n sacCer3chrVI -s "147510..149627" -o CEN6/${strain}.CEN6.fa
  faChooseSubseq.pl -i $file -n sacCer3chrVII -s "495920..498038" -o CEN7/${strain}.CEN7.fa
  faChooseSubseq.pl -i $file -n sacCer3chrVIII -s "104586..106703" -o CEN8/${strain}.CEN8.fa
  faChooseSubseq.pl -i $file -n sacCer3chrIX -s "354629..356745" -o CEN9/${strain}.CEN9.fa
  faChooseSubseq.pl -i $file -n sacCer3chrX -s "435307..437425" -o CEN10/${strain}.CEN10.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXI -s "439129..441246" -o CEN11/${strain}.CEN11.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXII -s "149828..151947" -o CEN12/${strain}.CEN12.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXIII -s "267031..269149" -o CEN13/${strain}.CEN13.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXIV -s "627758..629875" -o CEN14/${strain}.CEN14.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXV -s "325584..327702" -o CEN15/${strain}.CEN15.fa
  faChooseSubseq.pl -i $file -n sacCer3chrXVI -s "554957..557073" -o CEN16/${strain}.CEN16.fa
done
