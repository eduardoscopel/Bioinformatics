### running PROVEAN
module load BLAST+
module load CD-HIT

cd /scratch/es47540/rawdata/aneuploidy/SSD1/provean/mut_files

for file in *.txt
do
  /home/es47540/apps/provean/bin/provean.sh -q ref_prot.sl.fa -v $file --save_supporting_set ../output_files/${file%%.txt}.sss >> ../output_sp.txt
done
