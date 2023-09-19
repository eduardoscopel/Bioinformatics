#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N borneman
#PBS -l nodes=1:ppn=4
#PBS -l walltime=6:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae

cd /scratch/es47540/rawdata/aneuploidy/borneman

module load SAMtools/1.6-foss-2016b

# create a file with the number of reads at every position of chromosomes

for bamfile in *bam

do

samtools depth "$bamfile" > "$bamfile".coverage

done

# replace sacCer3 in the name of the chromosomes by the name of the strain

for covfile in *.coverage

do

STRNAME=$(echo $covfile | cut -d'.' -f 1)

sed -i "s/sacCer3/$STRNAME/g" $covfile

done
