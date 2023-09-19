#!/bin/bash
#SBATCH --job-name=contSeqs
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu

cd $SLURM_SUBMIT_DIR
cd fastq

module load seqtk/1.2-foss-2019b
module load BWA/0.7.17-GCC-8.3.0
module load SAMtools/1.6-foss-2019b
module load BCFtools/1.6-foss-2019b
module load R/3.5.0-foss-2019b

echo 'Contaminated strain name (use prefix of _x.fastq.gz):' $1
echo 'Contaminant strain name (use prefix of _x.fastq.gz):' $2
echo 'Contamination rate (use 0-1 interval):' $3
sampleReads=$(echo "x = 4000000 - 4000000*$3; x / 1" | bc)
contReads=$(echo "x = 4000000*$3; x / 1" | bc)
proportion=$(echo "x = 100*$3; x / 1" | bc)
echo 'Number of reads sampled from contaminant seq:' $contReads
echo 'Number of reads sampled from contaminated seq:' $sampleReads

seqtk sample -s100 $1_1.fastq.gz $sampleReads > $1_$2_${proportion}_sample_1.fastq.gz
seqtk sample -s100 $1_2.fastq.gz $sampleReads > $1_$2_${proportion}_sample_2.fastq.gz
seqtk sample -s100 $2_1.fastq.gz $contReads > $1_$2_${proportion}_cont_1.fastq.gz
seqtk sample -s100 $2_2.fastq.gz $contReads > $1_$2_${proportion}_cont_2.fastq.gz

sampleID=$(zcat $1_1.fastq.gz | head -n 1)
sampleID=${sampleID##*@}
sampleID=${sampleID%%.1*}

perl -p -i -e "s/$sampleID/$1_$2_${proportion}/g" $1_$2_${proportion}_sample_1.fastq.gz
perl -p -i -e "s/$sampleID/$1_$2_${proportion}/g" $1_$2_${proportion}_sample_2.fastq.gz

contID=$(zcat $2_1.fastq.gz | head -n 1)
contID=${contID##*@}
contID=${contID%%.1*}

perl -p -i -e "s/$contID/$1_$2_${proportion}/g" $1_$2_${proportion}_cont_1.fastq.gz
perl -p -i -e "s/$contID/$1_$2_${proportion}/g" $1_$2_${proportion}_cont_2.fastq.gz

cat $1_$2_${proportion}_sample_1.fastq.gz $1_$2_${proportion}_cont_1.fastq.gz > $1_$2_${proportion}_merged_1.fastq.gz
cat $1_$2_${proportion}_sample_2.fastq.gz $1_$2_${proportion}_cont_2.fastq.gz > $1_$2_${proportion}_merged_2.fastq.gz

bwa mem /home/es47540/sacCer3/sacCer3.mfa $1_$2_${proportion}_merged_1.fastq.gz $1_$2_${proportion}_merged_2.fastq.gz | samtools view -b - | samtools sort -o ../$1_$2_${proportion}.sorted.bam

samtools mpileup -I -d 100000 -uf /home/es47540/sacCer3/sacCer3.mfa ../$1_$2_${proportion}.sorted.bam | bcftools call -c > ../$1_$2_${proportion}.vcf

cd ../

sed -i "s/sacCer3/${1}_${2}_${proportion}/g" $1_$2_${proportion}.vcf

/home/es47540/scripts/vcf2alleleplot.pl -i $1_$2_${proportion}.vcf

vcfutils.pl vcf2fq $1_$2_${proportion}.vcf > $1_$2_${proportion}.fq

seqtk seq -q 40 -A $1_$2_${proportion}.fq > $1_$2_${proportion}.fa

/home/es47540/scripts/basecomp.pl -s -i $1_$2_${proportion}.fa > $1_$2_${proportion}.bc