#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N bwa2fq
#PBS -l nodes=1:ppn=4
#PBS -l walltime=2:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae

cd /scratch/es47540/wkdir

module load SAMtools/1.6-foss-2016b

input1="DB_D_20_CGAGGCTG-AGAGTAGA_L007"
input2="DB_D_20_CGAGGCTG-AGAGTAGA_L008"
output="DB_D_20"

samtools merge $output.bam $input1.bam $input2.bam
samtools mpileup -d 100000 -uf sacCer3.mfa $output.bam | bcftools call -c | vcfutils.pl vcf2fq > $output.merged.fq
