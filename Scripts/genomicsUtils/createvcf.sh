#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N db20vcf
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae

cd /scratch/es47540/wkdir

module load SAMtools/1.6-foss-2016b

samtools mpileup -uf sacCer3.mfa DB_D_20.bam | bcftools call -c > DB_D_20.vcf
