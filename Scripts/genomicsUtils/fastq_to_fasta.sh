#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N 4h00
#PBS -l walltime=24:00:00,mem=32gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -o /scratch/es47540/rawdata/contamination/scer/output
#PBS -e /scratch/es47540/rawdata/contamination/scer/error

cd /scratch/es47540/rawdata/contamination/scer/fastq

module load seqtk/1.2-foss-2016b
module load BWA/0.7.17-foss-2016b
module load SAMtools/1.6-foss-2016b
module load BCFtools/1.6-foss-2016b
module load R/3.5.0-foss-2018a-X11-20180131-GACRC

seqtk sample -s100 UCD_06-645	_1.fastq.gz 4000000 >	4h00_sample_1.fastq.gz
seqtk sample -s100 UCD_06-645	_2.fastq.gz 4000000 > 4h00_sample_2.fastq.gz
#seqtk sample -s100 CLIB219_2b_1.fastq.gz 400000 >	4h00_cont_1.fastq.gz
#seqtk sample -s100 CLIB219_2b_2.fastq.gz 400000 > 4h00_cont_2.fastq.gz

perl -p -i -e 's/ERR1309334/4h00/g' 4h00_sample_1.fastq.gz
perl -p -i -e 's/ERR1309334/4h00/g' 4h00_sample_2.fastq.gz

#perl -p -i -e 's/ERR1308618/4h00/g' 4h00_cont_1.fastq.gz
#perl -p -i -e 's/ERR1308618/4h00/g' 4h00_cont_2.fastq.gz


#cat 4h00_sample_1.fastq.gz 4h00_cont_1.fastq.gz > 4h00_merged_1.fastq.gz
#cat 4h00_sample_2.fastq.gz 4h00_cont_2.fastq.gz > 4h00_merged_2.fastq.gz

#bwa mem /scratch/es47540/rawdata/sacCer3/sacCer3.mfa 4h00_merged_1.fastq.gz 4h00_merged_2.fastq.gz | samtools view -b - | samtools sort -o ../tet_hap/4h00.sorted.bam
bwa mem /scratch/es47540/rawdata/sacCer3/sacCer3.mfa 4h00_sample_1.fastq.gz 4h00_sample_2.fastq.gz | samtools view -b - | samtools sort -o ../tet_hap/4h00.sorted.bam
cd /scratch/es47540/rawdata/contamination/scer/tet_hap
samtools mpileup -I -d 100000 -uf /scratch/es47540/rawdata/sacCer3/sacCer3.mfa 4h00.sorted.bam | bcftools call -c > 4h00.vcf
/home/es47540/scripts/vcf2alleleplot.pl -i 4h00.vcf
vcfutils.pl vcf2fq 4h00.vcf > 4h00.fq
seqtk seq -q 40 -A 4h00.fq > 4h00.fa
/home/es47540/scripts/basecomp.pl -i 4h00.fa
