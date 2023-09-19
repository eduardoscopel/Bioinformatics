#!/bin/bash
#SBATCH --job-name=master-fastq-to-vcf
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/hahn2021/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/hahn2021/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu


cd /scratch/es47540/rawdata/aneuploidy/hahn2021/fastq

for file in *_1.fastq.gz
do
        strain=${file%%_1.fastq.gz}
        cd /scratch/es47540/rawdata/aneuploidy/hahn2021

        echo '#!/bin/bash' > scripts/${strain}_map.sh
        echo "#SBATCH --job-name=${strain}_fastq-to-vcf" >> scripts/${strain}_map.sh
        echo "#SBATCH --partition=bensasson_p" >> scripts/${strain}_map.sh
        echo "#SBATCH --ntasks=1" >> scripts/${strain}_map.sh
        echo "#SBATCH --nodes=1" >> scripts/${strain}_map.sh
        echo "#SBATCH --mem=8G" >> scripts/${strain}_map.sh
        echo "#SBATCH --time=24:00:00" >> scripts/${strain}_map.sh
        echo "#SBATCH --mail-type=END,FAIL" >> scripts/${strain}_map.sh
        echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> scripts/${strain}_map.sh
        echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/hahn2021/output/%x_%j.out" >> scripts/${strain}_map.sh
        echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/hahn2021/error/%x_%j.err" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "cd /scratch/es47540/rawdata/aneuploidy/hahn2021/fastq" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load BWA/0.7.17-GCC-8.3.0" >> scripts/${strain}_map.sh
        echo "module load SAMtools/1.6-foss-2019b" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "bwa mem /scratch/es47540/rawdata/sacCer3/sacCer3.mfa ${strain}_1.fastq.gz ${strain}_2.fastq.gz | samtools view -b - | samtools sort -o ../bam/${strain}.sorted.bam" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load BCFtools/1.6-foss-2019b" >> scripts/${strain}_map.sh
        echo "samtools mpileup -I -d 100000 -uf /home/es47540/sacCer3/sacCer3.mfa ../bam/${strain}.sorted.bam | bcftools call -c | bgzip -c > ../vcf/${strain}.vcf.gz"  >> scripts/${strain}_map.sh
        echo "gunzip -c $file | vcfutils.pl vcf2fq - | seqtk seq -q 40 -A > fasta/${strain}.fa" >> scripts/${strain}_vcf2fa.sh
        echo "bcftools index ../vcf/${strain}.vcf.gz"  >> scripts/${strain}_map.sh

        chmod 755 scripts/${strain}_map.sh
        sbatch scripts/${strain}_map.sh

done
