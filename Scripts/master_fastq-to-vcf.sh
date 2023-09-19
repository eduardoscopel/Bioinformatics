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

cd /scratch/es47540/rawdata/aneuploidy/almeida/fastq/SE

for file in *.fastq.gz
do
        strain=${file%%.fastq.gz}

        echo "#PBS -S /bin/bash" > scripts/${strain}_map.sh
        echo "#PBS -q batch" >> scripts/${strain}_map.sh
        echo "#PBS -l walltime=24:00:00,mem=16gb,nodes=1:ppn=1" >> scripts/${strain}_map.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> scripts/${strain}_map.sh
        echo "#PBS -m ae" >> scripts/${strain}_map.sh
        echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/fastq/SE/error" >> scripts/${strain}_map.sh
        echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/fastq/SE/output" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "cd /scratch/es47540/rawdata/aneuploidy/almeida/fastq/SE" >> scripts/${strain}_map.sh
        echo "module load Java/1.8.0_144" >> scripts/${strain}_map.sh
        echo "module load FastQC/0.11.8-Java-1.8.0_144" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "mkdir fastqc.raw" >> scripts/${strain}_map.sh
        echo "mkdir fastqc.trimmed" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.raw/ ${strain}.fastq.gz" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load Trimmomatic/0.33-Java-1.8.0_144" >> scripts/${strain}_map.sh
        echo "java -jar /usr/local/apps/eb/Trimmomatic/0.33-Java-1.8.0_144/trimmomatic-0.33.jar SE ${strain}.fastq.gz ${strain}.fq ILLUMINACLIP:/scratch/es47540/rawdata/adapters/trimmomaticfasta.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.trimmed/ ${strain}.fq" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load BWA/0.7.17-foss-2016b" >> scripts/${strain}_map.sh
        echo "module load SAMtools/1.6-foss-2016b" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "bwa mem /scratch/es47540/rawdata/sacCer3/sacCer3.mfa ${strain}.fq | samtools view -b - | samtools sort -o ../../bam/${strain}.sorted.bam" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh

        chmod 755 scripts/${strain}_map.sh

        qsub scripts/${strain}_map.sh

done
