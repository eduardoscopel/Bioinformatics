#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N mastermaptoref
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/output
#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/error


cd /scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE

for file in *_1.fastq.gz
do
        strain=$(echo $file | cut -d'_' -f 1)

        echo "#PBS -S /bin/bash" > scripts/${strain}_map.sh
        echo "#PBS -q batch" >> scripts/${strain}_map.sh
        echo "#PBS -l walltime=4:00:00,mem=2gb,nodes=1:ppn=1" >> scripts/${strain}_map.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> scripts/${strain}_map.sh
        echo "#PBS -m ae" >> scripts/${strain}_map.sh
        echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE/error" >> scripts/${strain}_map.sh
        echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE/output" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "cd /scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE" >> scripts/${strain}_map.sh
        echo "module load Java/1.8.0_144" >> scripts/${strain}_map.sh
        echo "module load FastQC/0.11.8-Java-1.8.0_144" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "mkdir fastqc.raw" >> scripts/${strain}_map.sh
        echo "mkdir fastqc.trimmed" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.raw/ ${strain}_1.fastq.gz" >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.raw/ ${strain}_2.fastq.gz" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load Trimmomatic/0.33-Java-1.8.0_144" >> scripts/${strain}_map.sh
        echo "java -jar /usr/local/apps/eb/Trimmomatic/0.33-Java-1.8.0_144/trimmomatic-0.33.jar PE ${strain}_1.fastq.gz ${strain}_2.fastq.gz ${strain}_1P.fq ${strain}_1U.fq ${strain}_2P.fq ${strain}_2U.fq  ILLUMINACLIP:/scratch/es47540/rawdata/adapters/trimmomaticfasta.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.trimmed/ ${strain}_1P.fq" >> scripts/${strain}_map.sh
        echo "fastqc -o fastqc.trimmed/ ${strain}_2P.fq" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "module load BWA/0.7.17-foss-2016b" >> scripts/${strain}_map.sh
        echo "module load SAMtools/1.6-foss-2016b" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh
        echo "bwa mem /scratch/es47540/rawdata/sacCer3/sacCer3.mfa ${strain}_1P.fq ${strain}_2P.fq | samtools view -b - | samtools sort -o ../../bam/${strain}.sorted.bam" >> scripts/${strain}_map.sh
        echo " " >> scripts/${strain}_map.sh

        chmod 755 scripts/${strain}_map.sh

        qsub scripts/${strain}_map.sh

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
