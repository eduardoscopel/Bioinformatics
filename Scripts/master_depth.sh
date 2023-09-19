#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N masterdepthyue
#PBS -l walltime=1:00:00,mem=2gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae /scratch/es47540/rawdata/read_data/paradoxes/Yue/output_error

cd /scratch/es47540/rawdata/read_data/paradoxus/Yue/

for file in *_1.fastq.gz
do
  	strain_name=$(echo $file | cut -d'_' -f 1)
        ar+=($strain_name)
done

for strain in "${ar[@]}"
do
  	echo "#PBS -S /bin/bash" > ${strain}_master_depth.sh
        echo "#PBS -q bensasson_q" >> ${strain}_master_depth.sh
        echo "#PBS -l walltime=48:00:00,mem=64gb,nodes=1:ppn=1" >> ${strain}_master_depth.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_master_depth.sh
        echo "#PBS -m ae /scratch/es47540/rawdata/read_data/paradoxes/Yue/output_error" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "cd /scratch/es47540/rawdata/read_data/paradoxus/Yue/" >> ${strain}_master_depth.sh
        echo "module load Java/1.8.0_144" >> ${strain}_master_depth.sh
        echo "module load FastQC/0.11.8-Java-1.8.0_144" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "mkdir fastqc.raw" >> ${strain}_master_depth.sh
        echo "mkdir fastqc.trimmed" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "fastqc -o fastqc.raw/ ${strain}_1.fastq.gz" >> ${strain}_master_depth.sh
        echo "fastqc -o fastqc.raw/ ${strain}_2.fastq.gz" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "java -jar /usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar PE -baseout ${strain}.fq ${strain}_1.fastq ${strain}_2.fastq ILLUMINACLIP:/scratch/es47540/rawdata/read_data/paradoxus/trimmomatic/trimmomaticfasta.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36" >> \
        ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "fastqc -o fastqc.trimmed/ ${strain}_1P.fq" >> ${strain}_master_depth.sh
        echo "fastqc -o fastqc.trimmed/ ${strain}_2P.fq" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "module load BWA/0.7.15-foss-2016b" >> ${strain}_master_depth.sh
        echo "module load SAMtools/1.6-foss-2016b" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "bwa mem /scratch/es47540/rawdata/read_data/paradoxus/reference/CBS432/SparCBS432.fa ${strain}_1P.fq ${strain}_2P.fq | samtools view -b - | samtools sort - ${strain}" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "samtools depth ${strain}.bam > ${strain}.depth" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "sed -i 's/spar/${strain}/g' ${strain}.depth" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "module load R/3.4.4-foss-2016b-X11-20160819-GACRC" >> ${strain}_master_depth.sh
        echo " " >> ${strain}_master_depth.sh
        echo "R CMD BATCH /scratch/es47540/rawdata/read_data/paradoxus/Bergstrom/test/swsp.R" >> ${strain}_master_depth.sh
        echo "R CMD BATCH /scratch/es47540/rawdata/read_data/paradoxus/Bergstrom/test/mhist.R" >> ${strain}_master_depth.sh


        chmod 755 ${strain}_master_depth.sh

        qcheck=$(qsub ${strain}_master_depth.sh)

done
