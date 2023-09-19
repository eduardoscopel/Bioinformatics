#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N alm_vcf
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/error
#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/output

cd /scratch/es47540/rawdata/aneuploidy/almeida/bam
mkdir error
mkdir output
mkdir scripts


for file in *.bam
do
        strain=${file%%.sorted.bam}

        cd /scratch/es47540/rawdata/aneuploidy/almeida/bam/scripts

        echo "#PBS -S /bin/bash" > ${strain}_pileup.sh
        echo "#PBS -q bensasson_q" >> ${strain}_pileup.sh
        echo "#PBS -l walltime=10:00:00,mem=1gb,nodes=1:ppn=1" >> ${strain}_pileup.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_pileup.sh
        echo "#PBS -m ae" >> ${strain}_pileup.sh
        echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/bam/error" >> ${strain}_pileup.sh
        echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/bam/output" >> ${strain}_pileup.sh
        echo " " >> ${strain}_pileup.sh
        echo "cd /scratch/es47540/rawdata/aneuploidy/almeida/bam" >> ${strain}_pileup.sh
        echo " " >> ${strain}_pileup.sh
        echo "module load SAMtools/1.6-foss-2016b" >> ${strain}_pileup.sh
        echo "module load BCFtools/1.6-foss-2016b" >> ${strain}_pileup.sh
        echo " " >> ${strain}_pileup.sh
        echo "samtools mpileup -I -d 100000 -uf /scratch/es47540/rawdata/sacCer3/sacCer3.mfa ${strain}.sorted.bam | bcftools call -c > ../vcf/${strain}.vcf" >> ${strain}_pileup.sh

        qsub ${strain}_pileup.sh

done
