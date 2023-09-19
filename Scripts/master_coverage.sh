#!/bin/bash
#SBATCH --job-name=master-depth
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu

cd /work/dblab/escopel/Scer/peter/

for file in *.bam
do
        strain=${file%%.sorted.bam}

        echo '#!/bin/bash' > /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --job-name=${strain}_bam-to-depth" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --partition=bensasson_p" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --ntasks=1" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --nodes=1" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --mem=1G" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --time=24:00:00" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --mail-type=END,FAIL" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/coverage/output/%x_%j.out" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/coverage/error/%x_%j.err" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo " " >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "cd /work/dblab/escopel/Scer/peter/" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo " " >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "module load SAMtools/1.6-foss-2019b" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo " " >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
        echo "samtools depth -a -d 100000 ${strain}.sorted.bam > /scratch/es47540/rawdata/aneuploidy/peter/coverage/${strain}.depth" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
#        echo "module load R/3.5.0-foss-2019b" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh
#        echo "Rscript --vanilla /home/es47540//scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/RD_indiv_file.R ../coverage/${strain}.depth" >> /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh

        sbatch /scratch/es47540/rawdata/aneuploidy/peter/coverage/scripts/${strain}_depth.sh

done
