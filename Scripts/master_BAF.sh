#!/bin/bash
#SBATCH --job-name=master-BAF
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/duan/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/duan/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu

cd /scratch/es47540/rawdata/aneuploidy/duan/vcf

for file in *.vcf
do
  	strain=${file%%.vcf}

    echo '#!/bin/bash' > scripts/${strain}_BAF.sh
    echo "#SBATCH --job-name=${strain}_BAF" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --partition=bensasson_p" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --ntasks=1" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --nodes=1" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --mem=10G" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --time=24:00:00" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --mail-type=END,FAIL" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/duan/vcf/output/%x_%j.out" >> scripts/${strain}_BAF.sh
    echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/duan/vcf/error/%x_%j.err" >> scripts/${strain}_BAF.sh
    echo " " >> scripts/${strain}_BAF.sh
    echo "cd /scratch/es47540/rawdata/aneuploidy/duan/vcf" >> scripts/${strain}_BAF.sh
    echo " " >> scripts/${strain}_BAF.sh
    echo "module load Perl/5.26.1-GCCcore-6.4.0" >> scripts/${strain}_BAF.sh
    echo "module load R/3.5.0-foss-2019b" >> scripts/${strain}_BAF.sh
    echo " " >> scripts/${strain}_BAF.sh
    echo "vcf2alleleplot.pl -i ${strain}.vcf -o ../BAF/$strain" >> scripts/${strain}_BAF.sh
    echo "mv ${strain}.Rcmds.Rout ../BAF" >> scripts/${strain}_BAF.sh

    sbatch scripts/${strain}_BAF.sh

done
