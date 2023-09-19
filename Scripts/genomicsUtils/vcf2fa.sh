#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=vcf2fasta
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --mail-user=eduardoscopel@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9

cat alignment/vcf621.txt | while read line;
do
file="/work/dblab/escopel/Scer/1011peter/vcf/"$line
strain=${line%%.vcf.gz}

  echo "#!/bin/bash" > scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --partition=batch" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --time=1:00:00" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --mem=1G" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --nodes=1" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --ntasks=1" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --export=NONE" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --mail-type=END,FAIL" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err" >> scripts/${strain}_vcf2fa.sh
  echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out" >> scripts/${strain}_vcf2fa.sh
  echo "" >> scripts/${strain}_vcf2fa.sh
  echo "cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9" >> scripts/${strain}_vcf2fa.sh
  echo "module load seqtk/1.2-foss-2019b" >> scripts/${strain}_vcf2fa.sh
  echo "module load BCFtools/1.6-foss-2019b" >> scripts/${strain}_vcf2fa.sh
  echo "" >> scripts/${strain}_vcf2fa.sh
  echo "gunzip -c $file | vcfutils.pl vcf2fq - | seqtk seq -q 40 -A > fasta/${strain}.fa" >> scripts/${strain}_vcf2fa.sh


  chmod 755 scripts/${strain}_vcf2fa.sh
  sbatch scripts/${strain}_vcf2fa.sh

done
