#!/bin/bash
#SBATCH --job-name=exclude_inv
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu

cd /scratch/es47540/rawdata/GWAS/diploid

pts=/scratch/es47540/rawdata/GWAS/diploid/scripts
pto=/scratch/es47540/rawdata/GWAS/diploid/vcf

for file in $pto/*.vcf.gz
do

  tmp=${file##*/}
  strain=${tmp%%.vcf.gz}
  echo '#!/bin/bash' > $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --job-name=${strain}_fastq-to-vcf" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --partition=bensasson_p" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --ntasks=1" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --nodes=1" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --mem=1G" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --time=2:00:00" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --mail-type=END,FAIL" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/output/%x_%j.out" >> $pts/${strain}_exclude_inv.sh
  echo "#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/error/%x_%j.err" >> $pts/${strain}_exclude_inv.sh
  echo "cd /scratch/es47540/rawdata/GWAS/diploid/vcf" >> $pts/${strain}_exclude_inv.sh
  echo "module load BCFtools/1.6-foss-2019b" >> $pts/${strain}_exclude_inv.sh
  echo "bcftools view --min-ac 1 $file > ${strain}.snps.vcf" >> $pts/${strain}_exclude_inv.sh
  echo "bgzip ${strain}.snps.vcf" >> $pts/${strain}_exclude_inv.sh
  echo "bcftools index ${strain}.snps.vcf.gz" >> $pts/${strain}_exclude_inv.sh

  chmod 755 $pts/${strain}_exclude_inv.sh
  sbatch $pts/${strain}_exclude_inv.sh

done
bgzip -cd ZP_611.vcf.gz | awk '{ if ( ($1 =="sacCer3chrXII") && ($2 == 859486) ) {print} }' | less
sed -e 's/.*DP4=\(.*\);MQ.*/\1/' debug.sacCer3chrXII_859486.text

bcftools view --min-ac 1 diploids_g.vcf.gz -Oz -o diploids_gSNPs.vcf.gz
