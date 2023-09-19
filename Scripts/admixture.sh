#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/scripts


for K in {1..20}
do
  echo "#!/bin/bash" > scripts/K${K}_12345.sh
  echo "#SBATCH --job-name=${K}_12345" >> scripts/K${K}_12345.sh
  echo "#SBATCH --partition=bensasson_p" >> scripts/K${K}_12345.sh
  echo "#SBATCH --ntasks=1" >> scripts/K${K}_12345.sh
  echo "#SBATCH --nodes=1" >> scripts/K${K}_12345.sh
  echo "#SBATCH --mem=1G" >> scripts/K${K}_12345.sh
  echo "#SBATCH --time=2:00:00" >> scripts/K${K}_12345.sh
  echo "#SBATCH --mail-type=END,FAIL" >> scripts/K${K}_12345.sh
  echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> scripts/K${K}_12345.sh
  echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out" >> scripts/K${K}_12345.sh
  echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err"
  echo " " >> scripts/K${K}_12345.sh
  echo "cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/admixture" >> scripts/K${K}_12345.sh
  echo "ml ADMIXTURE/1.3.0" >> scripts/K${K}_12345.sh
  echo "admixture --cv -s 12345 621-merged-snps.bed ${K}" >> scripts/K${K}_12345.sh
  chmod 755 scripts/K${K}_12345.sh
  sbatch scripts/K${K}_12345.sh
done

for file in *.out
do
seed=${file%_*_*}; seed=${seed#*_}
 cat $file | grep 'K=\|^Loglikelihood' | column | awk '{print $5, $6, $2}' >> ../$seed.txt
done

sed -i 's/(K=//g' *.txt
sed -i 's/)://g' *.txt


for K in {1..20}
do
  echo "#!/bin/bash" > scripts/K${K}_9999.sh
  echo "#SBATCH --job-name=${K}_1" >> scripts/K${K}_9999.sh
  echo "#SBATCH --partition=bensasson_p" >> scripts/K${K}_9999.sh
  echo "#SBATCH --ntasks=1" >> scripts/K${K}_9999.sh
  echo "#SBATCH --nodes=1" >> scripts/K${K}_9999.sh
  echo "#SBATCH --mem=1G" >> scripts/K${K}_9999.sh
  echo "#SBATCH --time=24:00:00" >> scripts/K${K}_9999.sh
  echo "#SBATCH --mail-type=END,FAIL" >> scripts/K${K}_9999.sh
  echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >> scripts/K${K}_9999.sh
  echo "#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/admixture/output/%x_%j.out" >> scripts/K${K}_9999.sh
  echo "#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/admixture/error/%x_%j.err" >> scripts/K${K}_9999.sh
  echo " " >> scripts/K${K}_9999.sh
  echo "cd /scratch/es47540/rawdata/GWAS/diploid/admixture" >> scripts/K${K}_9999.sh
  echo "ml ADMIXTURE/1.3.0" >> scripts/K${K}_9999.sh
  echo "admixture --cv -s 9999 gain2n_gSNPs_q40.bed s9999/${K}" >> scripts/K${K}_9999.sh
  chmod 755 scripts/K${K}_9999.sh
  sbatch scripts/K${K}_9999.sh
done
