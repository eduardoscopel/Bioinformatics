awk '{print $2, $1}' study.txt > files.txt # reorganize columns
split -d -l 1 files.txt strain # split file into n single files, with the strain prefix

### the following set of sed commands is exclusive for single-end read data
### I use axel to download files because it is faster
### the -o option is the output name, which, in this case, will be the strain name on EBI
sed -i "s/ ftp/.fastq.gz ftp/g" strain*
sed -i 's/^/axel -c -a -o /g' strain*

chmod +x strain*


### the following code is to create and run individual scripts to download the FASTQ files
for i in strain*
do

echo '#!/bin/bash' > scripts/down/${i}_down.sh
echo "#SBATCH --job-name=${i}_down" >> scripts/down/${i}_down.sh
echo "#SBATCH --partition=bensasson_p" >> scripts/down/${i}_down.sh
echo "#SBATCH --ntasks=1" >>scripts/down/${i}_down.sh
echo "#SBATCH --nodes=1" >>scripts/down/${i}_down.sh
echo "#SBATCH --mem=1G" >>scripts/down/${i}_down.sh
echo "#SBATCH --time=6:00:00" >>scripts/down/${i}_down.sh
echo "#SBATCH --mail-type=END,FAIL" >>scripts/down/${i}_down.sh
echo "#SBATCH --mail-user=eduardoscopel@uga.edu" >>scripts/down/${i}_down.sh
echo "#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE/output/%x_%j.out" >>scripts/down/${i}_down.sh
echo "#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/almeida/fastq/PE/error/%x_%j.err" >>scripts/down/${i}_down.sh
echo " " >> scripts/down/${i}_down.sh
echo "cd `pwd`" >> scripts/down/${i}_down.sh
echo "./$i" >> scripts/down/${i}_down.sh

chmod +x scripts/down/${i}_down.sh

sbatch scripts/down/${i}_down.sh

done
