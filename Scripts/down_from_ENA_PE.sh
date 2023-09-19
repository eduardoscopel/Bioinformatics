awk '{print $2, $1}' study.txt > files.txt # reorganize columns
split -d -l 1 files.txt strain # split file into n single files, with the strain prefix

### the following set of sed commands is exclusive for paired-end read data
for file in strain*
  do perl -ne 'print $_ x 2' $file > dup_${file} # duplicate lines in each file
done


sed -i "1 s/;.*//" dup*
sed -i "2 s/f.*;//" dup*
sed -i "1 s/ ftp/_1.fastq.gz ftp/g" dup*
sed -i "2 s/ ftp/_2.fastq.gz ftp/g" dup*
sed -i 's/^/wget -t 0 --retry-connrefused -nc -c -O /g' dup* # for wget

chmod +x dup*


#sed -i 's/^/axel -c -a -o /g' dup* # for axel
### I use axel to download files because it is faster
### the -o option is the output name, which, in this case, will be the strain name on EBI



### the following code is to create and run individual scripts to download the FASTQ files
mkdir scripts
mkdir scripts/down/

for i in dup*
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
  echo "#SBATCH --output=output/%x_%j.out" >>scripts/down/${i}_down.sh
  echo "#SBATCH --error=error/%x_%j.err" >>scripts/down/${i}_down.sh
  echo " " >> scripts/down/${i}_down.sh
  echo "cd `pwd`" >> scripts/down/${i}_down.sh
  echo "./$i" >> scripts/down/${i}_down.sh
  chmod +x scripts/down/${i}_down.sh
  sbatch scripts/down/${i}_down.sh

done
