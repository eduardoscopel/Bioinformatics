
#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N chgsname
#PBS -l nodes=1:ppn=4
#PBS -l walltime=2:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae

cd /scratch/es47540/rawdata/zhu

for file in *.coverage
do

STRNAME=$(echo $file | cut -d'.' -f 1)

sed -i "s/sacCer3/$STRNAME/g" $file

done
