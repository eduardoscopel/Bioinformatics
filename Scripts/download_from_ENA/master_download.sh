#PBS -S /bin/bash
#PBS -q batch
#PBS -N downloadstrainsfromENA
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1:dbnode
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -o /path/to/dir/output
#PBS -e /path/to/dir/error

mkdir /path/to/dir/output
mkdir /path/to/dir/error

cd /path/to/dir/

for i in do_get_*
do

echo "#PBS -S /bin/bash" > ${i}_down.sh
echo "#PBS -q batch" >> ${i}_down.sh
echo "#PBS -l walltime=24:00:00,mem=1gb,nodes=1:ppn=1:dbnode" >> ${i}_down.sh
echo "#PBS -M eduardoscopel@uga.edu" >> ${i}_down.sh
echo "#PBS -m ae" >> ${i}_down.sh
echo "#PBS -o /path/to/dir/output" >> ${i}_down.sh
echo "#PBS -e /path/to/dir/error" >> ${i}_down.sh
echo " " >> ${i}_down.sh
echo "cd /path/to/dir/" >> ${i}_down.sh

echo "./$i" >> ${i}_down.sh

chmod 755 ${i}_down.sh

qcheck=$(qsub ${i}_down.sh)

done


