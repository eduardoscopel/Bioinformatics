#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N master_IQTree
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -l mem=1gb
#PBS -M es47540@uga.edu
#PBS -m ae

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/dipnonhetna/scripts

for i in {1..100}
do
  echo "#PBS -S /bin/bash" > IQTree$i.sh
  echo "#PBS -q highmem_q" >> IQTree$i.sh
  echo "#PBS -N IQTree$i" >> IQTree$i.sh
  echo "#PBS -l walltime=120:00:00,mem=72gb,nodes=1:ppn=9" >> IQTree$i.sh
  echo "#PBS -M eduardoscopel@uga.edu" >> IQTree$i.sh
  echo "#PBS -m ae" >> IQTree$i.sh
  echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/dipnonhetna/error/" >> IQTree$i.sh
  echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/dipnonhetna/output/" >> IQTree$i.sh
  echo "cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/dipnonhetna/alignment" >> IQTree$i.sh
  echo "module load IQ-TREE/1.6.5-omp" >> IQTree$i.sh
  echo "iqtree -nt 9 -s wg_align.fa -m GTR+G -seed $i -bo 1 -pre boot$i" >> IQTree$i.sh

  chmod 755 IQTree$i.sh

  qsub IQTree$i.sh

done
