#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N almeida_RD
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/error/
#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/output/

cd /scratch/es47540/rawdata/aneuploidy/almeida/coverage

mkdir error
mkdir output
mkdir scripts

for file in *.depth
do
  	strain=${file%%.depth}

        echo "#PBS -S /bin/bash" > scripts/${strain}_RD.sh
        echo "#PBS -q bensasson_q" >> scripts/${strain}_RD.sh
        echo "#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1" >> scripts/${strain}_RD.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> scripts/${strain}_RD.sh
        echo "#PBS -m ae" >> scripts/${strain}_RD.sh
        echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/almeida/coverage/error/" >> scripts/${strain}_RD.sh
        echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/almeida/coverage/output/" >> scripts/${strain}_RD.sh
        echo " " >> scripts/${strain}_RD.sh
        echo "cd /scratch/es47540/rawdata/aneuploidy/almeida/coverage" >> scripts/${strain}_RD.sh
        echo " " >> scripts/${strain}_RD.sh
        echo "module load R/3.5.0-foss-2019b" >> scripts/${strain}_RD.sh
        echo " " >> scripts/${strain}_RD.sh
        echo "Rscript --vanilla /scratch/es47540/scripts/plot_RD_DT.R ${strain}.depth" >> scripts/${strain}_RD.sh

        qsub scripts/${strain}_RD.sh

done
