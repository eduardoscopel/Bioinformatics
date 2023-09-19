#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N master_CF
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/
#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/

cd /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/config_files/nocontrol/

for file in *.txt
do
  	strain_name=${file%%.txt}
        ar+=($strain_name)
done

cd /scratch/es47540/rawdata/ploidy_review/Scer/scripts/

for strain in "${ar[@]}"
do

        echo "#PBS -S /bin/bash" > ${strain}_RunCF.sh
        echo "#PBS -q bensasson_q" >> ${strain}_RunCF.sh
        echo "#PBS -l walltime=24:00:00,mem=16gb,nodes=1:ppn=1" >> ${strain}_RunCF.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_RunCF.sh
        echo "#PBS -m ae" >> ${strain}_RunCF.sh
        echo "#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/" >> ${strain}_RunCF.sh
        echo "#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/" >> ${strain}_RunCF.sh
        echo " " >> ${strain}_RunCF.sh
        echo "module load SAMtools/1.6-foss-2016b" >> ${strain}_RunCF.sh
        echo "module load BEDTools/2.26.0-foss-2016b" >> ${strain}_RunCF.sh
        echo "module load R/3.4.4-foss-2016b-X11-20160819-GACRC" >> ${strain}_RunCF.sh
        echo " " >> ${strain}_RunCF.sh
        echo "/scratch/es47540/apps/ControlFreec/FREEC-11.4/src/freec -config /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/config_files/nocontrol/${strain}.txt" >> $$

        qsub ${strain}_RunCF.sh

done
