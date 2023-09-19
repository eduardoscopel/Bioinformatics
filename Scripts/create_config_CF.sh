#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N create_config.s_CF
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/
#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/

cd /scratch/es47540/rawdata/ploidy_review/Scer/

for file in *CF.bam
do
  	strain_name=${file%%.bam}
        ar+=($strain_name)
done

cd /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/config_files/

for strain in "${ar[@]}"
do

        echo "[general]" >> ${strain}_config.s.txt
        echo "chrLenFile = /scratch/es47540/rawdata/ploidy_review/Scer/ref/ref_cf/sacCer3.fa.fai" >> ${strain}_config.s.txt
        echo "ploidy = 1,2,3,4"	>> ${strain}_config.s.txt
        echo "breakPointThreshold = .8" >> ${strain}_config.s.txt
        echo "coefficientOfVariation = 0.01" >> ${strain}_config.s.txt
        echo "window = 5000" >> ${strain}_config.s.txt
        echo "step=1000" >> ${strain}_config.s.txt
        echo "chrFiles = /scratch/es47540/rawdata/ploidy_review/Scer/ref/ref_cf/" >> ${strain}_config.s.txt
        echo "minExpectedGC = 0.35" >> ${strain}_config.s.txt
        echo "maxExpectedGC = 0.55" >> ${strain}_config.s.txt
        echo "readCountTreshold = 10" >> ${strain}_config.s.txt
        echo "numberOfProceses	= 4" >> ${strain}_config.s.txt
        echo "outputDir = /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/output/" >> ${strain}_config.s.txt
        echo "contaminationAdjustment = FALSE"	>> ${strain}_config.s.txt
        echo "contamination = 0.0" >> ${strain}_config.s.txt
        echo "minMappabilityPerWindow = 0.85" >> ${strain}_config.s.txt
        echo "breakPointType = 2" >> ${strain}_config.s.txt
        echo "forceGCcontentNormalization = 0"	>> ${strain}_config.s.txt
        echo "BedGraphOutput =	FALSE" >> ${strain}_config.s.txt
        echo "telocentromeric = 5000" >> ${strain}_config.s.txt
        echo " " >> ${strain}_config.s.txt
        echo "[sample]" >> ${strain}_config.s.txt
        echo " " >> ${strain}_config.s.txt
        echo "mateFile = /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/input_files/${strain}.bam" >> ${strain}_config.s.txt
        echo "inputFormat = BAM" >> ${strain}_config.s.txt
        echo " " >> ${strain}_config.s.txt
        echo "matesOrientation =	0" >> ${strain}_config.s.txt
      	echo " " >> ${strain}_config.s.txt
      	echo "[control]" >> ${strain}_config.s.txt
      	echo "mateFile = /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/input_files/YJM1078.CF.bam" >> ${strain}_config.s.txt
      	echo "inputFormat = BAM" >> ${strain}_config.s.txt
        echo " " >> ${strain}_config.s.txt
        echo "matesOrientation =	0" >> ${strain}_config.s.txt
      	echo " " >> ${strain}_config.s.txt
      	echo "[BAF]" >> ${strain}_config.s.txt
      	echo " " >> ${strain}_config.s.txt
      	echo "makePileup = /scratch/es47540/rawdata/ploidy_review/Scer/ControlFreec/input_files/${strain}.vcf" >> ${strain}_config.s.txt
      	echo "fastaFile = /scratch/es47540/rawdata/ploidy_review/Scer/ref/ref_cf/sacCer3.fa" >> ${strain}_config.s.txt

done
