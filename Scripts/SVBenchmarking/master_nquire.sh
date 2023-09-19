#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N master_nQuire
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/
#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/

cd /scratch/es47540/rawdata/ploidy_review/Scer/

for file in *.bam
do
        strain_name=${file%%.bam}
        ar+=($strain_name)
done

cd /scratch/es47540/rawdata/ploidy_review/Scer/scripts/

#(head -1; tail -1) < sacCer3chr01.bed > test.bed
#tr -d '\n' < test.bed > test2.bed
#awk '{print $1,$2,$6,$7}' test2.bed > test3.bed

for strain in "${ar[@]}"
do

        echo "#PBS -S /bin/bash" > ${strain}_nQuire.sh
        echo "#PBS -q bensasson_q" >> ${strain}_nQuire.sh
        echo "#PBS -l walltime=24:00:00,mem=16gb,nodes=1:ppn=1" >> ${strain}_nQuire.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_nQuire.sh
        echo "#PBS -m ae" >> ${strain}_nQuire.sh
        echo "#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/" >> ${strain}_nQuire.sh
        echo "#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/" >> ${strain}_nQuire.sh
        echo " " >> ${strain}_nQuire.sh
	      echo "cd /scratch/es47540/rawdata/ploidy_review/Scer/" >> ${strain}_nQuire.sh
	      echo "for chromosome in sacCer3*.bed;" >> ${strain}_nQuire.sh
	      echo "do start=(\$(awk 'NR==1{print \$2}' \$chromosome))" >> ${strain}_nQuire.sh
	      echo "stop=(\$(awk 'END{print\$3}' \$chromosome))" >> ${strain}_nQuire.sh
	      echo "tmp=\${chromosome%%.bed}" >> ${strain}_nQuire.sh
	      echo "chr=(\$(echo \$tmp | cut -c 8-12))" >> ${strain}_nQuire.sh
	      echo "echo \$tmp \$start \$stop \$chr >> ${strain}.bed" >> ${strain}_nQuire.sh
	      echo "done" >> ${strain}_nQuire.sh
	      echo "/scratch/es47540/apps/nQuire/nQuire create -b ${strain}.bam -o nQuire/$strain -r ${strain}.bed" >> ${strain}_nQuire.sh

	      qsub ${strain}_nQuire.sh
done
