#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N master_VarScan
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

for strain in "${ar[@]}"
do

        echo "#PBS -S /bin/bash" > ${strain}_VarScan.sh
        echo "#PBS -q bensasson_q" >> ${strain}_VarScan.sh
        echo "#PBS -l walltime=24:00:00,mem=16gb,nodes=1:ppn=1" >> ${strain}_VarScan.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_VarScan.sh
        echo "#PBS -m ae" >> ${strain}_VarScan.sh
        echo "#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/" >> ${strain}_VarScan.sh
        echo "#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "cd /scratch/es47540/rawdata/ploidy_review/Scer/" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "module load SAMtools/1.6-foss-2016b" >> ${strain}_VarScan.sh
        echo "module load Java/1.8.0_144" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "mapped_reads_p_normal=( \$(awk 'FNR == 7 {print \$3}' peter_normal_stats.txt) )" >> ${strain}_VarScan.sh
        echo "average_length_p_normal=( \$(awk 'FNR == 23 {print \$3}' peter_normal_stats.txt) )" >> ${strain}_VarScan.sh
        echo "uq_mapped_reads_p_normal=\$((\$mapped_reads_p_normal * \$average_length_p_normal))" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "mapped_reads_s_normal=( \$(awk 'FNR == 7 {print \$3}' strope_normal_stats.txt) )" >> ${strain}_VarScan.sh
        echo "average_length_s_normal=( \$(awk 'FNR == 23 {print \$3}' strope_normal_stats.txt) )" >> ${strain}_VarScan.sh
        echo "uq_mapped_reads_s_normal=\$((\$mapped_reads_s_normal * \$average_length_s_normal))" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "samtools stats ${strain}.bam | grep '^SN' | cut -f 2- > ${strain}_stats.txt" >> ${strain}_VarScan.sh
        echo "mapped_reads_tumor=( \$(awk 'FNR == 7 {print \$3}' ${strain}_stats.txt) )" >> ${strain}_VarScan.sh
        echo "average_length_tumor=( \$(awk 'FNR == 23 {print \$3}' ${strain}_stats.txt) )" >> ${strain}_VarScan.sh
        echo "uq_mapped_reads_tumor=\$((\$mapped_reads_tumor * \$average_length_tumor))" >> ${strain}_VarScan.sh
        echo "data_ratio_p=\`echo \$uq_mapped_reads_p_normal / \$uq_mapped_reads_tumor | bc -l\`" >> ${strain}_VarScan.sh
        echo "data_ratio_s=\`echo \$uq_mapped_reads_s_normal / \$uq_mapped_reads_tumor | bc -l\`" >> ${strain}_VarScan.sh
        echo " " >> ${strain}_VarScan.sh
        echo "cd /scratch/es47540/rawdata/ploidy_review/Scer/copynumber/" >> ${strain}_VarScan.sh
        echo "samtools mpileup -q 1 -f /scratch/es47540/rawdata/ploidy_review/Scer/ref/sacCer3.fa /scratch/es47540/rawdata/ploidy_review/Scer/DBVPG1554.bam /scratch/es47540/rawdata/ploidy_review/Scer/${strain}.bam | java -jar /scratch/es47540/apps/VarScan.v2.3.9.jar copynumber --mpileup ${strain}.p --data-ratio \$data_ratio_p" >> ${strain}_VarScan.sh
        echo "samtools mpileup -q 1 -f /scratch/es47540/rawdata/ploidy_review/Scer/ref/sacCer3.fa /scratch/es47540/rawdata/ploidy_review/Scer/YJM1078.bam /scratch/es47540/rawdata/ploidy_review/Scer/${strain}.bam | java -jar /scratch/es47540/apps/VarScan.v2.3.9.jar copynumber --mpileup ${strain}.s --data-ratio \$data_ratio_s" >> ${strain}_VarScan.sh
        echo "sed -i 's/sacCer3chr//g' ${strain}.s.copynumber" >> ${strain}_VarScan.sh
        echo "sed -i 's/sacCer3chr//g' ${strain}.p.copynumber" >> ${strain}_VarScan.sh
        echo "java -jar /scratch/es47540/apps/VarScan.v2.3.9.jar copyCaller ${strain}.p.copynumber --output-file ${strain}.p.copynumber.called" >> ${strain}_VarScan.sh
        echo "java -jar /scratch/es47540/apps/VarScan.v2.3.9.jar copyCaller ${strain}.s.copynumber --output-file ${strain}.s.copynumber.called" >> ${strain}_VarScan.sh

        qsub ${strain}_VarScan.sh

done
