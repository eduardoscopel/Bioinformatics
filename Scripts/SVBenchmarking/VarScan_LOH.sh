#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N master_mpileup_VSsnp
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

        echo "#PBS -S /bin/bash" > ${strain}_pileup_VS_snp.sh
        echo "#PBS -q bensasson_q" >> ${strain}_pileup_VS_snp.sh
        echo "#PBS -l walltime=24:00:00,mem=16gb,nodes=1:ppn=1" >> ${strain}_pileup_VS_snp.sh
        echo "#PBS -M eduardoscopel@uga.edu" >> ${strain}_pileup_VS_snp.sh
        echo "#PBS -m ae" >> ${strain}_pileup_VS_snp.sh
        echo "#PBS -e /scratch/es47540/rawdata/ploidy_review/Scer/error/" >> ${strain}_pileup_VS_snp.sh
        echo "#PBS -o /scratch/es47540/rawdata/ploidy_review/Scer/output/" >> ${strain}_pileup_VS_snp.sh
        echo " " >> ${strain}_pileup_VS_snp.sh
        echo "cd /scratch/es47540/rawdata/ploidy_review/Scer/" >> ${strain}_pileup_VS_snp.sh
        echo " " >> ${strain}_pileup_VS_snp.sh
        echo "module load SAMtools/1.6-foss-2016b" >> ${strain}_pileup_VS_snp.sh
        echo "module load java" >> ${strain}_pileup_VS_snp.sh
        echo " " >> ${strain}_pileup_VS_snp.sh
        echo "samtools mpileup -f /scratch/es47540/rawdata/ploidy_review/Scer/ref/sacCer3/sacCer3.fa ${strain}.bam > ${strain}.mpileup" >> ${strain}_pileup_VS_snp.sh
        echo "java -jar /scratch/es47540/apps/VarScan.v2.3.9.jar mpileup2snp ${strain}.mpileup --min-freq-for-hom 0.85 > VarScan/${strain}.snp" >> ${strain}_pileup_VS_snp.sh
        echo "awk '{print \$1, \$2, \$3, \$4, \$5, \$7, \$8, \$9, \$10}' VarScan/${strain}.snp > VarScan/${strain}.summary.snp" >> ${strain}_pileup_VS_snp.sh

        qsub ${strain}_pileup_VS_snp.sh

done
