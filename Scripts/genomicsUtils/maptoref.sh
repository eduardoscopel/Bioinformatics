#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N bwatest
#PBS -l nodes=1:ppn=4:dbnode
#PBS -l walltime=2:00:00
#PBS -l mem=64gb

#PBS -M eduardoscopel@uga.edu
#PBS -m ae

cd /scratch/es47540/wkdir/

module load BWA/0.7.15-foss-2016b
module load SAMtools/1.6-foss-2016b

input1="/scratch/es47540/rawdata/read_data/DB_D_20_CGAGGCTG-AGAGTAGA_L007_R1_001.fastq.gz"
input2="/scratch/es47540/rawdata/read_data/DB_D_20_CGAGGCTG-AGAGTAGA_L007_R2_001.fastq.gz"
output="DB_D_20_CGAGGCTG-AGAGTAGA_L007"

bwa mem sacCer3.mfa $input1 $input2 | samtools view -b - | samtools sort - $output

end
