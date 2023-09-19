qsub -I -l walltime=24:00:00 -l nodes=1:ppn=1 #instead of qlogin
module load seqtk/1.2-foss-2016b
seqtk seq -q 40 -A SRR850113.fq > SRR850113.fa
