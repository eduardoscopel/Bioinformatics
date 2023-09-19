#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N subset_SSD1
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/error/
#PBS -o /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/output/

cd /work/dblab/escopel/Scer/1011peter/vcf

pts=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/scripts
ptv=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/vcf
ptfq=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/fq
ptfa=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/fa

for file in *.vcf.gz
do

  tmp=${file##*/}
  strain=${tmp%%.vcf.gz}
  echo "#PBS -S /bin/bash" > $pts/${strain}_subset_SSD1.sh
  echo "#PBS -q bensasson_q" >> $pts/${strain}_subset_SSD1.sh
  echo "#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1" >> $pts/${strain}_subset_SSD1.sh
  echo "#PBS -M eduardoscopel@uga.edu" >> $pts/${strain}_subset_SSD1.sh
  echo "#PBS -m ae" >> $pts/${strain}_subset_SSD1.sh
  echo "#PBS -e /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/error/" >> $pts/${strain}_subset_SSD1.sh
  echo "#PBS -o /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/ssd1/output/" >> $pts/${strain}_subset_SSD1.sh
  echo "cd /work/dblab/escopel/Scer/1011peter/vcf" >> $pts/${strain}_subset_SSD1.sh
  echo "module load BCFtools/1.6-foss-2016b" >> $pts/${strain}_subset_SSD1.sh
  echo "module load seqtk/1.2-foss-2016b" >> $pts/${strain}_subset_SSD1.sh
  echo "bcftools view -r sacCer3chrIV:1045640-1049392 ${file} > $ptv/${strain}.SSD1.vcf" >> $pts/${strain}_subset_SSD1.sh
  echo "vcfutils.pl vcf2fq $ptv/${strain}.SSD1.vcf > $ptfq/${strain}.SSD1.fq" >> $pts/${strain}_subset_SSD1.sh
  echo "seqtk seq -q 40 -A $ptfq/${strain}.SSD1.fq > $ptfa/${strain}.tmp.fa" >> $pts/${strain}_subset_SSD1.sh
  echo "/home/es47540/scripts/faChooseSubseq.pl -i $ptfa/${strain}.tmp.fa -n sacCer3chrIV -s \"1045640..1049392\" \
  -r -o $ptfa/${strain}.SSD1.fa" >> $pts/${strain}_subset_SSD1.sh
  echo "sed -i 's/sacCer3chrIV:1045640-1049392/$strain/g' $ptfa/${strain}.SSD1.fa" >> $pts/${strain}_subset_SSD1.sh
  echo "command rm $ptfa/${strain}.tmp.fa" >> $pts/${strain}_subset_SSD1.sh

  qsub $pts/${strain}_subset_SSD1.sh

done

sacCer3chrI:150465-152582
