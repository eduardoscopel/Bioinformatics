
for file in *.snps.vcf
do
  sed 's/\<sacCer3chrIV\>/4/g' $file
  bgzip $file
  bcftools index ${file}.gz
  tabix -h $file 4:1043435-1052211 > ${file}.SSD1.vcf.gz

done

#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N merge_vcfs
#PBS -l walltime=48:00:00,mem=32gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -e /scratch/es47540/rawdata/aneuploidy/PCA/error/
#PBS -o /scratch/es47540/rawdata/aneuploidy/PCA/output/


cd /work/dblab/escopel/Scer/

module load BCFtools/1.6-foss-2016b

bcftools merge --force-samples -l /scratch/es47540/rawdata/aneuploidy/PCA/files_list.txt -o /scratch/es47540/rawdata/aneuploidy/PCA/vcf/nat_dip_snps.vcf


sed 's/\<sacCer3chrI\>/1/g;
s/\<sacCer3chrII\>/2/g;
s/\<sacCer3chrIII\>/3/g;
s/\<sacCer3chrIV\>/4/g;
s/\<sacCer3chrV\>/5/g;
s/\<sacCer3chrVI\>/6/g;
s/\<sacCer3chrVII\>/7/g;
s/\<sacCer3chrVIII\>/8/g;
s/\<sacCer3chrIX\>/9/g;
s/\<sacCer3chrX\>/10/g;
s/\<sacCer3chrXI\>/11/g;
s/\<sacCer3chrXII\>/12/g;
s/\<sacCer3chrXIII\>/13/g;
s/\<sacCer3chrXIV\>/14/g;
s/\<sacCer3chrXV\>/15/g;
s/\<sacCer3chrXVI\>/16/g;'
nat_dip_snps.vcf > test.vcf

grep -Ev '^(sacCer3chrM)' test.vcf > gsnpsnatdip.vcf

vcftools --vcf gsnpsnatdip.vcf --minQ 40 --plink --out natdipgenomic.snps

module load EIGENSOFT/7.2.1-foss-2016b
