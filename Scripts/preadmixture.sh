#!/bin/bash
#SBATCH --job-name=merge_vcf
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=30:00:00
#SBATCH --export=NONE
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=es47540@uga.edu

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/admixture

ml SAMtools/1.6-foss-2019b
ml BCFtools/1.6-foss-2019b
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
ml PLINK/1.9b_5-x86_64

input=${}

bcftools merge --file-list vcf621.txt --force-samples -o 621-merged.vcf

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
s/\<sacCer3chrXVI\>/16/g;
' 621-merged.vcf > 621-merged-temp.vcf

grep -Ev '^(sacCer3chrM)' 621-merged-temp.vcf > 621-merged.vcf

bcftools view --min-ac 1 621-merged.vcf > 621-merged-snps.vcf

vcftools --vcf 621-merged-snps.vcf --thin 1000 --recode --out 621-merged-snps

vcftools --vcf 621-merged-snps.recode.vcf  --minQ 40 --plink --out 621-merged-snps

bcftools stats 621-merged-snps.recode.vcf > 621-merged-snps.stats

plink --ped 621-merged-snps.ped -map 621-merged-snps.map --make-bed --out 621-merged-snps
