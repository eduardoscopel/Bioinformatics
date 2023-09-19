#!/bin/bash
#SBATCH --job-name=merge_vcfs
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu


cd /scratch/es47540/rawdata/GWAS/diploid/vcf/snps

module load BCFtools/1.6-foss-2019b
bcftools merge -l diploids-snps-list.txt -o diploids_SNPs.vcf.gz

gunzip < diploids.vcf.gz | sed -e 's/\<sacCer3chrI\>/1/g;
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
s/\<sacCer3chrXVI\>/16/g;' | gzip -c > temp.gz

zgrep -Ev '^(sacCer3chrM)' diploids_SNPs.vcf.gz > diploids_gSNPs.vcf.gz
