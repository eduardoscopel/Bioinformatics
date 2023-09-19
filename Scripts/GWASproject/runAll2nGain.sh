#!/bin/bash
#SBATCH --job-name=runAll
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/gainOnly/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/gainOnly/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu


cd /scratch/es47540/rawdata/GWAS/diploid/vcf/all

ml BCFtools/1.6-foss-2019b
ml PLINK/1.9b_5-x86_64
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
ml EIGENSOFT/7.2.1-foss-2019b

### merge vcfs based on a list of file names for each strain
bcftools merge -l gain2n-strain.list -Oz -o temp.vcf.gz

### change name of strains, if necessary, based on a list of names called newIDs.txt
bcftools reheader -s newIDs.txt -o gain2n.vcf.gz temp.vcf.gz

### change chromosome names, if necessary, so plink can read it (new chr names in newChrNames.txt)
bcftools annotate --rename-chrs newChrNames.txt -Oz -o gain2n.vcf.gz temp.vcf.gz

### remove mitochondrial chromosome
bcftools view -e 'CHROM=="sacCer3chrM"' -Oz -o gain2n_g.vcf.gz gain2n.vcf.gz
#grep -Ev '^(sacCer3chrM)' gain2n.vcf.gz > gain2n_g.vcf.gz

### copy vcf to folder
cp gain2n_g.vcf.gz /scratch/es47540/rawdata/GWAS/diploid/gainOnly

### change folder to current dir
cd /scratch/es47540/rawdata/GWAS/diploid/gainOnly

### Exclude invariant sites
bcftools view --min-ac 1 -Oz -o gain2n_gSNPs.vcf.gz gain2n_g.vcf.gz

### Create plink format files (ped, map, bim)
vcftools --gzvcf gain2n_gSNPs.vcf.gz --minQ 40 --plink --out gain2n_gSNPs

### Convert ped to bed to save space/time (bed, fam)
plink --file gain2n_gSNPs --make-bed --out gain2n_gSNPs

### QC1. delete SNPs with missingness > 0.1
plink --bfile gain2n_gSNPs --geno 0.1 --make-bed --out gain2n_gSNPs_QC1

### QC2. delete individuals with missingness > 0.1
plink --bfile gain2n_gSNPs_QC1 --mind 0.1 --make-bed --out gain2n_gSNPs_QC2

### QC3. delete SNPs with low MAF (MAF < 0.2%)
plink --bfile gain2n_gSNPs_QC2 --maf 0.002 --make-bed --out gain2n_gSNPs_QC3

### QC4. delete highly-correlated SNPs (recomb hot-spots & R^2 > 0.2 in 50:5 bp)
plink --bfile gain2n_gSNPs_QC3 --exclude hot_spots_total.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile gain2n_gSNPs_QC3 --extract indepSNP.prune.in --make-bed --out gain2n_gSNPs_QC4

### convert .bed to .ped
plink --bfile gain2n_gSNPs_QC4 --recode tab --out gain2n_gSNPs_QC4

### run PCA based on the parameters file (no outlier removal here)
smartpca -p gain2n_gSNPs_QC4.par > PCAgain2n_gSNPs_QC4.log
smartpca -p gain2n_gSNPs_noOut_QC4.par > PCAgain2n_gSNPs_noOut_QC4.log
