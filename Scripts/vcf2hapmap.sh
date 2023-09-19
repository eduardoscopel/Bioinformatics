#!/bin/bash
#SBATCH --job-name=vcf2tassel
#SBATCH --partition=bensasson_p
#SBATCH --ntasks=1
#SBATCH --mem=36G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/es47540/rawdata/GWAS/diploid/output/%x_%j.out
#SBATCH --error=/scratch/es47540/rawdata/GWAS/diploid/error/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eduardoscopel@uga.edu


cd /scratch/es47540/rawdata/GWAS/diploid/vcf/snps


module load TASSEL/5.2.44-Java-1.8.0_144
run_pipeline.pl  -Xms20g -Xmx30g -fork1 -vcf diploids_gSNPs.vcf -export diploids_gSNPs.hapmap -exportType Hapmap
