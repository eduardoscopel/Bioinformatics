#!/bin/bash
#SBATCH --partition=bensasson_p
#SBATCH --job-name=raxml-cen1
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --mail-user=eduardoscopel@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/error/%x_%j.err
#SBATCH --output=/scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/output/%x_%j.out

cd /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree9/CEN-tree/cen1
ml RAxML/8.2.11-foss-2019b-pthreads-sse

raxmlHPC-PTHREADS-SSE3 -T 18 \
-f a \
-x 12345 \
-p 12345 \
-m GTRGAMMA \
-o EM14S01_3B,EN14S01,GE14S01_7B \
-N 100 \
-s wg_align.fa \
-n RAxMLtree \
-w /scratch/es47540/rawdata/aneuploidy/peter/peter_tree/tree4/RAxML-tree

raxmlHPC-PTHREADS-SSE3 \
-f a \
-x 12345 \
-p 12345 \
-m GTRGAMMA \
-o EM14S01_3B.CEN1,EN14S01.CEN1,GE14S01_7B.CEN1 \
-N 100 \
-s 621.cen1.fa \
-n cen1

module load IQ-TREE/1.6.5-omp

iqtree -nt 4 \
-s wg_align.fa \
-m GTR+G \
-seed 12345 \
-bo 1 \
-pre boot1
