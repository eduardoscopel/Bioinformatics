### CLUSTALW2

clustalw2 -INFILE=wg_align.fa \
-TREE \
-SEED=1 \
-KIMURA \
-BOOTLABELS=node \
-CLUSTERING=NJ \
-BOOTSTRAP=10 \
-OUTPUTTREE=nexus

### MEGA
# Create .mao file using MEGA GUI > prototype > phylogeny > construct NJ tree > change settings > save

megacc -a infer_NJ_nucleotide.mao -d ../wg_align.fa -o test
