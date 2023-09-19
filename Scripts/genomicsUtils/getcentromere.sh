cat GCF_000146045.2_R64_genomic.gff | awk '$1 {print $0; } ' | grep Note=CEN\.%3B # 1 to 9

cat GCF_000146045.2_R64_genomic.gff | awk '$1 {print $0; } ' | grep Note=CEN\..%3B # >= 10 
