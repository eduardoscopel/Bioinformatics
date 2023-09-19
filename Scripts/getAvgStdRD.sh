ml SAMtools/1.6-foss-2019b
for file in *.bam
do
  strain=${file%%.sorted.bam}
  tot=$(samtools view -H $file | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
  #samtools depth -ao cov_$1 $1
  sum=$(awk '{sum+=$3} END {print sum}' ../coverage/$strain.depth)
  avg=$(echo "$sum/$tot" | bc -l)
  avg-std=$(awk '{++x; y+=$3;z+=$3*$3}END{y/=x;z=sqrt(z/x-y*y);print y,z}' ../coverage/$strain.depth)
  echo $strain $avg >> almeidaRD.txt
done

for file in *.depth
do
  avgstd=$(awk '$3 < 1000 {++x; y+=$3;z+=$3*$3}END{y/=x;z=sqrt(z/x-y*y);print y,z}' $file)
  echo $strain $avgstd >> almeidaRDavgStd.txt
done
