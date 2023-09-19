# to download sequences from ENA in parallel, see master_download.sh script 
# get a text file with the run acession number in the first column and the strain name in the second column, named strainslist.txt

awk 'print $4,$6' study.txt > strainslist.txt

# batch change filenames from downloaded sequences from ENA 

# ONLY IF THERE ARE MULTIPLE FASTQ FILES FOR THE SAME STRAIN, DO THIS BEFORE:
# awk '$1 in a {$1=$1 "-Part-" ++a[$1]}{a[$1];print}' Leducq_strain_run.txt

# for file in *-2.bam
# do strain_name=${file%%-2.bam}
# mkdir $strain_name
# cp $strain_name*.bam $strain_name/
# done



runlist=`awk '{print $1}' ./strainslist.txt`
strainlist=`awk '{print $2}' ./strainslist.txt`

i=0

for filename in *_1.fastq.gz
do 
strain_name=$(echo $filename | cut -d'_' -f 1)
if (( $strain_name==${runlist[$i]} ))
then mv -v $filename ${strainlist[$i]}_1.fastq.gz
fi
i=$((i+1))
done

i=0

for filename in *_2.fastq.gz
do 
strain_name=$(echo $filename | cut -d'_' -f 1)
if (( $strain_name==${runlist[$i]} ))
then mv -v $filename ${strainlist[$i]}_2.fastq.gz
fi
i=$((i+1))
done

