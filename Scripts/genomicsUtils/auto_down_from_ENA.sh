awk '{print $2, $1}' study.txt > files.txt
split -d -l 1 files.txt fay
for file in fay*; do perl -ne 'print $_ x 2' $file > ${file}_dup; done
sed -i "1 s/;.*//" *dup
sed -i "2 s/f.*;//" *dup
sed -i "1 s/ ftp/_1.fastq.gz ftp/g" *dup
sed -i "2 s/ ftp/_2.fastq.gz ftp/g" *dup
sed -i 's/^/axel -o /g' *dup
