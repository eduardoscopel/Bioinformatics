# extract ftp link to download a batch of strains from ENA, using a text file with the format 
# study_accession	sample_accession	secondary_sample_accession	sample_alias	experiment_accession	run_accession	scientific_name	instrument_model	lib# rary_layout	library_source	library_selection	read_count	base_count	experiment_title	fastq_ftp

awk '{for(i=1;i<=NF;i++){if($i~/^ftp/){a=$i}} print $4,$6,a}' study.txt > study_summary.txt
sed -i 's/;/ /g' study_summary.txt
awk '{print $3, $4}' study_summary.txt > study_files.txt
sed -i 's/ftp/"ftp/g' study_files.txt
sed -i 's/gz/gz"/g' study_files.txt
sed -i $'s/ /\\\n/g' study_files.txt
sed -i 's/^/wget -t 0/g' study_files.txt
tail -n +3 study_files.txt > do_get_files

# split do_get_files in m files with n lines, where m=(number of lines in do_get_files)/n, where each split file will have a specific prefix

split -d -l n do_get_files prefix

