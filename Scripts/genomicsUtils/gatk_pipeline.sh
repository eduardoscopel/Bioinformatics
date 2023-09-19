#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N gatk_hh0
#PBS -l walltime=1:00:00,mem=1gb,nodes=1:ppn=1
#PBS -M eduardoscopel@uga.edu
#PBS -m ae
#PBS -o /scratch/es47540/rawdata/contamination/scer/output
#PBS -e /scratch/es47540/rawdata/contamination/scer/error

cd /scratch/es47540/rawdata/contamination/scer/gatk/hap_hap
mkdir error
mkdir output
mkdir scripts

for file in *.sorted.bam
do

  mix=${file%%.sorted.bam}

  cd /scratch/es47540/rawdata/contamination/scer/gatk/hap_hap/scripts

  echo "#PBS -S /bin/bash" >> ${mix}_gatk.sh
  echo "#PBS -q bensasson_q" >> ${mix}_gatk.sh
  echo "#PBS -N ${mix}_gatk" >> ${mix}_gatk.sh
  echo "#PBS -l walltime=12:00:00,mem=32gb,nodes=1:ppn=4" >> ${mix}_gatk.sh
  echo "#PBS -M eduardoscopel@uga.edu" >> ${mix}_gatk.sh
  echo "#PBS -m ae" >> ${mix}_gatk.sh
  echo "#PBS -o /scratch/es47540/rawdata/contamination/scer/gatk/hap_hap/output" >> ${mix}_gatk.sh
  echo "#PBS -e /scratch/es47540/rawdata/contamination/scer/gatk/hap_hap/error" >> ${mix}_gatk.sh
  echo "cd /scratch/es47540/rawdata/contamination/scer/gatk/hap_hap" >> ${mix}_gatk.sh
  echo "module load picard/2.16.0-Java-1.8.0_144" >> ${mix}_gatk.sh
  echo "module load SAMtools/1.6-foss-2016b" >> ${mix}_gatk.sh
  echo "module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8" >> ${mix}_gatk.sh
  echo "java -jar picard.jar ValidateSamFile \
          I=${mix}.sorted.bam \
          IGNORE_WARNINGS=true \
          MODE=VERBOSE" >> ${mix}_gatk.sh
  echo "java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups \
          I=$file \
          O=${mix}.corrected.bam \
          RGLB=lib1 \
          RGPL=illumina \
          RGPU=unit1 \
          RGSM=$mix" >> ${mix}_gatk.sh
  echo "samtools index ${mix}.corrected.bam" >> ${mix}_gatk.sh
  echo "gatk --java-options "-Xmx32g" HaplotypeCaller \
          -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta -ERC GVCF \
          -I ${mix}.corrected.bam \
          -O ${mix}.annotated.g.vcf.gz \
          -A AlleleFraction \
          -A Coverage \
          -A DepthPerAlleleBySample \
          -A DepthPerSampleHC \
          -A GenotypeSummaries \
          -A MappingQuality \
          -A MappingQualityZero \
          -A RMSMappingQuality \
          -A QualByDepth \
          -A UniqueAltReadCount \
          --max-reads-per-alignment-start 0 \
          -mbq 13 \
          --output-mode EMIT_ALL_SITES \
          --standard-min-confidence-threshold-for-calling 0 \
          --max-alternate-alleles 6" >> ${mix}_gatk.sh
  echo "gatk --java-options "-Xmx32g" SelectVariants \
            -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta \
            -V ${mix}.annotated.g.vcf.gz \
            --select-type-to-exclude INDEL \
            -O ${mix}.snps.vcf" >> ${mix}_gatk.sh
  echo "gatk VariantFiltration \
            -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta \
            --variant ${mix}.snps.vcf \
            --filter-name "snpsfilter" \
            --filter-expression "QD<2.0 || MQ<40.0 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0" \
            --output ${mix}.snps.tagged.vcf &> ${mix}.varfilter.out &" >> ${mix}_gatk.sh

  qsub ${mix}_gatk.sh

  done


###  gatk --java-options "-Xmx32g" HaplotypeCaller \
###          -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta -ERC BP_RESOLUTION \
###          -I input.bam \
###          -O output.annotated.vcf.gz \
###          --max-reads-per-alignment-start 0 \
###          -mbq 13 \
###          --output-mode EMIT_ALL_SITES \
###          -ploidy 1 \
###          --standard-min-confidence-threshold-for-calling 0 \
###          --enable-all-annotations TRUE \
###          --all-site-pls TRUE \
###          --max-alternate-alleles 6


###  gatk --java-options "-Xmx4g" GenotypeGVCFs \
###          -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta \
###          -V output.g.vcf.gz \
###          -O genotyped.vcf.gz \
###          -A AlleleFraction \
###          -A Coverage \
###          -A DepthPerAlleleBySample \
###          -A DepthPerSampleHC \
###          -A GenotypeSummaries \
###          -A MappingQuality \
###          -A MappingQualityZero \
###          -A RMSMappingQuality \
###          -A QualByDepth \
###          -A UniqueAltReadCount \
###          -all-sites TRUE

###  gatk FastaAlternateReferenceMaker \
###          -R /scratch/es47540/rawdata/sacCer3/sacCer3.fasta \
###          -O consensus.fa \
###          -V snps.vcf \
