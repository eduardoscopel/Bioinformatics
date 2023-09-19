SAMPLES = ["A", "B", "C"]

rule all:
    input:
        "plots/quals.svg"


rule bwa:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
        #fq=expand("data/samples/{sample}_{rep}.fastq.gz", sample=SAMPLES, rep=[1,2])
    output:
        "mapped_reads/{sample}.bam"
    shell:
        """
        ml BWA/0.7.17-GCC-8.3.0
        ml SAMtools/1.6-foss-2019b
        bwa mem {input} | samtools view -Sb - > {output}
        """

rule sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.sorted.bam"
    shell:
        """
        ml SAMtools/1.6-foss-2019b
        samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}
        """

rule samtools_index:
    input:
        "sorted_reads/{sample}.sorted.bam"
    output:
        "sorted_reads/{sample}.sorted.bam.bai"
    shell:
        """
        ml SAMtools/1.6-foss-2019b
        samtools index {input}
        """


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.sorted.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.sorted.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        """
        ml SAMtools/1.6-foss-2019b
        ml BCFtools/1.6-foss-2019b
        samtools mpileup -g -f {input.fa} {input.bam} | bcftools call -mv - > {output}
        """

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    shell:
        """
        ml matplotlib/3.1.1-intel-2019b-Python-3.7.4
        ml Pysam/0.16.0.1-GCC-8.3.0
        ml SciPy-bundle/2020.03-foss-2019b-Python-3.8.2
        ./scripts/plot-quals.py
        """
