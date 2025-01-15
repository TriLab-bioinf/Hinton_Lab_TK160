# vim: set ft=python:

# RIPseq Single-end workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a bowtie2 index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os
import glob

configfile: "config/config.yaml"
samples = config["samples"].keys()
genome = config["reference"]["genome_file"]
annotation = config["reference"]["ensembl_gtf"]
adapters = config["reference"]["adapters"]

if config["remove_duplicated_reads"]:
    rm_dup_flag =  "--ignoreDup"
else:
    rm_dup_flag = ""

# Functions
def get_fq1(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00-reads/{wildcards.sample}.fastq.gz")]

# Set what rules to run locally
localrules: all #,
            #build_abundant_db

rule all:
    # IMPORTANT: output file for all rules has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input: "results/06-multiqc/multiqc_report.html"  

#rule trimming:
#    input:  fq1 = get_fq1,
#            adapt = f"{adapters}"
#    output: fqU = "results/01-trim/{sample}.trimmed.fastq.gz"
#    resources:
#        cpus_per_task = 8,
#        partition = "quick",
#        time = "4:00:00"
#    threads: 8 
#    log:    log1 = "results/01-trim/{sample}.log",
#            report = "results/01-trim/{sample}.html",
#            reportjson = "results/01-trim/{sample}.json"
#    benchmark:
#            "benchmarks/trim/{sample}.tsv"
#    shell:
#        """
#        module load fastp 
#
#        fastp -i {input.fq1} \
#            -o {output.fqU} \
#            --html {log.report} \
#            --json {log.reportjson} \
#            --qualified_quality_phred 20 \
#            --length_required 25 \
#            --adapter_fasta {input.adapt} \
#            --cut_right --cut_mean_quality 20 --cut_window_size 5 \
#            --thread {threads} > {log.log1} 2>&1
#        """ # ADDED --trim_front1 2 due to low quality values for first two bases of reads  

rule trimming:
    input:  fq1 = get_fq1
    output:
            fq1P = "results/01-trim/{sample}.trimmed.fastq.gz" 
    params: jar="/usr/local/apps/trimmomatic/0.39/trimmomatic-0.39.jar",
            read ="SE",
            adapter=f"ILLUMINACLIP:{adapters}:2:30:12",
            leading="LEADING:0",
            trailing="TRAILING:0",
            slidingwindow="SLIDINGWINDOW:4:20",
            minlen="MINLEN:20",
            headcrop="HEADCROP:0"
    resources:
        cpus_per_task = 16,
        partition = "quick",
        time = "4:00:00"
    threads: 16
    log:    log1 = "results/01-trim/{sample}.log"
    benchmark:
            "benchmarks/trim/{sample}.tsv"
    shell:
        """
        trimmomatic {params.read} \
                -threads {threads} \
                {input.fq1} \
                {output.fq1P} \
                {params.adapter} \
                {params.leading} \
                {params.headcrop} \
                {params.trailing} \
                {params.slidingwindow} \
                {params.minlen} 2>{log.log1}
        """


rule make_bowtie_index:
    input: gen = f"{genome}"
    output: db = "data/00-ref/bowtie_db",
    resources: 
        mem_mb = 32000,
        cpus_per_task = 16,
        partition = "quick",
        time = "3:00:00",
        gres = "lscratch:20"
    threads: 8
    params: 
    shell:
        """
            bowtie2-build --threads {threads} {input.gen} {output.db} 
            
            touch {output.db}
        """

rule map_reads:
    input: fq1 = "results/01-trim/{sample}.trimmed.fastq.gz",
           db = "data/00-ref/bowtie_db"
    output: bam = "results/03-map_reads/{sample}.bam",
            bai = "results/03-map_reads/{sample}.bam.bai"
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 64000,
        gres = "lscratch:20"
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    params: genome_dir = "data/00-ref",
            prefix = "results/03-map_reads/{sample}",
            db_base = "data/00-ref/bowtie_db"
    log: "results/03-map_reads/{sample}.log"        
    shell:
        """
        bowtie2 -x {input.db} \
        --threads {threads} \
        --phred33 \
        --rg-id {params.prefix} --rg SM:{params.prefix} --rg LB:library --rg PL:ILLUMINA \
        -U {input.fq1} \
        --sensitive | samtools sort -@ 8 -O BAM -o {output.bam} -
        
        samtools index -@ 8 {output.bam}
        """

rule filter_reads:
    input: bam = "results/03-map_reads/{sample}.bam"
    output: filtered_bam = "results/04-filtered_bam/{sample}.filter.bam"
    log: "results/04-filtered_bam/{sample}.filter.log"
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "8:00:00"
    shell:
        """
        sambamba view -t {threads} -h -f bam -F "[XS] == null and not unmapped  and not duplicate" \
                {input.bam} > {output.filtered_bam} 2>{log}
        """

rule make_bigwig:
    input: "results/04-filtered_bam/{sample}.filter.bam" 
    output: "results/05-bigwig/{sample}.bw"
    params: "--binSize 10 --normalizeUsing BPM" #  + "--filterRNAstrand [forward/reverse]" to plot strand-specific data
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "14:00:00"
    shell:
        """
        samtools index {input}

        bamCoverage -p {threads} -b {input} -o {output} {params}
        """

rule fastqc:
    input: raw = get_fq1,
           trim = "results/01-trim/{sample}.trimmed.fastq.gz"
    output: 
            o1 = "results/06-fastqc_raw/{sample}_fastqc.html",
            o2 = "results/06-fastqc_trim/{sample}.trimmed_fastqc.html",
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000
    params:
        "--quiet"
    shell:
        """
        fastqc {params} -t {threads} -o results/06-fastqc_raw {input.raw}
        fastqc {params} -t {threads} -o results/06-fastqc_trim {input.trim}
        """
#results/03-map_reads/RNAseLow_plMotB1.2.fastq.gz.Log.final.out
rule multiqc:
    input: 
           i1 = expand("results/01-trim/{s}.log", s = samples),
           i3 = expand("results/03-map_reads/{s}.log", s = samples),
           i5 = expand("results/04-filtered_bam/{s}.filter.log", s = samples),
           i7 = expand("results/05-bigwig/{s}.bw", s = samples),
           i8 = expand("results/06-fastqc_raw/{s}_fastqc.html", s = samples),
           i9 = expand("results/06-fastqc_trim/{s}.trimmed_fastqc.html", s = samples)
    output: "results/06-multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000                    
    shell:
        """
        multiqc -f -d -o results/06-multiqc results/ 
        """

