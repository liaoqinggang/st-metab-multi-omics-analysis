import re
import os

dir_path = "/pathway/to/01DNA-reseq/"

# find C0xx_1.fq.gz files
file_pattern = r'C0(\d+)_1\.fq\.gz'
filenames = [f for f in os.listdir(dir_path) if re.fullmatch(file_pattern, f)]
samples = [re.sub(r'_1\.fq\.gz$', '', f) for f in filenames]
samples = sorted(set(samples))  # 

rule all:
    input:
        expand(os.path.join(dir_path, "{sample}.mem.mk.sorted.bam"), sample=samples),
        expand(os.path.join(dir_path, "{sample}.mem.metrics"), sample=samples)

rule bwa_mem:
    input:
        fq1=os.path.join(dir_path, "{sample}_1.fq.gz"),
        fq2=os.path.join(dir_path, "{sample}_2.fq.gz"),
        genome="/public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/03_genome/DM8.1_genome.fasta"
    output:
        bam=os.path.join(dir_path, "{sample}.mem.mk.sorted.bam"),     
        metrics=os.path.join(dir_path, "{sample}.mem.metrics")         
    conda:
        "pathway/to/gatk4.yml"
    threads: 8
    shell:
        """
        # 1. bwa mem align + sort
        bwa mem -t {threads} {input.genome} {input.fq1} {input.fq2} | \
        samtools sort -@4 -o {output.bam}

        # 2. Read Group and MarkDuplicates）
        java -jar /pathway/to/00_software/picard.jar AddOrReplaceReadGroups \
            I={output.bam} \
            O=/tmp/{wildcards.sample}.rg.bam \
            RGID={wildcards.sample} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={wildcards.sample}

        
        java -jar /pathway/to/00_software/picard.jar MarkDuplicates \
            I=/tmp/{wildcards.sample}.rg.bam \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000

        # 3. create index
        samtools index {output.bam}

        # 4. clean
        rm -f /tmp/{wildcards.sample}.rg.bam
        """
