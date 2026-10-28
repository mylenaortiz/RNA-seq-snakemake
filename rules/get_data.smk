rule get_data:
    output:
        "data/fastq/{sample}_1.fastq",
        "data/fastq/{sample}_2.fastq"
    envmodules:
        "SRA-Toolkit/3.0.10-gompi-2023a"
    threads: 30
    params:
        dir= "data/fastq",
        temp= "data/temp",
        sra= "data/sra"
    shell:
        """
        mkdir -p {params.dir} {params.temp} {params.sra}

        prefetch {wildcards.sample} -O {params.sra} --max-size u

        fasterq-dump {params.sra}/{wildcards.sample} --split-files \
            -O {params.dir} -t {params.temp} -e {threads}

        rm -rf {params.sra}/{wildcards.sample}

        """

rule fastqc:
    input: 
        expand("data/fastq/{{sample}}_{read}.fastq", sample=SAMPLES, read=READS)
    output:
        expand("data/fastq/fastqc/{{sample}}_{read}_fastqc.html",sample=SAMPLES, read=READS)
    envmodules:
        "FastQC/0.12.1-Java-11"
    threads: 40
    params:
        dir= "data/fastq/fastqc"
    shell:
        """
        mkdir -p {params.dir}

        fastqc -t {threads} {input} -o {params.dir}
        """
