rule trimmomatic:
    input:
        expand("data/fastq/{{sample}}_{read}.fastq", sample=SAMPLES, read=READS)
    output:
        "data/Trimmed/{sample}.ForwardPaired.fastq.gz",
        "data/Trimmed/{sample}.ForwardUnpaired.fastq.gz",
        "data/Trimmed/{sample}.ReversePaired.fastq.gz",
        "data/Trimmed/{sample}.ReverseUnpaired.fastq.gz"
    envmodules:
        "Trimmomatic/0.39-Java-11"
    threads: 80
    params:
        outdir= "data/Trimmed",
        illuminaclip= config["Trimmomatic"]["illuminaclip"]
    shell:
        """
        mkdir -p {params.outdir}

        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads {threads} -phred33 \
        {input} \
        {output} \
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/{params.illuminaclip} LEADING:15 TRAILING:15 MINLEN:36
        """

rule fastqc_2:
    input: 
        "data/Trimmed/{sample}.ForwardPaired.fastq.gz",
        "data/Trimmed/{sample}.ForwardUnpaired.fastq.gz",
        "data/Trimmed/{sample}.ReversePaired.fastq.gz",
        "data/Trimmed/{sample}.ReverseUnpaired.fastq.gz"
    output:
        "data/Trimmed/fastqc/{sample}.ForwardPaired_fastqc.html",
        "data/Trimmed/fastqc/{sample}.ForwardUnpaired_fastqc.html",
        "data/Trimmed/fastqc/{sample}.ReversePaired_fastqc.html",
        "data/Trimmed/fastqc/{sample}.ReverseUnpaired_fastqc.html"
    envmodules:
        "FastQC/0.12.1-Java-11"
    threads: 40
    params:
        dir= "data/Trimmed/fastqc"
    shell:
        """
        mkdir -p {params.dir}

        fastqc -t {threads} {input} -o {params.dir}
        """

    