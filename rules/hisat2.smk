rule hisat:
    input:
        r1="data/Trimmed/{sample}.ForwardPaired.fastq.gz",
        r2="data/Trimmed/{sample}.ReversePaired.fastq.gz"
    output:
        "data/hisat2/{sample}.hisat2_aln.sam"
    envmodules:
        "HISAT2/2.2.1-gompi-2023a"
    threads: 30
    params: 
        genome_index= config["HISAT2"]["genome_index_hisat"],
        outdir= "data/hisat2/"
    shell:
        """
        mkdir -p {params.outdir}

        hisat2 -p {threads} -x {params.genome_index} -1 {input.r1} -2 {input.r2} -S {output}

        echo "HISAT2 alignment done!"

        """
    
rule samtools_2:
    input:
        "data/hisat2/{sample}.hisat2_aln.sam"
    output:
        "data/hisat2/{sample}.hisat2_aln.bam"
    envmodules: 
        "SAMtools/1.23.1-GCC-14.3.0"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        samtools index -@ {threads} {output}

        echo "indexing done!"    
        """
