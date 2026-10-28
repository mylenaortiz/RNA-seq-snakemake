rule star_index:
    input:
        fa= config["STAR"]["genome_fasta"],
        gtf= config["STAR"]["genome_gtf"]
    output:
        directory("data/genome/genome_index")
    envmodules:
        "STAR/2.7.11b-GCC-13.2.0"
    threads: 80
    shell:
        """
        mkdir -p {output}

        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99

        echo "genome indexing done!"

        """
    
rule star:
    input:
        r1="data/Trimmed/{sample}.ForwardPaired.fastq.gz",
        r2="data/Trimmed/{sample}.ReversePaired.fastq.gz",
        index="data/genome/genome_index"
    output:
        "data/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    envmodules:
        "STAR/2.7.11b-GCC-13.2.0"
    threads: 80
    params: 
        gtf= config["STAR"]["genome_gtf"],
        outdir= "data/STAR/{sample}"
    shell:
        """
        mkdir -p {params.outdir}

        STAR \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --readFilesCommand zcat \
        --runThreadN {threads} \
        --sjdbGTFfile {params.gtf} \
        --outReadsUnmapped Fastx \
        --outMultimapperOrder Random \
        --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix {params.outdir}/

        echo "STAR alignment done!"

        """

rule samtools:
    input:
        "data/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "data/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    envmodules: 
            "SAMtools/1.23.1-GCC-14.3.0"
    threads: 10
    shell:
        """
        samtools index -@ {threads} {input}

        echo "indexing done!"    
        """


