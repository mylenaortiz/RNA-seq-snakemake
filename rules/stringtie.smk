rule stringtie:
    input:
        bam="data/hisat2/{sample}.hisat2_aln.bam",
        gtf= config["HISAT2"]["genome_gtf_hisat"]
    output:
        "data/stringtie/{sample}/transcripts.gtf"
    envmodules:
        "StringTie/2.2.1-GCC-11.2.0-Python-2.7.18"
    threads: 20
    params: 
        outdir= "data/stringtie/{sample}",
        mode= config["stringtie"]
    shell:
        """
        mkdir -p {params.outdir}

        if [ "{params.mode}" == "expression_estimation" ]; then
            stringtie {input.bam} -e -B -p {threads} -G {input.gtf} -o {output} -A {params.outdir}/gene_abundances.tsv
        elif [ "{params.mode}" == "isoforms" ]; then
            stringtie {input.bam} -p {threads} -G {input.gtf} -l {wildcards.sample} -o {output}
        else
            echo "Unknown stringtie mode: {params.mode}" >&2
            exit 1
        fi

        """

