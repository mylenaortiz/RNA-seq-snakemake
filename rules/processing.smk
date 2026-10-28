#   expression_estimation -> straight from rule stringtie
#   isoforms              -> re-estimated against the merged transcriptome
#   isoforms only: merge the assemblies of all samples into one transcriptome (runs once)
if config["stringtie"] == "isoforms":
    QUANTDIR = "data/stringtie_merged"
    rule stringtie_merge:
        input:
            gtfs= expand("data/stringtie/{sample}/transcripts.gtf", sample=SAMPLES),
            ref= config["HISAT2"]["genome_gtf_hisat"]
        output:
            merged= QUANTDIR + "/stringtie_merged.gtf",
            gtf_list= QUANTDIR + "/assembly_GTF_list.txt"
        envmodules:
            "StringTie/2.2.1-GCC-11.2.0-Python-2.7.18"
        threads: 20
        params:
            outdir= QUANTDIR,
            mode= config["stringtie"]
        shell:
            """
            mkdir -p {params.outdir}
            
            ls -1 {input.gtfs} > {output.gtf_list}
            cat {output.gtf_list}
            stringtie --merge -p {threads} -G {input.ref} -o {output.merged} {output.gtf_list}
            """

    rule stringtie_requant:
        input:
            bam= "data/hisat2/{sample}.hisat2_aln.bam",
            merged= QUANTDIR + "/stringtie_merged.gtf"
        output:
            QUANTDIR + "/{sample}/transcripts.gtf"
        envmodules:
            "StringTie/2.2.1-GCC-11.2.0-Python-2.7.18"
        threads: 20
        params:
            outdir= QUANTDIR + "/{sample}",
            mode= config["stringtie"]
        shell:
            """
            mkdir -p {params.outdir}
            
            stringtie {input.bam} -e -B -p {threads} -G {input.merged} -o {output} 
            """
elif config["stringtie"] == "expression_estimation":
    QUANTDIR = "data/stringtie"    

rule stringtie_prepde:
    input:
        expand(QUANTDIR + "/{sample}/transcripts.gtf", sample=SAMPLES)
    output:
        gene= QUANTDIR + "/gene_count_matrix.csv",
        transcript= QUANTDIR + "/transcript_count_matrix.csv"
    params:
        quantdir= QUANTDIR
    shell:
        """
        python3 scripts/prepDE.py3 -g {output.gene} -t {output.transcript} -i {params.quantdir}
        """
