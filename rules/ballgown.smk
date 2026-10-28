N = config["ballgown"]["n_plots"]
PLOTS = list(range(1, N+1))

rule ballgown:
    input: 
        counts= expand("data/stringtie_merged/{sample}/e_data.ctab", sample=SAMPLES),
        metadata= config["ballgown"]["metadata"]
    
    output: 
        results_unfiltered= "data/ballgown/isoformtranscript_results.csv",
        results_filtered= "data/ballgown/isoformtranscript_results_filtered.csv",
        bg= "data/ballgown/bg.rda",
        environment = "data/ballgown/environment.RData",
        sig_transcripts= "data/ballgown/signature_genes.csv",
        isoforms_pdf= expand("data/ballgown/{plots}_isoforms.pdf", plots = PLOTS),
        fpkm_boxp= "data/ballgown/fpkm_boxplot.tiff"
    conda:
        "../envs/ballgown.yaml"
    
    log:
        "logs/ballgown.log"

    params:
        outdir= "data/ballgown",
        workdir= "data/stringtie_merged/",
        covariate= config["ballgown"]["grouping_variable"],
        reference_level= config["ballgown"]["reference_level"],
        replicates= config["ballgown"]["replicates"],
        source = config["ballgown"]["source"],
        samples = config["sample_table"],
        n_plots= config["ballgown"]["n_plots"]

    script:
        "../scripts/ballgown.R"