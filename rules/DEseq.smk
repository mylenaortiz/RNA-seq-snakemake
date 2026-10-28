DATASETS = config["DESeq"]["fgsea_datasets"]

rule DESeq:
    input: 
        counts="data/stringtie/gene_count_matrix.csv",
        metadata= config["DESeq"]["metadata"],
        fgsea= expand(config["DESeq"]["fgsea"], dataset=DATASETS)
    
    output: 
        results_1= "data/DESeq/DeSeq-results-with-normalized.csv",
        results_2= "data/DESeq/DeSeq-replaceoutliers-results.csv",
        rlog= "data/DESeq/DESeq-rlog-transformed-counts.txt",
        vsd = "data/DESeq/DESeq-vst-transformed-counts.txt",
        environment = "data/DESeq/environment.RData",
        PCA = "data/DESeq/PCA_normalizedcounts.png",
        MA = "data/DESeq/DESeq2-MAplot_initial_analysis.tiff",
        heatmap = "data/DESeq/heatmap_top30.png", 
        volcano = "data/DESeq/volcano_DEG.tiff",
        sig_genes= "data/DESeq/signature_genes.csv",
        fgsea = expand("data/DESeq/{dataset}_fgsea.csv", dataset=DATASETS),
        dotplot = expand("data/DESeq/{dataset}_dotplot.tiff", dataset=DATASETS)
    
    conda:
        "../envs/DEseq.yaml"
    
    log:
        "logs/DESEQ.log"

    params:
        outdir= "data/DESeq",
        org_database= config["DESeq"]["org_database"],
        fgsea_datasets= config["DESeq"]["fgsea_datasets"],
        reference_level= config["DESeq"]["reference_level"],
        contrast_level= config["DESeq"]["contrast_level"],
        grouping_variable= config["DESeq"]["grouping_variable"],
        replicates= config["DESeq"]["replicates"],
        source = config["DESeq"]["source"],
        batch= config["DESeq"]["batch"],
        design= config["DESeq"]["design_formula"]
    
    script:
        "../scripts/DESeq.R"