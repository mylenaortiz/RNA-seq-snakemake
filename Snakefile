configfile: "config.yml"
READS= config["reads"]

import pandas as pd
SAMPLES= pd.read_csv(config["sample_table"], header=None).iloc[:, 0].tolist()

if config["stringtie"] == "expression_estimation":
    QUANTDIR = "data/stringtie"
elif config["stringtie"] == "isoforms":
    QUANTDIR = "data/stringtie_merged"
else:
    raise ValueError("config['stringtie'] must be 'isoforms' or 'expression_estimation'")

OUTPUTS = []
if config["stringtie"] == "expression_estimation":
    DATASETS = config["DESeq"]["fgsea_datasets"]
    OUTPUTS = [
        "data/DESeq/DeSeq-results-with-normalized.csv",
        "data/DESeq/DeSeq-replaceoutliers-results.csv",
        "data/DESeq/DESeq-rlog-transformed-counts.txt",
        "data/DESeq/DESeq-vst-transformed-counts.txt",
        "data/DESeq/environment.RData",
        "data/DESeq/DESeq2-MAplot_initial_analysis.tiff",
        "data/DESeq/PCA_normalizedcounts.png",
        "data/DESeq/heatmap_top30.png",
        "data/DESeq/volcano_DEG.tiff",
        "data/DESeq/signature_genes.csv",
    ] + expand("data/DESeq/{dataset}_fgsea.csv", dataset=DATASETS) \
      + expand("data/DESeq/{dataset}_dotplot.tiff", dataset=DATASETS)
elif config["stringtie"] == "isoforms":
    N = config["ballgown"]["n_plots"]
    PLOTS = list(range(1, N+1))
    OUTPUTS = [
        "data/ballgown/isoformtranscript_results.csv",
        "data/ballgown/isoformtranscript_results_filtered.csv",
        "data/ballgown/bg.rda",
        "data/ballgown/environment.RData",
        "data/ballgown/signature_genes.csv",
        expand("data/ballgown/{plots}_isoforms.pdf", plots = PLOTS),
        "data/ballgown/fpkm_boxplot.tiff"
    ]

rule all:
    input: 
        expand("data/fastq/fastqc/{sample}_{read}_fastqc.html",
               sample=SAMPLES, read=READS),
        expand(["data/Trimmed/fastqc/{sample}.ForwardPaired_fastqc.html",
                "data/Trimmed/fastqc/{sample}.ForwardUnpaired_fastqc.html",
                "data/Trimmed/fastqc/{sample}.ReversePaired_fastqc.html",
                "data/Trimmed/fastqc/{sample}.ReverseUnpaired_fastqc.html"],sample=SAMPLES),
        "data/genome/genome_index",
        expand("data/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        QUANTDIR + "/gene_count_matrix.csv",
        OUTPUTS

include: "rules/get_data.smk"
include: "rules/trimmomatic.smk"
include: "rules/star.smk"
include: "rules/hisat2.smk"
include: "rules/stringtie.smk"
include: "rules/processing.smk"

if config["stringtie"] == "expression_estimation":
    include: "rules/DEseq.smk"
elif config["stringtie"] == "isoforms":
    include: "rules/ballgown.smk"
