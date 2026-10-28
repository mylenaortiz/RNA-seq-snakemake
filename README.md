# RNA-Seq Snakemake Pipeline

RNA-seq analysis using data retrieved from SRA, aligned with HISAT2 (STAR outputs are also produced but not used in the downstream pipeline), and abundances estimated with StringTie. Following, users have the option to perform isoform abundance analysis using Ballgown or DEGs with DESeq2. 
The pipeline has been converted into a Snakemake workflow for the HPCC (MSU). To run the pipeline, make sure that Snakemake has been downloaded and the conda environment is available in run_snakemake.sh. 

The envmodules listed in the .smk rules (available in the 'rules' directory) can be modified to reflect the version available in the HPCC. Modules can also be downloaded through a conda environment (please note that in this case, the rule will need to be modified accordingly).

All inputs and parameters requested by this pipeline can be found in the config.yml. Importantly, set up the stringtie parameter to "expression_estimation" if DEG analysis through DESeq is intended, or to "isoforms" if you desire analysis through Ballgown.

The pipeline in the current configuration is designed to reproduce the figures included in the manuscript by Ortiz et al (2026). However, it can be adapted to RNA-seq analysis using any dataset containing FASTQ files. Inputs should be modified in the config.yml accordingly, as described in this file.

To submit this script to the cluster:

```bash
sbatch run_snakemake.sh
```
