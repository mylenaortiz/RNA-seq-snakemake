#!/bin/bash --login
#SBATCH --time=90:00:00                 # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                      # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=80               # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=300G                       # memory required per node - amount of memory (in bytes)
#SBATCH --job-name rnaseq_smk    #can give your job a name for easier identification (same as -J)
#SBATCH --partition=general-long
#SBATCH --mail-user=olivei28@msu.edu
#SBATCH --mail-type=FAIL,BEGIN,END

module purge
module load Miniforge3
conda activate snakemake

cd /mnt/gs21/scratch/olivei28/rnaseq/

snakemake all -p --cores 80 --use-envmodules --use-conda --resources mem_gb=300 --rerun-incomplete
