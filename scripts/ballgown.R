#Code from https://rnabio.org/module-03-expression/0003/03/01/Differential_Expression/
#Adapted by MO on 03/14/2024
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(R6)
library(tidyverse)
library(desc)
library(usethis)
library(ggplot2)

#loading the stringtie output
dir.create(snakemake@params[["outdir"]])

dir<-snakemake@params[["workdir"]]

covariate<- snakemake@params[["covariate"]]
ref<- snakemake@params[["reference_level"]]

phenoData<- read.csv(snakemake@input[["metadata"]]) #you create the phenodata with information about your samples
phenoData[,covariate] <- relevel(factor(phenoData[,covariate]), ref = ref)

replicate<- snakemake@params[["replicates"]]
source<- snakemake@params[["source"]]

phenoData$unique_id<- paste0(phenoData[,source]," ", phenoData[,covariate], " ", phenoData[,replicate])

samples<- read_table(snakemake@params[["samples"]], col_names = FALSE)[[1]] #one SRR ID per line, no header
samples<- file.path(dir, samples) #ballgown needs full paths to each sample folder

bg = ballgown(samples= samples, pData=phenoData, meas='all')

## DE
# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[,9:10])
bg_transcript_names = unique(bg_table[,c(1,6)])

# Save the ballgown object to a file for later use
save(bg, file=snakemake@output[["bg"]])

# Perform differential expression (DE) analysis with no filtering, at both gene and transcript level
results_transcripts = stattest(bg, feature="transcript", covariate= covariate, meas="FPKM", getFC = TRUE)
results_transcripts = merge(results_transcripts, bg_transcript_names, by.x="id", by.y="t_id")

#Subset q< 0.1

DEres<- results_transcripts %>%
  filter(qval < 0.1)

# Save a tab delimited file for both the transcript and gene results
write.csv(results_transcripts, snakemake@output[["results_unfiltered"]])

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
bg_filt_transcript_names = unique(bg_filt_table[,c(1,6)])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate=covariate, meas="FPKM", getFC = TRUE)
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
results_transcripts = merge(results_transcripts, bg_filt_transcript_names, by.x="id", by.y="t_id")

# Output the filtered list of genes and transcripts and save to tab delimited files
write.csv(results_transcripts, snakemake@output[["results_filtered"]])

# Identify the significant genes with qval < 0.1
sig_transcripts = subset(results_transcripts, results_transcripts$qval<0.1)
sig_transcripts<- sig_transcripts %>%
  arrange(qval)

# Output the significant gene results to a pair of tab delimited files
write.csv(sig_transcripts, snakemake@output[["sig_transcripts"]])

##visualization

n_plots<- snakemake@params[["n_plots"]]
genes<- sig_transcripts[1:n_plots,"geneIDs"]
plots_names<- snakemake@output[["isoforms_pdf"]]

for (i in seq_along(genes)){
  pdf(plots_names[i], width = 8, height = 5)
  plotMeans(genes[i], bg_filt, groupvar = covariate, meas = "FPKM", colorby = "transcript")
  dev.off()
}

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg_filt,meas="FPKM")

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

colnames(fpkm)<-phenoData$unique_id

# Create boxplots to display summary statistics for the FPKM values for each sample
tiff(snakemake@output[["fpkm_boxp"]], res=300, height=2000, width=2000)
boxplot(fpkm,col=as.numeric(phenoData[,covariate])+1,las=2,ylab='log2(FPKM+1)')
dev.off()

save.image(snakemake@output[["environment"]])
