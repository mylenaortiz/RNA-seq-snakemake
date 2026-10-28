log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library("DESeq2")
library("sva") 
library("ggplot2")
library("gridExtra")
library("UpSetR")
library("tidyverse")
library("RColorBrewer")
library("gplots")
library("ComplexHeatmap")
library("circlize")
library("fgsea")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("ggrepel")

dir.create(snakemake@params[["outdir"]])

countData_gene <- as.matrix(read.csv(snakemake@input[["counts"]], row.names="gene_id"))
colData <- read.csv(snakemake@input[["metadata"]], row.names = 1)

all (rownames(colData) %in% colnames(countData_gene))
colnames(countData_gene)<- rownames(colData)
countData_gene <- countData_gene[, rownames(colData)]

replicate<- snakemake@params[["replicates"]]
grouping_variable<- snakemake@params[["grouping_variable"]]
source<- snakemake@params[["source"]]

colData$unique_id<- paste0(colData[,source]," ", colData[,grouping_variable], " ", colData[,replicate])
colData$group_2<- paste0(colData[,source]," ", colData[,grouping_variable])

index<- match(colnames(countData_gene),rownames(colData))
colnames(countData_gene)<- colData$unique_id[index]
rownames(colData) <- colData$unique_id

#pca_plot in the raw data

pca_deseq2<- prcomp(countData_gene)  #to calculate values for PCs
pca_deseq2_df<- as.data.frame(pca_deseq2[2]$rotation)
rownames(pca_deseq2_df)<- colnames(countData_gene)
pca_deseq2_df[,"labels"]<-rownames(pca_deseq2_df)

p1 = ggplot(data=pca_deseq2_df, aes(x=PC1, y=PC2, color=labels))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA (uncorrected data)", color="labels")

#Batch correction (for explanation on how to set up your function, view https://rnabio.org/module-03-expression/0003/06/01/Batch-Correction/)

batch <- colData [,snakemake@params[["batch"]]] #assign different numbers for samples processed in different batches, same number if processed in the same
colData$ComBat_group<- as.integer(factor(colData[, grouping_variable])) #unique number for each level of the grouping variable
ComBat_group<- colData$ComBat_group
corrected_data = ComBat_seq(counts = countData_gene, batch = batch, group = ComBat_group)

#pca_plot in the corrected data
corrected_pca_deseq2<- prcomp(corrected_data)  #to calculate values for PCs
corrected_pca_deseq2_df<- as.data.frame(corrected_pca_deseq2[2]$rotation)
corrected_pca_deseq2_df[,"labels"]<-rownames(corrected_pca_deseq2_df)

p2 = ggplot(data=corrected_pca_deseq2_df, aes(x=PC1, y=PC2, color=labels))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA (corrected data)", color="labels")

#Create a DESeqDataSet from count matrix and labels
dds_corrected <- DESeqDataSetFromMatrix (countData = corrected_data,
                                         colData = colData, design = as.formula(snakemake@params[["design"]]))

#Run DESeq2 and generate results table
ref<- snakemake@params[["reference_level"]]
contrast<- snakemake@params[["contrast_level"]]

dds_corrected[[grouping_variable]] <- relevel(
  dds_corrected[[grouping_variable]], ref = ref)
dds_corrected <- DESeq(dds_corrected)
res_corrected <- results(dds_corrected, contrast = c(grouping_variable, contrast, ref))

#Sort by adjusted p-value and display
res_corrected= subset(res_corrected, padj<0.05)
res_corrected <- res_corrected[order(res_corrected$padj), ]

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res_corrected), as.data.frame(counts(dds_corrected,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = snakemake@output[["results_1"]])

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds_corrected)
ddsClean <- DESeq(ddsClean)
resClean <- results(ddsClean, contrast = c(grouping_variable, contrast, ref))
resClean = subset(resClean, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = snakemake@output[["results_2"]])

#SIGNATURE GENES

sign<- as.data.frame(resClean)
sign$symbol<- unlist(rownames(sign), use.names = FALSE)
sign<- sign %>%
  separate(symbol, into = c("ensembl_gene_id", "gene_symbol"), sep = "\\|")
sign<- sign[,-7]
sign<-na.omit(sign)

sign.top <- sign[ (sign$baseMean > 50) & (abs(sign$log2FoldChange) > 2),]
sign.top <- sign.top[order(sign.top$log2FoldChange, decreasing = TRUE),]
write.csv(sign.top, file=snakemake@output[["sig_genes"]], col.names = TRUE)

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(ddsClean, blind=T)
vsd <- varianceStabilizingTransformation(ddsClean, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld)),file = snakemake@output[["rlog"]], sep = '\t')
write.table(as.data.frame(assay(vsd)),file = snakemake@output[["vsd"]], sep = '\t')

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored blue

tiff(snakemake@output[["MA"]], res = 300, width = 2000, height = 2000)
plotMA(ddsClean, ylim=c(-8,8),main = "MA Plot")
dev.off()

#PCA on the normalized counts
vsdata <- vst(ddsClean, blind=FALSE) 

PCA_PLOT<-plotPCA(vsdata, intgroup = "group_2")

png(snakemake@output[["PCA"]], res = 300, width = 2000, height = 2000)
print(PCA_PLOT+ theme_classic())
dev.off()

#heatmap

mat<-assay(rld)[rownames(sign.top), rownames(colData)] #sig genes x samples
colnames(mat) <- colData$unique_id
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
num_keep <- 30
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )
l2_val <- as.matrix(sign.top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"
mean <- as.matrix(sign.top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = sign.top$gene_symbol[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = sign.top$gene_symbol[rows_keep], 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3

tiff(snakemake@output[["heatmap"]], res = 300, width = 5000, height = 6200)
print(h)
dev.off()

#volcano plot

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

DE_genes<-as.data.frame(resClean)
DE_genes$gene<- unlist(rownames(DE_genes), use.names = FALSE)
DE_genes<- DE_genes %>%
  tidyr::separate(gene, into = c("ensembl_gene_id", "gene_symbol"), sep = "\\|")
ids <- sub("\\..*$", "", DE_genes$ensembl_gene_id)

org<-get(snakemake@params[["org_database"]])
DE_genes$entrezgene_id <- mapIds(org,
                                 keys = ids,
                                 keytype = "ENSEMBL",
                                 column = "ENTREZID",
                                 multiVals = "first")

DE_genes_clean <-DE_genes[!is.na(DE_genes$entrezgene_id),]

## Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
DE_genes_clean$diffexpressed <- "NO"
DE_genes_clean$diffexpressed[DE_genes_clean$log2FoldChange > 2 & DE_genes_clean$padj < 0.05] <- "UP"
DE_genes_clean$diffexpressed[DE_genes_clean$log2FoldChange < -2 & DE_genes_clean$padj < 0.05] <- "DOWN"
table(DE_genes_clean$diffexpressed)

DE_genes_clean2 <- DE_genes_clean[order(DE_genes_clean$padj), ]

volcano<- ggplot(DE_genes_clean2, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +  # Corrected yintercept
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "gray", "red3"),
                     labels = c("Down", "Non enriched", "Up")) +
  labs(color = "Expression", 
       x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) +  # Customize the breaks in the x axis
  ggtitle(paste0("Enriched genes in \n", contrast, " vs ", ref)) +  # Plot title 
  coord_cartesian(xlim = c(-20, 20), ylim = c(0, 350))  # Set the x and y axis limits

tiff(snakemake@output[["volcano"]], height = 2000, width = 3200, res=300)
print(volcano)
dev.off()

## DE Pathway Analysis

De_gene.fc <- DE_genes_clean$log2FoldChange
names(De_gene.fc) <- DE_genes_clean$entrezgene_id
gene_list<-na.omit(De_gene.fc)
gene_list = sort(gene_list, decreasing = TRUE)
duplicated<- duplicated(names(gene_list)) #fgsea won't run if there is duplications
gene_list<- gene_list[!duplicated]

gene_sets<- snakemake@input[["fgsea"]]
dataset_names<- snakemake@params[["fgsea_datasets"]]

dotplots <- list()
gsea_list<- list()

for (n in seq_along(gene_sets)){
    set<- gmtPathways(gene_sets[n])
    gsea<- fgsea(pathways = set,
                        stats = gene_list,
                        minSize = 15,
                        maxSize = 500)
    
    gsea_list[[dataset_names[n]]] <- gsea

    df<- as.data.frame(gsea)
    df$leadingEdge<- sapply(df$leadingEdge, paste, collapse = ";") #list column can't be written to csv
    write.csv(df, file = snakemake@output[["fgsea"]][n], row.names = FALSE)

    significant <- subset(gsea, padj < 0.05)  
    significant<- significant[order(significant$NES, decreasing = TRUE), ]  # Order by NES
    significant$pathway <- factor(significant$pathway, levels = significant$pathway)

    dotplot<-ggplot(significant, aes(pathway,NES, size= abs(NES),color= padj, fill= padj))+
        geom_point(stat="identity")+
        coord_flip()+
        labs(title = paste0(dataset_names[n]," Significant Pathways"))+
        theme_minimal()+
        theme(
            axis.title.y = element_text(size = 28),  # Y-axis label
            axis.title.x = element_text(size = 28),  # X-axis label
            axis.text = element_text(size = 28),     # Axis tick labels
            plot.title = element_text(size = 34, face = "bold"),  # Title,
            legend.position = "right",
            legend.text = element_text(size=28),
            legend.title = element_text(size = 28),      # Legend title
            legend.key.size = unit(2, "cm")              # Size of the legend boxes
        )
        
    dotplots[[dataset_names[n]]] <- dotplot

    tiff(snakemake@output[["dotplot"]][n], res = 300, width = 11500, height = 10000)
    print(dotplot)
    dev.off()
}


save.image(snakemake@output[["environment"]])


