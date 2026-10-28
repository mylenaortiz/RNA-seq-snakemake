library("DESeq2")
library("sva") 
library("ggplot2")
library("gridExtra")
library("UpSetR")
library("tidyverse")
library("fgsea")
library(AnnotationDbi)
library(org.Hs.eg.db)

#FILTER GO BP

load("~/DESeq/environment.RData")

c5_gsea<- gsea_list[["c5"]]

significant_c5_bp<- significant_c5[grepl("^GOBP",significant_c5$pathway),]
significant_c5_bp$pathway<- gsub("^GOBP_","", significant_c5_bp$pathway)
significant_c5_bp$pathway <- factor(significant_c5_bp$pathway, levels = significant_c5_bp$pathway)

dot_c5<-ggplot(significant_c5_bp, aes(pathway,NES, size= abs(NES),color= padj, fill= padj))+
  geom_point(stat="identity")+
  coord_flip()+
  labs(title = "GO BP \n Significant Pathways")+
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

dot_c5

tiff("C5 BP_dot_2.tiff", res = 300, width = 7800, height = 7700)
print(dot_c5)
dev.off()

#expression of phosphatase genes (https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md)

##extract PTPs that are DE  
res_table<- resClean

res_table <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table<- res_table %>%
  separate(gene, into = c("ensembl_gene_id", "gene"), sep = "\\|")

res_table<- res_table[,-1]

phosphatases <- res_table[str_detect(res_table$gene, "PTP"),]
phosphatases<- na.omit(phosphatases)
phosphatases_names<-as.vector(phosphatases$gene)
exclude<-c("NAMPTP1","CPTP","PTPRG-AS1")
phosphatases_names<- phosphatases_names[!phosphatases_names %in% exclude]
phosphatases_names<-as.character(phosphatases_names)

normalized_counts<- as.data.frame(counts(dds_corrected),normalized=T)
normalized_counts$symbol<- unlist(rownames(normalized_counts), use.names = FALSE)
normalized_counts<- normalized_counts %>%
  separate(symbol, into = c("ensembl_gene_id", "gene_symbol"), sep = "\\|")
normalized_counts<- na.omit(normalized_counts)
rownames(normalized_counts)<-make.unique(normalized_counts$gene_symbol)
normalized_counts<-normalized_counts[,-16:-17]
rownames(normalized_counts)<- as.character(rownames(normalized_counts))

phosphatase_counts<- normalized_counts[rownames(normalized_counts) %in% phosphatases_names,]
colnames(phosphatase_counts)<-colData$unique_id

gathered_phosphatase <- phosphatase_counts %>%
  gather(colnames(phosphatase_counts), key = "sample", value = "normalized_counts")

phosphatases_names<-as.vector(as.character(rownames(phosphatase_counts)))
gathered_phosphatase$gene<- rep(phosphatases_names,nrow(gathered_phosphatase)/length(phosphatases_names))
gathered_phosphatase<- inner_join(colData, gathered_phosphatase, by= c("unique_id"="sample"))

phosphatase_plot<-ggplot(gathered_phosphatase) +
  geom_point(aes(x = gene, y = normalized_counts, color = genotype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Differentially expressed protein tyrosine phosphatases") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text= element_text(size = 18)
  )

phosphatase_plot

tiff("./phosphatase_plot.tiff", res = 300, width = 5000, height = 2500)
print(phosphatase_plot)
dev.off()

#expression of TKs

TK<- read.csv("tyrosine kinase list.csv")
TK_names<- as.vector(as.character(TK$Approved.symbol))
res_table<- as.data.frame(res_table)
TK_DE <- as.vector(as.character(res_table[res_table$gene %in% TK_names,1]))

TK_counts<- normalized_counts[rownames(normalized_counts) %in% TK_DE,]
colnames(TK_counts)<-colData$unique_id

gathered_TK <- TK_counts %>%
  gather(colnames(TK_counts), key = "sample", value = "normalized_counts")

TK_DE<-as.vector(as.character(rownames(TK_counts)))

gathered_TK$gene<- rep(TK_DE,nrow(gathered_TK)/length(TK_DE))
gathered_TK<- inner_join(colData, gathered_TK, by=c("unique_id"="sample"))

TK_plot<-ggplot(gathered_TK) +
  geom_point(aes(x = gene, y = normalized_counts, color = genotype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Differentially expressed tyrosine kinases") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text= element_text(size = 18)
  )
TK_plot

tiff("./TK_plot.tiff", res = 300, width = 5000, height = 2500)
print(TK_plot)
dev.off()

#gene expression of leading edge genes (cell adhesion)

c5_gsea<-gsea_list[["c5"]] 

leading_edge<- c5_gsea[grepl("GOBP_HOMOPHILIC_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",c5_gsea$pathway),"leadingEdge"]
leading_edge<-unlist(leading_edge, use.names = FALSE)
entrez2gene<-mapIds(org.Hs.eg.db,
                    keys = leading_edge,
                    keytype = "ENTREZID",
                    column = "SYMBOL",
                    multiVals = "first")
adhesion_counts<- normalized_counts[rownames(normalized_counts) %in% entrez2gene,]
colnames(adhesion_counts)<-colData$unique_id

adhesion_counts$gene<- rownames(adhesion_counts)
levels<- res_table %>%
  filter(gene %in% entrez2gene) %>%
  arrange(log2FoldChange) %>%
  dplyr::select(gene)
levels<- levels
levels<-factor(levels[["gene"]])
adhesion_counts$gene<-factor(adhesion_counts$gene, levels=levels)

gathered_adhesion<- pivot_longer(adhesion_counts, -"gene", names_to = "unique_id", values_to = "normalized_counts")
gathered_adhesion<- inner_join(colData, gathered_adhesion, by= "unique_id")

adhesion_plot<-ggplot(gathered_adhesion) +
  geom_point(aes(x = gene, y = normalized_counts, color = genotype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Leading edge: GOBP_HOMOPHILIC_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    title = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text= element_text(size = 16)
  )

adhesion_plot

tiff("./homophilic cell adhesion_plot.tiff", res = 300, width = 5500, height = 1000)
print(adhesion_plot)
dev.off()

#ECM

leading_edge<- c5_gsea[grepl("GOBP_POSITIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",c5_gsea$pathway),"leadingEdge"]
leading_edge<-unlist(leading_edge, use.names = FALSE)
entrez2gene<-mapIds(org.Hs.eg.db,
                    keys = leading_edge,
                    keytype = "ENTREZID",
                    column = "SYMBOL",
                    multiVals = "first")

ECM_counts<- normalized_counts[rownames(normalized_counts) %in% entrez2gene,]
colnames(ECM_counts)<-colData$unique_id
ECM_counts$gene<- rownames(ECM_counts)
levels<- res_table %>%
  filter(gene %in% entrez2gene) %>%
  arrange(log2FoldChange) %>%
  dplyr::select(gene)
levels<- levels
levels<-factor(levels[["gene"]])
ECM_counts$gene<-factor(ECM_counts$gene, levels=levels)

gathered_ECM<- pivot_longer(ECM_counts, -"gene", names_to = "unique_id", values_to = "normalized_counts")
gathered_ECM<- inner_join(colData, gathered_ECM, by="unique_id")

ECM_plot<-ggplot(gathered_ECM) +
  geom_point(aes(x = gene, y = normalized_counts, color = genotype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Leading edge: POSITIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    title = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text= element_text(size = 16)
  )

ECM_plot

tiff("./cell substrate adhesion_plot.tiff", res = 300, width = 5500, height = 1000)
print(ECM_plot)
dev.off()


