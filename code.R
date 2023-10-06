install.packages("BiocManager")

BiocManager::install("ggplot2", "dplyr", "pheatmap", "DESeq2", "AnnotationHub", 
                     "clusterProfiler", "tximport")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(AnnotationHub)
library(clusterProfiler)
library(tximport)


samples<-read.csv("samples.csv", header = T, stringsAsFactors = F, sep="\t")
rownames(samples)<-samples$sample_name

files<-dir("gene", full.names = T, recursive = T)
names(files)<-files
names(files)<-gsub("gene/(14|21)/", "", names(files))
names(files)<-gsub(".genes.results", "", names(files))

names(files) %in% samples$samplename


rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE, countsFromAbundance = "no", abundanceCol = "TPM", 
                 geneIdCol = "gene_id", txIdCol = "transcript_id(s)", countsCol = "expected_count", 
                 lengthCol = "length")

rsem$counts

ah <- AnnotationHub()
query(ah, "EnsDb")
ahDb <- query(ah, pattern = c("Mus Musculus", "EnsDb", 95))
ahEdb<-ahDb[[1]]

annotations<-as.data.frame(genes(ahEdb, columns = c("gene_id", "gene_name", "gene_biotype")))

orgDb <- query(ah, pattern = c("Mus Musculus", "OrgDb"))
orgDb<-orgDb[[1]]


counts<-round(rsem$counts)
tpms<-rsem$abundance

expressed<-rowSums(tpms > 1) >= 1
tpms<-tpms[expressed, ]
counts<-counts[expressed, ]


deseq_genotype<-DESeqDataSetFromMatrix(counts, colData = samples, design = ~ genotype)
deseq_time<-DESeqDataSetFromMatrix(counts, colData = samples, design = ~ time)
deseq_both<-DESeqDataSetFromMatrix(counts, colData = samples, design = ~ genotype+time+genotype*time)

samples$condition<-paste0(samples$genotype, "_", samples$time)
deseq_condition<-DESeqDataSetFromMatrix(counts, colData = samples, design = ~ condition)

de<-DESeq(deseq_condition)

vsd <- vst(de, blind=FALSE, nsub = 10000)
pca<-plotPCA(vsd, intgroup=c("genotype"), ntop=10000, returnData=T)

ggplot(pca, aes(x=PC1, y=PC2, color=genotype, label=name))+geom_text()+theme_minimal()+
  scale_color_brewer(palette = "Set1")

wt_vs_ko_21<-as.data.frame(results(de, contrast = c("condition", "KO_21", "WT_21")))

wt_samples<-samples$samplename[samples$condition=="WT_21"]
ko_samples<-samples$samplename[samples$condition=="KO_21"]

wt_tpm<-rowMeans(log2(tpms[, wt_samples]+1))
ko_tpm<-rowMeans(log2(tpms[, ko_samples]+1))

wt_vs_ko_21$WT<-wt_tpm
wt_vs_ko_21$KO<-ko_tpm

wt_vs_ko_21$gene_id<-rownames(wt_vs_ko_21)
wt_vs_ko_21<-left_join(wt_vs_ko_21, annotations, by="gene_id")
wt_vs_ko_21$sig<-wt_vs_ko_21$padj<0.05

ggplot(wt_vs_ko_21, aes(x=WT, y=KO, color=sig))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+geom_abline(slope = 1, intercept = 0, color="blue")+
  ggtitle("WT vs KO")+xlab("log2(WT_TPM +1)")+ylab("log2(KO_TPM +1)")


ggplot(wt_vs_ko_21, aes(x=log2FoldChange, y=-log10(padj), color=sig))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+ ggtitle("WT vs KO")


sig_genes<-na.omit(wt_vs_ko_21$gene_name[wt_vs_ko_21$sig])
universe<-wt_vs_ko_21$gene_name

go<-enrichGO(gene=sig_genes, universe=universe, OrgDb = orgDb, ont="ALL", keyType = "SYMBOL")
dots<-dotplot(go)
concept<-cnetplot(go, showCategory = 100)


head(as.data.frame(go))

















