# FastQC
FastQC is an open source quality control (QC) tool used to assess the quality of high-throughput sequencing data. FastQC needs to be conducted in order to check the quality of provided files and to trim, or to remove any low quality bases which may be found on the provided files.

``` bash
$ srun -c 4 --time=10:00:00 -N 1 --mem=32G -A BLUMBERG_LAB --pty /bin/bash -i
$ module load fastqc
$ cd Z13_MSC_RNASeq / 
$ fastqc -t 6 *.gz 
```

# STAR Aligner
STAR aligner is a software tool used to align high-throughput RNA-sequencing reads to a reference genome. In this case, a genome is to be generated and then each of  the reads are to be aligned to the generated genome.

``` bash
#load the module
$ module load star/2.7.10a  

#request an interactive node
$ srun -c 4 --time=10:00:00 -N 1 --mem=32G -A BLUMBERG_LAB --pty /bin/bash -i

#generating the human genome
$ STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /dfs8/pub/byukimti/reRun/genomeDirect/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf --sjdbOverhang 99

#aligning to the human genome
$ STAR --runMode alignReads --genomeDir GenomeDir/ --outSAMtype BAM Unsorted SortedByCoordinate --readFilesCommand zcat --readFilesIn /dfs6/pub/byukimti/trimmedFiles/nR203-L1-G5-P305-CCGCTGTT-TGCGGTCT-READ1-Sequences.txt.gz_val_1.fq.gz /dfs6/pub/byukimti/trimmedFiles/nR203-L1-G5-P305-CCGCTGTT-TGCGGTCT-READ2-Sequences.txt.gz_val_2.fq.gz --runThreadN 10 --outFileNamePrefix /dfs6/pub/byukimti/alignGenome/finishedAligning/

```


# Rsubread
Rsubread is an R package designed for analysis of RNA-seq data. In this case, the featureCounts function is being used in order to count the number of RNA-seq reads aligned to specific genomic features. 

``` r
#install the package and call the library
install.packages("Rsubread")
library("Rsubread")

#set input and output directory
dir_input <- "/dfs6/pub/byukimti/alignGenome/finishedAligning/sortedBams/"
dir_output <- "/dfs6/pub/byukimti/alignGenome/finishedAligning/rsubread/"

#set working directory 
setwd(dir_output)

#assign the gtf file
gtf.file <- "/dfs6/pub/byukimti/Mus_musculus.GRCm39.108.gtf"

#assign variables to each .bam file
ZO04 <- paste0(dir_input, "304Aligned.sortedByCoord.out.bam")
ZO05 <- paste0(dir_input, "305Aligned.sortedByCoord.out.bam")
ZO06 <- paste0(dir_input, "306Aligned.sortedByCoord.out.bam")
IO07 <- paste0(dir_input, "307Aligned.sortedByCoord.out.bam")
IO08 <- paste0(dir_input, "308Aligned.sortedByCoord.out.bam")
IO09 <- paste0(dir_input, "309Aligned.sortedByCoord.out.bam")
ZI10 <- paste0(dir_input, "310Aligned.sortedByCoord.out.bam")
ZI11 <- paste0(dir_input, "311Aligned.sortedByCoord.out.bam")
ZI12 <- paste0(dir_input, "312Aligned.sortedByCoord.out.bam")
II13 <- paste0(dir_input, "313Aligned.sortedByCoord.out.bam")
II14 <- paste0(dir_input, "314Aligned.sortedByCoord.out.bam")
II15 <- paste0(dir_input, "315Aligned.sortedByCoord.out.bam")
bam.files <- c(ZO04,ZO05,ZO06,IO07,IO08,IO09,ZI10,ZI11,ZI12,II13,II14,II15)

```
# Differential Expression Analysis

Differential expression analysis allows for the identification of genes which are differentially expressed between different conditions of provided sample groups. In this case, differential expression analysis is to be conducted in order to identify differences in upregulated and downregulated genes. 

```r
#load the data from the HPC into R studio
load("/Users/UCI/Downloads/NormalizationData.RData")

#install packages
BiocManager::install("edgeR")
install.packages("rlang")
BiocManager::install("tidyverse")
BiocManager::install("DESeq2")
install.packages("devtools")
BiocManager::install("Rsubread")
install.packages("stats")
install.packages("statmod")
install.packages("tibble")
install.packages("dplyr")
BiocManager::install(“GenomicFeatures”)

#load libraries
library(edgeR)
library(rlang)
library(tidyverse)
library(DESeq2)
library(devtools)
library(Rsubread)
library(stats)
library(statmod)
library(tibble)
library(dplyr)
library(GenomicFeatures)

group_all <- factor(c(rep("ZO", 3), rep("IO", 3), rep("ZI", 3), rep("II", 3)))

y_all <- DGEList(counts=countsRedo$counts,group=group_all)

#filter out genes which have no expression 
keepy_all <- filterByExpr(y_all)
y_all <- y_all[keepy_all,,keep.lib.sizes=FALSE]

#normalize the library size using TMM normalization
y_all <- calcNormFactors(y_all, method = "TMM")

#the design matrix: to determine your control and treatment
design_all <- model.matrix(~ 0 + group_all)

colnames(design_all) <- levels(group_all)

#estimate the dispersion
y_all <- estimateDisp(y_all, design_all, robust=TRUE)

#QL dispersions
fit_all <- glmQLFit(y_all, design_all, robust=TRUE)

#differential expression analysis
# left = treatment group, right = control group
ZOvsIO_filtered <- makeContrasts(ZO - IO, levels=design_all)
ZOvsII_filtered <- makeContrasts(ZO - II, levels=design_all)
ZOvsZI_filtered <- makeContrasts(ZO - ZI, levels=design_all)
ZIvsII_filtered <- makeContrasts(ZI - II, levels=design_all)
ZIvsIO_filtered <- makeContrasts(ZI - IO, levels=design_all)
IOvsII_filtered <- makeContrasts(IO - II, levels=design_all)

#performs quasi-likelihood (QL) ratio tests for differential expression analysis 
res.ZOvsIO <- glmQLFTest(fit_all, contrast=ZOvsIO)
res.ZOvsZI <- glmQLFTest(fit_all, contrast=ZOvsZI)
res.ZOvsII <- glmQLFTest(fit_all, contrast=ZOvsII) 
res.ZIvsII <- glmQLFTest(fit_all, contrast=ZIvsII)
res.ZIvsIO <- glmQLFTest(fit_all, contrast=ZIvsIO)
res.IOvsII <- glmQLFTest(fit_all, contrast=IOvsII)

summary(decideTests(res.ZOvsIO))
summary(decideTests(res.ZOvsZI))
summary(decideTests(res.ZOvsII))
summary(decideTests(res.ZIvsII))
summary(decideTests(res.ZIvsIO))

baseComp <- decideTests(res.ZOvsIO))

baseComp <- as.data.frame(baseComp)
colnames(baseComp) <- "expression"

baseComp$gene <- rownames(baseComp)
baseComp <- baseComp %>% filter(expression != 0) 

```
# Heatmaps
Heatmaps are visual representations of data, and represent patterns, relationships or trends in large data sets as a display matrix of colored cells. They are used to display gene expression patterns using colored cells. In this case, the heatmap is to be created in order to compare the differentially expressed genes between all provided conditions.

```r
#install the reference genome and call the libraries
BiocManager::install("Homo.sapiens")
library("Homo.sapiens")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
require(org.Hs.eg.db)

GeneSymbol <- mapIds(org.Hs.eg.db, keys = rownames(y_all), keytype = "ENSEMBL", column = "ENTREZID")

y_all$gene <- data.frame(ENTREZID = GeneSymbol)
head(y_all$gene)
SYMBOL <- rownames(y_all$gene)
ENTREZID <- y_all$gene[,1] 
y_all$gene <- data.frame(Symbol = SYMBOL, ENTREZID = ENTREZID)

Hs_genes <- transcriptsBy(Homo.sapiens, by="gene", columns=c("SYMBOL", "ENTREZID", "TXCHROM", "TXSTRAND"))

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")
gene.length <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.refGene, by = "gene")

gene.length_2 <- unlist(gene.length)
gene.length_2$ENTREZID <- names(gene.length_2)

names(gene.length_2) <- gene.length_2$tx_name
gene.length <- relist(gene.length_2, gene.length)
gene.length.df <- as.data.frame(gene.length)
gene.length.df
gene.length.df <- gene.length.df[ -c(1:2) ]

gene.length.df.2 = gene.length.df %>% group_by(ENTREZID) %>% top_n(n = 1, wt = width) %>% distinct(ENTREZID, .keep_all = TRUE)
gene.length.df.2
gene.length.df.2$length = gene.length.df.2$width

y_all$gene = y_all$gene %>% left_join(dplyr::select(gene.length.df.2, c("length", "ENTREZID")), by = "ENTREZID")

length(unique(y_all$gene$Symbol))

head(y_all$gene)

install.packages("pheatmap")
library(pheatmap)
install.packages("RColorBrewer")
library(RColorBrewer)

# to plot the heatmap, you need to use edgeR cpm() to scale the data so the difference in reads will be indicated by log(cpm)
DGEList.cpm.2021 <- cpm(y_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

#genes to plot
genes.x <- baseComp$gene

length(genes.x)

df.cpm = as.data.frame(DGEList.cpm.2021)
df.cpm$Symbol <- y_all$gene$Symbol
length(df.cpm$Symbol)

df.cpm.subset <- df.cpm %>% filter(Symbol %in% genes.x)

length(df.cpm.subset$Symbol)
mex.cpm.subset = as.matrix(df.cpm.subset[, 1:12])
rownames(mex.cpm.subset) = df.cpm.subset$Symbol
colnames(mex.cpm.subset) = group_all
colnames(mex.cpm.subset)

install.packages("pheatmap")
pdf("Heatmap_Final.pdf", width = 10, height = 10)
pheatmap(mex.cpm.subset, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         scale = "none", 
         display_numbers = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, 
         main = NA)
```
