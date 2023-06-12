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
