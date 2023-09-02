# :hammer: Scan
**Scanning sample-specific miRNA regulation from bulk and single-cell RNA-sequencing data**

## :boom: Background
In the gene regulation field, the miRNA regulation has attracted broad attentions on account of its potential clinical translation. Therefore, inferring and characterizing miRNA regulation is a central problem for revealing miRNA regulatory mechanisms in systems biology. Just like dynamic biological system in a biological sample (cell or tissue), the miRNA regulation is also dynamic. For example, by using bulk and single-cell transcriptomics data, existing studies have found that the miRNA regulatory networks tend to be condition-specific, tissue-specific and cell-specific. Moreover, it is widely known that the biological samples are characterized by heterogeneity, and are shown to be unique. Thus, to understand the miRNA regulation specific to biological samples (cancer tissues or cells in particular), it is essential to investigate sample-specific miRNA regulation.

To model the dynamic regulatory processes of miRNAs at the single-sample level, in this work, we propose a **S**ample-spe**c**ific miRNA regul**a**tio**n** (**Scan**) framework to scan sample-specific miRNA regulation in bulk and single-cell RNA-sequencing data. Given bulk or single-cell RNA-sequencing data with or without putative miRNA-mRNA interactions, Scan provides 27 well-used computational methods to calculate the strength of the relationships between miRNAs and mRNAs. In addition, Scan applies two strategies: statistical perturbation and linear interpolation to infer sample-specific miRNA regulatory networks. In the context of multiple network inference methods and in two transcriptomics datasets (bulk data from breast cancer tissues and single-cell data from chronic myelogenous leukemia cells), we demonstrate the accuracy, scalability and practicability of Scan.

A schematic illustration of **Scan** is shown in the folowing.

<p align="center">
  <img src="https://github.com/zhangjunpeng411/Scan/blob/master/Scan_schematic_illustration.png" alt="Schematic illustration of Scan" border="0.1">
</p>

Given bulk or single-cell RNA-sequencing data with or without priori information of miRNA-mRNA interactions, Scan can apply 27 computational methods covering 7 types (correlation, distance, information, regression, bayesian, proportionality and causality) to construct miRNA-mRNA relation matrix. By using one of 27 computational methods, Scan constructs two miRNA-mRNA relation matrices (one for all samples and the other for all samples except *k*). For each sample (cell or tissue) *k*, Scan conducts sample-specific network inference from the two constructed miRNA-mRNA relation matrices to infer a miRNA regulatory network specific to the sample *k*. In total, Scan can identify *m* sample-specific miRNA regulatory networks across *m* samples (one network for one sample). In terms of accuracy and scalability, Scan further evaluates the constructed sample-specific miRNA regulatory networks.

## :book: Description of each file in R and GenMiR++ folders
- **Scan.interp.R**: Utility functions for scanning sample-specific miRNA regulation using a linear interpolation strategy.

- **Scan.perturb.R**: Utility functions for scanning sample-specific miRNA regulation using a statistical perturbation strategy.

- **Case_study.R**: Case study for scanning sample-specific miRNA regulation.

- **evalE.m, GenMiR.m, GenMiR_evalE.m, GenMiR_generic.m, GenMiR_VBEStep.m, GenMiR_VBMStep.m**: Default functions of GenMiR++.

- **GenMiR_SSN.m**: Utility function of scanning sample-specific miRNA regulation using a statistical perturbation strategy, the network inference method is GenMiR++. 

- **matrixzcore.m**: Utility function of calculating zscore of a matrix.

- **Main_GenMiR_Scan.interp.m, Main_GenMiR_Scan.perturb.m, GenMiR.R**: Case study for scanning sample-specific miRNA regulation, the network inference method is GenMiR++.

## :gear: The usage of SCOM
Paste all files into a single folder (set the folder as the directory of R environment). The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("Case_study.R")
```

## :zap: Quick example to use Scan
For scanning sample-specific miRNA regulation, users should prepare matched miRNA and mRNA expression data and putative miRNA-target interactions (optional). Paste the datasets and our source file (**Scan.interp.R** and **Scan.perturb.R**) into a single folder (set the folder as the directory of R environment), users can use the following scripts to scan sample-specific miRNA regulation. For convenience, our bulk and single-cell transcriptomics datasets prepared for users are from [here](https://drive.google.com/file/d/1OUkOJW9TDnGDi0lsntR4by8J5-oZsgU1/view?usp=share_link).

```{r echo=FALSE, results='hide', message=FALSE}
## The dataset "K562_19_single-cell_matched_miR_mR.RData" is used in this quick example
# Load packages
library(pracma)
library(WGCNA)
library(igraph)
library(energy)
library(Hmisc)
library(parmigene)
library(minerva)
library(glmnet)
library(pcalg)
library(doParallel)
library(philentropy)
library(StatMatch)
library(propr)
library(gtools)
library(pbapply)
library(pcaPP)

## Preprocess the single-cell sequencing data including log2(x+1), compute the average expression values of duplicate genes
## and remove genes with constant expression values in all cells

# Transformation using log2(x+1)
miRNA_scRNA_norm <- log2(miRNA_scRNA_raw+1)
mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1)

# Compute the average expression values of duplicate genes
miRNA_scRNA_norm_average <- Averg_Duplicate(miRNA_scRNA_norm)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)
    
# Remove genes with zero expression values in all cells
miRNA_scRNA_norm_mean <- unlist(lapply(seq(dim(miRNA_scRNA_norm_average)[2]), function(i) mean(miRNA_scRNA_norm_average[, i])))
miRNA_scRNA_norm_zero <- miRNA_scRNA_norm_average[, which(miRNA_scRNA_norm_mean > 0)]
mRNA_scRNA_norm_mean <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) mean(mRNA_scRNA_norm_average[, i])))
mRNA_scRNA_norm_zero <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_mean > 0)]
    
# Reserve genes with higher mean expression values in all cells
miRNA_scRNA_norm_mean_update <- unlist(lapply(seq(dim(miRNA_scRNA_norm_zero)[2]), function(i) mean(miRNA_scRNA_norm_zero[, i])))
miRNA_scRNA_norm_filter <- miRNA_scRNA_norm_zero[, which(miRNA_scRNA_norm_mean_update > median(miRNA_scRNA_norm_mean_update))]
mRNA_scRNA_norm_mean_update <- unlist(lapply(seq(dim(mRNA_scRNA_norm_zero)[2]), function(i) mean(mRNA_scRNA_norm_zero[, i])))
mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_zero[, which(mRNA_scRNA_norm_mean_update > median(mRNA_scRNA_norm_mean_update))]

# No priori information
Scan.interp_Pearson_timestart <- Sys.time()
Scan.interp_Pearson_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Pearson")
Scan.interp_Pearson_timeend <- Sys.time()
Scan.interp_Pearson_runningtime_NULL <- Scan.interp_Pearson_timeend - Scan.interp_Pearson_timestart

```    
