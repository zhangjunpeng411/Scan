######################################################################################################################################################################################
########################################################################### Scan.interp application in K562 dataset ##################################################################
######################################################################################################################################################################################

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


set.seed(123)
ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

# No prior information
Scan.interp_Pearson_timestart <- Sys.time()
Scan.interp_Pearson_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Pearson")
Scan.interp_Pearson_timeend <- Sys.time()
Scan.interp_Pearson_runningtime_NULL <- Scan.interp_Pearson_timeend - Scan.interp_Pearson_timestart
Scan.interp_Pearson_runningtime_NULL_numeric <- as.numeric(Scan.interp_Pearson_runningtime_NULL)

Scan.interp_Spearman_timestart <- Sys.time()
Scan.interp_Spearman_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Spearman")
Scan.interp_Spearman_timeend <- Sys.time()
Scan.interp_Spearman_runningtime_NULL <- Scan.interp_Spearman_timeend - Scan.interp_Spearman_timestart

Scan.interp_Kendall_timestart <- Sys.time()
Scan.interp_Kendall_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Kendall")
Scan.interp_Kendall_timeend <- Sys.time()
Scan.interp_Kendall_runningtime_NULL <- Scan.interp_Kendall_timeend - Scan.interp_Kendall_timestart

Scan.interp_Dcor_timestart <- Sys.time()
Scan.interp_Dcor_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Dcor")
Scan.interp_Dcor_timeend <- Sys.time()
Scan.interp_Dcor_runningtime_NULL <- Scan.interp_Dcor_timeend - Scan.interp_Dcor_timestart

Scan.interp_RDC_timestart <- Sys.time()
Scan.interp_RDC_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "RDC")
Scan.interp_RDC_timeend <- Sys.time()
Scan.interp_RDC_runningtime_NULL <- Scan.interp_RDC_timeend - Scan.interp_RDC_timestart

Scan.interp_Hoeffd_timestart <- Sys.time()
Scan.interp_Hoeffd_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Hoeffd")
Scan.interp_Hoeffd_timeend <- Sys.time()
Scan.interp_Hoeffd_runningtime_NULL <- Scan.interp_Hoeffd_timeend - Scan.interp_Hoeffd_timestart

Scan.interp_Zscore_timestart <- Sys.time()
Scan.interp_Zscore_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Zscore")
Scan.interp_Zscore_timeend <- Sys.time()
Scan.interp_Zscore_runningtime_NULL <- Scan.interp_Zscore_timeend - Scan.interp_Zscore_timestart

Scan.interp_Biweight_timestart <- Sys.time()
Scan.interp_Biweight_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Biweight")
Scan.interp_Biweight_timeend <- Sys.time()
Scan.interp_Biweight_runningtime_NULL <- Scan.interp_Biweight_timeend - Scan.interp_Biweight_timestart

Scan.interp_Weighted_rank_timestart <- Sys.time()
Scan.interp_Weighted_rank_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Weighted_rank")
Scan.interp_Weighted_rank_timeend <- Sys.time()
Scan.interp_Weighted_rank_runningtime_NULL <- Scan.interp_Weighted_rank_timeend - Scan.interp_Weighted_rank_timestart

Scan.interp_Cosine_timestart <- Sys.time()
Scan.interp_Cosine_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Cosine")
Scan.interp_Cosine_timeend <- Sys.time()
Scan.interp_Cosine_runningtime_NULL <- Scan.interp_Cosine_timeend - Scan.interp_Cosine_timestart

Scan.interp_Euclidean_timestart <- Sys.time()
Scan.interp_Euclidean_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Euclidean")
Scan.interp_Euclidean_timeend <- Sys.time()
Scan.interp_Euclidean_runningtime_NULL <- Scan.interp_Euclidean_timeend - Scan.interp_Euclidean_timestart

Scan.interp_Manhattan_timestart <- Sys.time()
Scan.interp_Manhattan_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Manhattan")
Scan.interp_Manhattan_timeend <- Sys.time()
Scan.interp_Manhattan_runningtime_NULL <- Scan.interp_Manhattan_timeend - Scan.interp_Manhattan_timestart

Scan.interp_Canberra_timestart <- Sys.time()
Scan.interp_Canberra_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Canberra")
Scan.interp_Canberra_timeend <- Sys.time()
Scan.interp_Canberra_runningtime_NULL <- Scan.interp_Canberra_timeend - Scan.interp_Canberra_timestart

Scan.interp_Chebyshev_timestart <- Sys.time()
Scan.interp_Chebyshev_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Chebyshev")
Scan.interp_Chebyshev_timeend <- Sys.time()
Scan.interp_Chebyshev_runningtime_NULL <- Scan.interp_Chebyshev_timeend - Scan.interp_Chebyshev_timestart

Scan.interp_Dice_timestart <- Sys.time()
Scan.interp_Dice_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Dice")
Scan.interp_Dice_timeend <- Sys.time()
Scan.interp_Dice_runningtime_NULL <- Scan.interp_Dice_timeend - Scan.interp_Dice_timestart

Scan.interp_Jaccard_timestart <- Sys.time()
Scan.interp_Jaccard_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Jaccard")
Scan.interp_Jaccard_timeend <- Sys.time()
Scan.interp_Jaccard_runningtime_NULL <- Scan.interp_Jaccard_timeend - Scan.interp_Jaccard_timestart

Scan.interp_Mahalanobis_timestart <- Sys.time()
Scan.interp_Mahalanobis_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Mahalanobis")
Scan.interp_Mahalanobis_timeend <- Sys.time()
Scan.interp_Mahalanobis_runningtime_NULL <- Scan.interp_Mahalanobis_timeend - Scan.interp_Mahalanobis_timestart

Scan.interp_MI_timestart <- Sys.time()
Scan.interp_MI_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "MI")
Scan.interp_MI_timeend <- Sys.time()
Scan.interp_MI_runningtime_NULL <- Scan.interp_MI_timeend - Scan.interp_MI_timestart

Scan.interp_MIC_timestart <- Sys.time()
Scan.interp_MIC_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "MIC")
Scan.interp_MIC_timeend <- Sys.time()
Scan.interp_MIC_runningtime_NULL <- Scan.interp_MIC_timeend - Scan.interp_MIC_timestart

Scan.interp_Lasso_timestart <- Sys.time()
Scan.interp_Lasso_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Lasso")
Scan.interp_Lasso_timeend <- Sys.time()
Scan.interp_Lasso_runningtime_NULL <- Scan.interp_Lasso_timeend - Scan.interp_Lasso_timestart

Scan.interp_Elastic_timestart <- Sys.time()
Scan.interp_Elastic_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Elastic")
Scan.interp_Elastic_timeend <- Sys.time()
Scan.interp_Elastic_runningtime_NULL <- Scan.interp_Elastic_timeend - Scan.interp_Elastic_timestart

Scan.interp_Ridge_timestart <- Sys.time()
Scan.interp_Ridge_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Ridge")
Scan.interp_Ridge_timeend <- Sys.time()
Scan.interp_Ridge_runningtime_NULL <- Scan.interp_Ridge_timeend - Scan.interp_Ridge_timestart

Scan.interp_Phit_timestart <- Sys.time()
Scan.interp_Phit_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Phit")
Scan.interp_Phit_timeend <- Sys.time()
Scan.interp_Phit_runningtime_NULL <- Scan.interp_Phit_timeend - Scan.interp_Phit_timestart

Scan.interp_Phis_timestart <- Sys.time()
Scan.interp_Phis_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Phis")
Scan.interp_Phis_timeend <- Sys.time()
Scan.interp_Phis_runningtime_NULL <- Scan.interp_Phis_timeend - Scan.interp_Phis_timestart

Scan.interp_Rhop_timestart <- Sys.time()
Scan.interp_Rhop_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Rhop")
Scan.interp_Rhop_timeend <- Sys.time()
Scan.interp_Rhop_runningtime_NULL <- Scan.interp_Rhop_timeend - Scan.interp_Rhop_timestart

Scan.interp_IDA_timestart <- Sys.time()
Scan.interp_IDA_NULL_res <- Scan.interp(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "IDA", pcmethod = "stable")
Scan.interp_IDA_timeend <- Sys.time()
Scan.interp_IDA_runningtime_NULL <- Scan.interp_IDA_timeend - Scan.interp_IDA_timestart


# The prior information is TargetScan
Scan.interp_Pearson_TargetScan_res <- lapply(seq(Scan.interp_Pearson_NULL_res), function(i) Scan.interp_Pearson_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Spearman_TargetScan_res <- lapply(seq(Scan.interp_Spearman_NULL_res), function(i) Scan.interp_Spearman_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Kendall_TargetScan_res <- lapply(seq(Scan.interp_Kendall_NULL_res), function(i) Scan.interp_Kendall_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Dcor_TargetScan_res <- lapply(seq(Scan.interp_Dcor_NULL_res), function(i) Scan.interp_Dcor_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_RDC_TargetScan_res <- lapply(seq(Scan.interp_RDC_NULL_res), function(i) Scan.interp_RDC_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Hoeffd_TargetScan_res <- lapply(seq(Scan.interp_Hoeffd_NULL_res), function(i) Scan.interp_Hoeffd_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Zscore_TargetScan_res <- lapply(seq(Scan.interp_Zscore_NULL_res), function(i) Scan.interp_Zscore_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Biweight_TargetScan_res <- lapply(seq(Scan.interp_Biweight_NULL_res), function(i) Scan.interp_Biweight_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Weighted_rank_TargetScan_res <- lapply(seq(Scan.interp_Weighted_rank_NULL_res), function(i) Scan.interp_Weighted_rank_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Cosine_TargetScan_res <- lapply(seq(Scan.interp_Cosine_NULL_res), function(i) Scan.interp_Cosine_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Euclidean_TargetScan_res <- lapply(seq(Scan.interp_Euclidean_NULL_res), function(i) Scan.interp_Euclidean_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Manhattan_TargetScan_res <- lapply(seq(Scan.interp_Manhattan_NULL_res), function(i) Scan.interp_Manhattan_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Canberra_TargetScan_res <- lapply(seq(Scan.interp_Canberra_NULL_res), function(i) Scan.interp_Canberra_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Chebyshev_TargetScan_res <- lapply(seq(Scan.interp_Chebyshev_NULL_res), function(i) Scan.interp_Chebyshev_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Dice_TargetScan_res <- lapply(seq(Scan.interp_Dice_NULL_res), function(i) Scan.interp_Dice_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Jaccard_TargetScan_res <- lapply(seq(Scan.interp_Jaccard_NULL_res), function(i) Scan.interp_Jaccard_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Mahalanobis_TargetScan_res <- lapply(seq(Scan.interp_Mahalanobis_NULL_res), function(i) Scan.interp_Mahalanobis_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_MI_TargetScan_res <- lapply(seq(Scan.interp_MI_NULL_res), function(i) Scan.interp_MI_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_MIC_TargetScan_res <- lapply(seq(Scan.interp_MIC_NULL_res), function(i) Scan.interp_MIC_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Lasso_TargetScan_res <- lapply(seq(Scan.interp_Lasso_NULL_res), function(i) Scan.interp_Lasso_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Elastic_TargetScan_res <- lapply(seq(Scan.interp_Elastic_NULL_res), function(i) Scan.interp_Elastic_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Ridge_TargetScan_res <- lapply(seq(Scan.interp_Ridge_NULL_res), function(i) Scan.interp_Ridge_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Phit_TargetScan_res <- lapply(seq(Scan.interp_Phit_NULL_res), function(i) Scan.interp_Phit_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Phis_TargetScan_res <- lapply(seq(Scan.interp_Phis_NULL_res), function(i) Scan.interp_Phis_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_Rhop_TargetScan_res <- lapply(seq(Scan.interp_Rhop_NULL_res), function(i) Scan.interp_Rhop_NULL_res[[i]] %s% TargetScan_graph)
Scan.interp_IDA_TargetScan_res <- lapply(seq(Scan.interp_IDA_NULL_res), function(i) Scan.interp_IDA_NULL_res[[i]] %s% TargetScan_graph)

# The prior information is ENCORI
Scan.interp_Pearson_ENCORI_res <- lapply(seq(Scan.interp_Pearson_NULL_res), function(i) Scan.interp_Pearson_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Spearman_ENCORI_res <- lapply(seq(Scan.interp_Spearman_NULL_res), function(i) Scan.interp_Spearman_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Kendall_ENCORI_res <- lapply(seq(Scan.interp_Kendall_NULL_res), function(i) Scan.interp_Kendall_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Dcor_ENCORI_res <- lapply(seq(Scan.interp_Dcor_NULL_res), function(i) Scan.interp_Dcor_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_RDC_ENCORI_res <- lapply(seq(Scan.interp_RDC_NULL_res), function(i) Scan.interp_RDC_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Hoeffd_ENCORI_res <- lapply(seq(Scan.interp_Hoeffd_NULL_res), function(i) Scan.interp_Hoeffd_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Zscore_ENCORI_res <- lapply(seq(Scan.interp_Zscore_NULL_res), function(i) Scan.interp_Zscore_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Biweight_ENCORI_res <- lapply(seq(Scan.interp_Biweight_NULL_res), function(i) Scan.interp_Biweight_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Weighted_rank_ENCORI_res <- lapply(seq(Scan.interp_Weighted_rank_NULL_res), function(i) Scan.interp_Weighted_rank_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Cosine_ENCORI_res <- lapply(seq(Scan.interp_Cosine_NULL_res), function(i) Scan.interp_Cosine_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Euclidean_ENCORI_res <- lapply(seq(Scan.interp_Euclidean_NULL_res), function(i) Scan.interp_Euclidean_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Manhattan_ENCORI_res <- lapply(seq(Scan.interp_Manhattan_NULL_res), function(i) Scan.interp_Manhattan_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Canberra_ENCORI_res <- lapply(seq(Scan.interp_Canberra_NULL_res), function(i) Scan.interp_Canberra_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Chebyshev_ENCORI_res <- lapply(seq(Scan.interp_Chebyshev_NULL_res), function(i) Scan.interp_Chebyshev_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Dice_ENCORI_res <- lapply(seq(Scan.interp_Dice_NULL_res), function(i) Scan.interp_Dice_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Jaccard_ENCORI_res <- lapply(seq(Scan.interp_Jaccard_NULL_res), function(i) Scan.interp_Jaccard_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Mahalanobis_ENCORI_res <- lapply(seq(Scan.interp_Mahalanobis_NULL_res), function(i) Scan.interp_Mahalanobis_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_MI_ENCORI_res <- lapply(seq(Scan.interp_MI_NULL_res), function(i) Scan.interp_MI_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_MIC_ENCORI_res <- lapply(seq(Scan.interp_MIC_NULL_res), function(i) Scan.interp_MIC_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Lasso_ENCORI_res <- lapply(seq(Scan.interp_Lasso_NULL_res), function(i) Scan.interp_Lasso_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Elastic_ENCORI_res <- lapply(seq(Scan.interp_Elastic_NULL_res), function(i) Scan.interp_Elastic_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Ridge_ENCORI_res <- lapply(seq(Scan.interp_Ridge_NULL_res), function(i) Scan.interp_Ridge_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Phit_ENCORI_res <- lapply(seq(Scan.interp_Phit_NULL_res), function(i) Scan.interp_Phit_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Phis_ENCORI_res <- lapply(seq(Scan.interp_Phis_NULL_res), function(i) Scan.interp_Phis_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_Rhop_ENCORI_res <- lapply(seq(Scan.interp_Rhop_NULL_res), function(i) Scan.interp_Rhop_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_IDA_ENCORI_res <- lapply(seq(Scan.interp_IDA_NULL_res), function(i) Scan.interp_IDA_NULL_res[[i]] %s% ENCORI_graph)


# Number of predicted sample-specific interactions
Scan.interp_Pearson_NULL_res_num <- unlist(lapply(seq(Scan.interp_Pearson_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Pearson_NULL_res[[i]] ))))
Scan.interp_Pearson_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Pearson_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Pearson_TargetScan_res[[i]] ))))
Scan.interp_Pearson_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Pearson_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Pearson_ENCORI_res[[i]] ))))

Scan.interp_Spearman_NULL_res_num <- unlist(lapply(seq(Scan.interp_Spearman_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Spearman_NULL_res[[i]] ))))
Scan.interp_Spearman_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Spearman_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Spearman_TargetScan_res[[i]] ))))
Scan.interp_Spearman_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Spearman_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Spearman_ENCORI_res[[i]] ))))

Scan.interp_Kendall_NULL_res_num <- unlist(lapply(seq(Scan.interp_Kendall_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Kendall_NULL_res[[i]] ))))
Scan.interp_Kendall_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Kendall_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Kendall_TargetScan_res[[i]] ))))
Scan.interp_Kendall_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Kendall_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Kendall_ENCORI_res[[i]] ))))

Scan.interp_Dcor_NULL_res_num <- unlist(lapply(seq(Scan.interp_Dcor_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Dcor_NULL_res[[i]] ))))
Scan.interp_Dcor_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Dcor_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Dcor_TargetScan_res[[i]] ))))
Scan.interp_Dcor_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Dcor_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Dcor_ENCORI_res[[i]] ))))

Scan.interp_RDC_NULL_res_num <- unlist(lapply(seq(Scan.interp_RDC_NULL_res), function(i) nrow(as_data_frame(Scan.interp_RDC_NULL_res[[i]] ))))
Scan.interp_RDC_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_RDC_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_RDC_TargetScan_res[[i]] ))))
Scan.interp_RDC_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_RDC_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_RDC_ENCORI_res[[i]] ))))

Scan.interp_Hoeffd_NULL_res_num <- unlist(lapply(seq(Scan.interp_Hoeffd_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_NULL_res[[i]] ))))
Scan.interp_Hoeffd_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Hoeffd_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_TargetScan_res[[i]] ))))
Scan.interp_Hoeffd_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Hoeffd_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_ENCORI_res[[i]] ))))

Scan.interp_Zscore_NULL_res_num <- unlist(lapply(seq(Scan.interp_Zscore_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Zscore_NULL_res[[i]] ))))
Scan.interp_Zscore_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Zscore_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Zscore_TargetScan_res[[i]] ))))
Scan.interp_Zscore_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Zscore_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Zscore_ENCORI_res[[i]] ))))

Scan.interp_Biweight_NULL_res_num <- unlist(lapply(seq(Scan.interp_Biweight_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Biweight_NULL_res[[i]] ))))
Scan.interp_Biweight_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Biweight_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Biweight_TargetScan_res[[i]] ))))
Scan.interp_Biweight_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Biweight_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Biweight_ENCORI_res[[i]] ))))

Scan.interp_Weighted_rank_NULL_res_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_NULL_res[[i]] ))))
Scan.interp_Weighted_rank_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_TargetScan_res[[i]] ))))
Scan.interp_Weighted_rank_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_ENCORI_res[[i]] ))))

Scan.interp_Cosine_NULL_res_num <- unlist(lapply(seq(Scan.interp_Cosine_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Cosine_NULL_res[[i]] ))))
Scan.interp_Cosine_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Cosine_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Cosine_TargetScan_res[[i]] ))))
Scan.interp_Cosine_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Cosine_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Cosine_ENCORI_res[[i]] ))))

Scan.interp_Euclidean_NULL_res_num <- unlist(lapply(seq(Scan.interp_Euclidean_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Euclidean_NULL_res[[i]] ))))
Scan.interp_Euclidean_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Euclidean_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Euclidean_TargetScan_res[[i]] ))))
Scan.interp_Euclidean_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Euclidean_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Euclidean_ENCORI_res[[i]] ))))

Scan.interp_Manhattan_NULL_res_num <- unlist(lapply(seq(Scan.interp_Manhattan_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Manhattan_NULL_res[[i]] ))))
Scan.interp_Manhattan_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Manhattan_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Manhattan_TargetScan_res[[i]] ))))
Scan.interp_Manhattan_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Manhattan_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Manhattan_ENCORI_res[[i]] ))))

Scan.interp_Canberra_NULL_res_num <- unlist(lapply(seq(Scan.interp_Canberra_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Canberra_NULL_res[[i]] ))))
Scan.interp_Canberra_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Canberra_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Canberra_TargetScan_res[[i]] ))))
Scan.interp_Canberra_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Canberra_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Canberra_ENCORI_res[[i]] ))))

Scan.interp_Chebyshev_NULL_res_num <- unlist(lapply(seq(Scan.interp_Chebyshev_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_NULL_res[[i]] ))))
Scan.interp_Chebyshev_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Chebyshev_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_TargetScan_res[[i]] ))))
Scan.interp_Chebyshev_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Chebyshev_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_ENCORI_res[[i]] ))))

Scan.interp_Dice_NULL_res_num <- unlist(lapply(seq(Scan.interp_Dice_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Dice_NULL_res[[i]] ))))
Scan.interp_Dice_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Dice_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Dice_TargetScan_res[[i]] ))))
Scan.interp_Dice_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Dice_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Dice_ENCORI_res[[i]] ))))

Scan.interp_Jaccard_NULL_res_num <- unlist(lapply(seq(Scan.interp_Jaccard_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Jaccard_NULL_res[[i]] ))))
Scan.interp_Jaccard_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Jaccard_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Jaccard_TargetScan_res[[i]] ))))
Scan.interp_Jaccard_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Jaccard_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Jaccard_ENCORI_res[[i]] ))))

Scan.interp_Mahalanobis_NULL_res_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_NULL_res[[i]] ))))
Scan.interp_Mahalanobis_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_TargetScan_res[[i]] ))))
Scan.interp_Mahalanobis_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_ENCORI_res[[i]] ))))

Scan.interp_MI_NULL_res_num <- unlist(lapply(seq(Scan.interp_MI_NULL_res), function(i) nrow(as_data_frame(Scan.interp_MI_NULL_res[[i]] ))))
Scan.interp_MI_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_MI_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_MI_TargetScan_res[[i]] ))))
Scan.interp_MI_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_MI_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_MI_ENCORI_res[[i]] ))))

Scan.interp_MIC_NULL_res_num <- unlist(lapply(seq(Scan.interp_MIC_NULL_res), function(i) nrow(as_data_frame(Scan.interp_MIC_NULL_res[[i]] ))))
Scan.interp_MIC_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_MIC_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_MIC_TargetScan_res[[i]] ))))
Scan.interp_MIC_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_MIC_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_MIC_ENCORI_res[[i]] ))))

Scan.interp_Lasso_NULL_res_num <- unlist(lapply(seq(Scan.interp_Lasso_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Lasso_NULL_res[[i]] ))))
Scan.interp_Lasso_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Lasso_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Lasso_TargetScan_res[[i]] ))))
Scan.interp_Lasso_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Lasso_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Lasso_ENCORI_res[[i]] ))))

Scan.interp_Elastic_NULL_res_num <- unlist(lapply(seq(Scan.interp_Elastic_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Elastic_NULL_res[[i]] ))))
Scan.interp_Elastic_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Elastic_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Elastic_TargetScan_res[[i]] ))))
Scan.interp_Elastic_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Elastic_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Elastic_ENCORI_res[[i]] ))))

Scan.interp_Ridge_NULL_res_num <- unlist(lapply(seq(Scan.interp_Ridge_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Ridge_NULL_res[[i]] ))))
Scan.interp_Ridge_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Ridge_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Ridge_TargetScan_res[[i]] ))))
Scan.interp_Ridge_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Ridge_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Ridge_ENCORI_res[[i]] ))))

Scan.interp_Phit_NULL_res_num <- unlist(lapply(seq(Scan.interp_Phit_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Phit_NULL_res[[i]] ))))
Scan.interp_Phit_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Phit_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Phit_TargetScan_res[[i]] ))))
Scan.interp_Phit_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Phit_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Phit_ENCORI_res[[i]] ))))

Scan.interp_Phis_NULL_res_num <- unlist(lapply(seq(Scan.interp_Phis_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Phis_NULL_res[[i]] ))))
Scan.interp_Phis_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Phis_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Phis_TargetScan_res[[i]] ))))
Scan.interp_Phis_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Phis_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Phis_ENCORI_res[[i]] ))))

Scan.interp_Rhop_NULL_res_num <- unlist(lapply(seq(Scan.interp_Rhop_NULL_res), function(i) nrow(as_data_frame(Scan.interp_Rhop_NULL_res[[i]] ))))
Scan.interp_Rhop_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_Rhop_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_Rhop_TargetScan_res[[i]] ))))
Scan.interp_Rhop_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_Rhop_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_Rhop_ENCORI_res[[i]] ))))

Scan.interp_IDA_NULL_res_num <- unlist(lapply(seq(Scan.interp_IDA_NULL_res), function(i) nrow(as_data_frame(Scan.interp_IDA_NULL_res[[i]] ))))
Scan.interp_IDA_TargetScan_res_num <- unlist(lapply(seq(Scan.interp_IDA_TargetScan_res), function(i) nrow(as_data_frame(Scan.interp_IDA_TargetScan_res[[i]] ))))
Scan.interp_IDA_ENCORI_res_num <- unlist(lapply(seq(Scan.interp_IDA_ENCORI_res), function(i) nrow(as_data_frame(Scan.interp_IDA_ENCORI_res[[i]] ))))


# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.interp_Pearson_NULL_res_validated <- lapply(seq(Scan.interp_Pearson_NULL_res), function(i) as_data_frame(Scan.interp_Pearson_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Pearson_TargetScan_res_validated <- lapply(seq(Scan.interp_Pearson_TargetScan_res), function(i) as_data_frame(Scan.interp_Pearson_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Pearson_ENCORI_res_validated <- lapply(seq(Scan.interp_Pearson_ENCORI_res), function(i) as_data_frame(Scan.interp_Pearson_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Spearman_NULL_res_validated <- lapply(seq(Scan.interp_Spearman_NULL_res), function(i) as_data_frame(Scan.interp_Spearman_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Spearman_TargetScan_res_validated <- lapply(seq(Scan.interp_Spearman_TargetScan_res), function(i) as_data_frame(Scan.interp_Spearman_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Spearman_ENCORI_res_validated <- lapply(seq(Scan.interp_Spearman_ENCORI_res), function(i) as_data_frame(Scan.interp_Spearman_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Kendall_NULL_res_validated <- lapply(seq(Scan.interp_Kendall_NULL_res), function(i) as_data_frame(Scan.interp_Kendall_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Kendall_TargetScan_res_validated <- lapply(seq(Scan.interp_Kendall_TargetScan_res), function(i) as_data_frame(Scan.interp_Kendall_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Kendall_ENCORI_res_validated <- lapply(seq(Scan.interp_Kendall_ENCORI_res), function(i) as_data_frame(Scan.interp_Kendall_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Dcor_NULL_res_validated <- lapply(seq(Scan.interp_Dcor_NULL_res), function(i) as_data_frame(Scan.interp_Dcor_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dcor_TargetScan_res_validated <- lapply(seq(Scan.interp_Dcor_TargetScan_res), function(i) as_data_frame(Scan.interp_Dcor_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dcor_ENCORI_res_validated <- lapply(seq(Scan.interp_Dcor_ENCORI_res), function(i) as_data_frame(Scan.interp_Dcor_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_RDC_NULL_res_validated <- lapply(seq(Scan.interp_RDC_NULL_res), function(i) as_data_frame(Scan.interp_RDC_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_RDC_TargetScan_res_validated <- lapply(seq(Scan.interp_RDC_TargetScan_res), function(i) as_data_frame(Scan.interp_RDC_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_RDC_ENCORI_res_validated <- lapply(seq(Scan.interp_RDC_ENCORI_res), function(i) as_data_frame(Scan.interp_RDC_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Hoeffd_NULL_res_validated <- lapply(seq(Scan.interp_Hoeffd_NULL_res), function(i) as_data_frame(Scan.interp_Hoeffd_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Hoeffd_TargetScan_res_validated <- lapply(seq(Scan.interp_Hoeffd_TargetScan_res), function(i) as_data_frame(Scan.interp_Hoeffd_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Hoeffd_ENCORI_res_validated <- lapply(seq(Scan.interp_Hoeffd_ENCORI_res), function(i) as_data_frame(Scan.interp_Hoeffd_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Zscore_NULL_res_validated <- lapply(seq(Scan.interp_Zscore_NULL_res), function(i) as_data_frame(Scan.interp_Zscore_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Zscore_TargetScan_res_validated <- lapply(seq(Scan.interp_Zscore_TargetScan_res), function(i) as_data_frame(Scan.interp_Zscore_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Zscore_ENCORI_res_validated <- lapply(seq(Scan.interp_Zscore_ENCORI_res), function(i) as_data_frame(Scan.interp_Zscore_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Biweight_NULL_res_validated <- lapply(seq(Scan.interp_Biweight_NULL_res), function(i) as_data_frame(Scan.interp_Biweight_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Biweight_TargetScan_res_validated <- lapply(seq(Scan.interp_Biweight_TargetScan_res), function(i) as_data_frame(Scan.interp_Biweight_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Biweight_ENCORI_res_validated <- lapply(seq(Scan.interp_Biweight_ENCORI_res), function(i) as_data_frame(Scan.interp_Biweight_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Weighted_rank_NULL_res_validated <- lapply(seq(Scan.interp_Weighted_rank_NULL_res), function(i) as_data_frame(Scan.interp_Weighted_rank_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Weighted_rank_TargetScan_res_validated <- lapply(seq(Scan.interp_Weighted_rank_TargetScan_res), function(i) as_data_frame(Scan.interp_Weighted_rank_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Weighted_rank_ENCORI_res_validated <- lapply(seq(Scan.interp_Weighted_rank_ENCORI_res), function(i) as_data_frame(Scan.interp_Weighted_rank_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Cosine_NULL_res_validated <- lapply(seq(Scan.interp_Cosine_NULL_res), function(i) as_data_frame(Scan.interp_Cosine_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Cosine_TargetScan_res_validated <- lapply(seq(Scan.interp_Cosine_TargetScan_res), function(i) as_data_frame(Scan.interp_Cosine_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Cosine_ENCORI_res_validated <- lapply(seq(Scan.interp_Cosine_ENCORI_res), function(i) as_data_frame(Scan.interp_Cosine_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Euclidean_NULL_res_validated <- lapply(seq(Scan.interp_Euclidean_NULL_res), function(i) as_data_frame(Scan.interp_Euclidean_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Euclidean_TargetScan_res_validated <- lapply(seq(Scan.interp_Euclidean_TargetScan_res), function(i) as_data_frame(Scan.interp_Euclidean_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Euclidean_ENCORI_res_validated <- lapply(seq(Scan.interp_Euclidean_ENCORI_res), function(i) as_data_frame(Scan.interp_Euclidean_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Manhattan_NULL_res_validated <- lapply(seq(Scan.interp_Manhattan_NULL_res), function(i) as_data_frame(Scan.interp_Manhattan_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Manhattan_TargetScan_res_validated <- lapply(seq(Scan.interp_Manhattan_TargetScan_res), function(i) as_data_frame(Scan.interp_Manhattan_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Manhattan_ENCORI_res_validated <- lapply(seq(Scan.interp_Manhattan_ENCORI_res), function(i) as_data_frame(Scan.interp_Manhattan_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Canberra_NULL_res_validated <- lapply(seq(Scan.interp_Canberra_NULL_res), function(i) as_data_frame(Scan.interp_Canberra_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Canberra_TargetScan_res_validated <- lapply(seq(Scan.interp_Canberra_TargetScan_res), function(i) as_data_frame(Scan.interp_Canberra_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Canberra_ENCORI_res_validated <- lapply(seq(Scan.interp_Canberra_ENCORI_res), function(i) as_data_frame(Scan.interp_Canberra_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Chebyshev_NULL_res_validated <- lapply(seq(Scan.interp_Chebyshev_NULL_res), function(i) as_data_frame(Scan.interp_Chebyshev_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Chebyshev_TargetScan_res_validated <- lapply(seq(Scan.interp_Chebyshev_TargetScan_res), function(i) as_data_frame(Scan.interp_Chebyshev_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Chebyshev_ENCORI_res_validated <- lapply(seq(Scan.interp_Chebyshev_ENCORI_res), function(i) as_data_frame(Scan.interp_Chebyshev_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Dice_NULL_res_validated <- lapply(seq(Scan.interp_Dice_NULL_res), function(i) as_data_frame(Scan.interp_Dice_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dice_TargetScan_res_validated <- lapply(seq(Scan.interp_Dice_TargetScan_res), function(i) as_data_frame(Scan.interp_Dice_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dice_ENCORI_res_validated <- lapply(seq(Scan.interp_Dice_ENCORI_res), function(i) as_data_frame(Scan.interp_Dice_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Jaccard_NULL_res_validated <- lapply(seq(Scan.interp_Jaccard_NULL_res), function(i) as_data_frame(Scan.interp_Jaccard_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Jaccard_TargetScan_res_validated <- lapply(seq(Scan.interp_Jaccard_TargetScan_res), function(i) as_data_frame(Scan.interp_Jaccard_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Jaccard_ENCORI_res_validated <- lapply(seq(Scan.interp_Jaccard_ENCORI_res), function(i) as_data_frame(Scan.interp_Jaccard_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Mahalanobis_NULL_res_validated <- lapply(seq(Scan.interp_Mahalanobis_NULL_res), function(i) as_data_frame(Scan.interp_Mahalanobis_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Mahalanobis_TargetScan_res_validated <- lapply(seq(Scan.interp_Mahalanobis_TargetScan_res), function(i) as_data_frame(Scan.interp_Mahalanobis_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Mahalanobis_ENCORI_res_validated <- lapply(seq(Scan.interp_Mahalanobis_ENCORI_res), function(i) as_data_frame(Scan.interp_Mahalanobis_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_MI_NULL_res_validated <- lapply(seq(Scan.interp_MI_NULL_res), function(i) as_data_frame(Scan.interp_MI_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MI_TargetScan_res_validated <- lapply(seq(Scan.interp_MI_TargetScan_res), function(i) as_data_frame(Scan.interp_MI_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MI_ENCORI_res_validated <- lapply(seq(Scan.interp_MI_ENCORI_res), function(i) as_data_frame(Scan.interp_MI_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_MIC_NULL_res_validated <- lapply(seq(Scan.interp_MIC_NULL_res), function(i) as_data_frame(Scan.interp_MIC_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MIC_TargetScan_res_validated <- lapply(seq(Scan.interp_MIC_TargetScan_res), function(i) as_data_frame(Scan.interp_MIC_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MIC_ENCORI_res_validated <- lapply(seq(Scan.interp_MIC_ENCORI_res), function(i) as_data_frame(Scan.interp_MIC_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Lasso_NULL_res_validated <- lapply(seq(Scan.interp_Lasso_NULL_res), function(i) as_data_frame(Scan.interp_Lasso_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Lasso_TargetScan_res_validated <- lapply(seq(Scan.interp_Lasso_TargetScan_res), function(i) as_data_frame(Scan.interp_Lasso_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Lasso_ENCORI_res_validated <- lapply(seq(Scan.interp_Lasso_ENCORI_res), function(i) as_data_frame(Scan.interp_Lasso_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Elastic_NULL_res_validated <- lapply(seq(Scan.interp_Elastic_NULL_res), function(i) as_data_frame(Scan.interp_Elastic_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Elastic_TargetScan_res_validated <- lapply(seq(Scan.interp_Elastic_TargetScan_res), function(i) as_data_frame(Scan.interp_Elastic_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Elastic_ENCORI_res_validated <- lapply(seq(Scan.interp_Elastic_ENCORI_res), function(i) as_data_frame(Scan.interp_Elastic_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Ridge_NULL_res_validated <- lapply(seq(Scan.interp_Ridge_NULL_res), function(i) as_data_frame(Scan.interp_Ridge_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Ridge_TargetScan_res_validated <- lapply(seq(Scan.interp_Ridge_TargetScan_res), function(i) as_data_frame(Scan.interp_Ridge_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Ridge_ENCORI_res_validated <- lapply(seq(Scan.interp_Ridge_ENCORI_res), function(i) as_data_frame(Scan.interp_Ridge_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Phit_NULL_res_validated <- lapply(seq(Scan.interp_Phit_NULL_res), function(i) as_data_frame(Scan.interp_Phit_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phit_TargetScan_res_validated <- lapply(seq(Scan.interp_Phit_TargetScan_res), function(i) as_data_frame(Scan.interp_Phit_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phit_ENCORI_res_validated <- lapply(seq(Scan.interp_Phit_ENCORI_res), function(i) as_data_frame(Scan.interp_Phit_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Phis_NULL_res_validated <- lapply(seq(Scan.interp_Phis_NULL_res), function(i) as_data_frame(Scan.interp_Phis_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phis_TargetScan_res_validated <- lapply(seq(Scan.interp_Phis_TargetScan_res), function(i) as_data_frame(Scan.interp_Phis_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phis_ENCORI_res_validated <- lapply(seq(Scan.interp_Phis_ENCORI_res), function(i) as_data_frame(Scan.interp_Phis_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Rhop_NULL_res_validated <- lapply(seq(Scan.interp_Rhop_NULL_res), function(i) as_data_frame(Scan.interp_Rhop_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Rhop_TargetScan_res_validated <- lapply(seq(Scan.interp_Rhop_TargetScan_res), function(i) as_data_frame(Scan.interp_Rhop_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Rhop_ENCORI_res_validated <- lapply(seq(Scan.interp_Rhop_ENCORI_res), function(i) as_data_frame(Scan.interp_Rhop_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_IDA_NULL_res_validated <- lapply(seq(Scan.interp_IDA_NULL_res), function(i) as_data_frame(Scan.interp_IDA_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_IDA_TargetScan_res_validated <- lapply(seq(Scan.interp_IDA_TargetScan_res), function(i) as_data_frame(Scan.interp_IDA_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_IDA_ENCORI_res_validated <- lapply(seq(Scan.interp_IDA_ENCORI_res), function(i) as_data_frame(Scan.interp_IDA_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.interp_Pearson_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_NULL_res), function(i) 100*nrow(Scan.interp_Pearson_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_NULL_res[[i]]))))
Scan.interp_Pearson_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_TargetScan_res), function(i) 100*nrow(Scan.interp_Pearson_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_TargetScan_res[[i]]))))
Scan.interp_Pearson_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_ENCORI_res), function(i) 100*nrow(Scan.interp_Pearson_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_ENCORI_res[[i]]))))

Scan.interp_Spearman_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_NULL_res), function(i) 100*nrow(Scan.interp_Spearman_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_NULL_res[[i]]))))
Scan.interp_Spearman_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_TargetScan_res), function(i) 100*nrow(Scan.interp_Spearman_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_TargetScan_res[[i]]))))
Scan.interp_Spearman_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_ENCORI_res), function(i) 100*nrow(Scan.interp_Spearman_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_ENCORI_res[[i]]))))

Scan.interp_Kendall_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_NULL_res), function(i) 100*nrow(Scan.interp_Kendall_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_NULL_res[[i]]))))
Scan.interp_Kendall_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_TargetScan_res), function(i) 100*nrow(Scan.interp_Kendall_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_TargetScan_res[[i]]))))
Scan.interp_Kendall_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_ENCORI_res), function(i) 100*nrow(Scan.interp_Kendall_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_ENCORI_res[[i]]))))

Scan.interp_Dcor_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_NULL_res), function(i) 100*nrow(Scan.interp_Dcor_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_NULL_res[[i]]))))
Scan.interp_Dcor_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_TargetScan_res), function(i) 100*nrow(Scan.interp_Dcor_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_TargetScan_res[[i]]))))
Scan.interp_Dcor_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_ENCORI_res), function(i) 100*nrow(Scan.interp_Dcor_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_ENCORI_res[[i]]))))

Scan.interp_RDC_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_RDC_NULL_res), function(i) 100*nrow(Scan.interp_RDC_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_NULL_res[[i]]))))
Scan.interp_RDC_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_RDC_TargetScan_res), function(i) 100*nrow(Scan.interp_RDC_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_TargetScan_res[[i]]))))
Scan.interp_RDC_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_RDC_ENCORI_res), function(i) 100*nrow(Scan.interp_RDC_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_ENCORI_res[[i]]))))

Scan.interp_Hoeffd_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_NULL_res), function(i) 100*nrow(Scan.interp_Hoeffd_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_NULL_res[[i]]))))
Scan.interp_Hoeffd_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_TargetScan_res), function(i) 100*nrow(Scan.interp_Hoeffd_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_TargetScan_res[[i]]))))
Scan.interp_Hoeffd_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_ENCORI_res), function(i) 100*nrow(Scan.interp_Hoeffd_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_ENCORI_res[[i]]))))

Scan.interp_Zscore_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_NULL_res), function(i) 100*nrow(Scan.interp_Zscore_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_NULL_res[[i]]))))
Scan.interp_Zscore_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_TargetScan_res), function(i) 100*nrow(Scan.interp_Zscore_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_TargetScan_res[[i]]))))
Scan.interp_Zscore_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_ENCORI_res), function(i) 100*nrow(Scan.interp_Zscore_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_ENCORI_res[[i]]))))

Scan.interp_Biweight_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_NULL_res), function(i) 100*nrow(Scan.interp_Biweight_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_NULL_res[[i]]))))
Scan.interp_Biweight_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_TargetScan_res), function(i) 100*nrow(Scan.interp_Biweight_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_TargetScan_res[[i]]))))
Scan.interp_Biweight_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_ENCORI_res), function(i) 100*nrow(Scan.interp_Biweight_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_ENCORI_res[[i]]))))

Scan.interp_Weighted_rank_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_NULL_res), function(i) 100*nrow(Scan.interp_Weighted_rank_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_NULL_res[[i]]))))
Scan.interp_Weighted_rank_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_TargetScan_res), function(i) 100*nrow(Scan.interp_Weighted_rank_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_TargetScan_res[[i]]))))
Scan.interp_Weighted_rank_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_ENCORI_res), function(i) 100*nrow(Scan.interp_Weighted_rank_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_ENCORI_res[[i]]))))

Scan.interp_Cosine_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_NULL_res), function(i) 100*nrow(Scan.interp_Cosine_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_NULL_res[[i]]))))
Scan.interp_Cosine_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_TargetScan_res), function(i) 100*nrow(Scan.interp_Cosine_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_TargetScan_res[[i]]))))
Scan.interp_Cosine_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_ENCORI_res), function(i) 100*nrow(Scan.interp_Cosine_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_ENCORI_res[[i]]))))

Scan.interp_Euclidean_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_NULL_res), function(i) 100*nrow(Scan.interp_Euclidean_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_NULL_res[[i]]))))
Scan.interp_Euclidean_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_TargetScan_res), function(i) 100*nrow(Scan.interp_Euclidean_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_TargetScan_res[[i]]))))
Scan.interp_Euclidean_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_ENCORI_res), function(i) 100*nrow(Scan.interp_Euclidean_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_ENCORI_res[[i]]))))

Scan.interp_Manhattan_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_NULL_res), function(i) 100*nrow(Scan.interp_Manhattan_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_NULL_res[[i]]))))
Scan.interp_Manhattan_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_TargetScan_res), function(i) 100*nrow(Scan.interp_Manhattan_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_TargetScan_res[[i]]))))
Scan.interp_Manhattan_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_ENCORI_res), function(i) 100*nrow(Scan.interp_Manhattan_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_ENCORI_res[[i]]))))

Scan.interp_Canberra_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_NULL_res), function(i) 100*nrow(Scan.interp_Canberra_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_NULL_res[[i]]))))
Scan.interp_Canberra_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_TargetScan_res), function(i) 100*nrow(Scan.interp_Canberra_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_TargetScan_res[[i]]))))
Scan.interp_Canberra_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_ENCORI_res), function(i) 100*nrow(Scan.interp_Canberra_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_ENCORI_res[[i]]))))

Scan.interp_Chebyshev_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_NULL_res), function(i) 100*nrow(Scan.interp_Chebyshev_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_NULL_res[[i]]))))
Scan.interp_Chebyshev_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_TargetScan_res), function(i) 100*nrow(Scan.interp_Chebyshev_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_TargetScan_res[[i]]))))
Scan.interp_Chebyshev_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_ENCORI_res), function(i) 100*nrow(Scan.interp_Chebyshev_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_ENCORI_res[[i]]))))

Scan.interp_Dice_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Dice_NULL_res), function(i) 100*nrow(Scan.interp_Dice_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_NULL_res[[i]]))))
Scan.interp_Dice_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Dice_TargetScan_res), function(i) 100*nrow(Scan.interp_Dice_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_TargetScan_res[[i]]))))
Scan.interp_Dice_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Dice_ENCORI_res), function(i) 100*nrow(Scan.interp_Dice_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_ENCORI_res[[i]]))))

Scan.interp_Jaccard_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_NULL_res), function(i) 100*nrow(Scan.interp_Jaccard_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_NULL_res[[i]]))))
Scan.interp_Jaccard_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_TargetScan_res), function(i) 100*nrow(Scan.interp_Jaccard_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_TargetScan_res[[i]]))))
Scan.interp_Jaccard_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_ENCORI_res), function(i) 100*nrow(Scan.interp_Jaccard_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_ENCORI_res[[i]]))))

Scan.interp_Mahalanobis_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_NULL_res), function(i) 100*nrow(Scan.interp_Mahalanobis_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_NULL_res[[i]]))))
Scan.interp_Mahalanobis_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_TargetScan_res), function(i) 100*nrow(Scan.interp_Mahalanobis_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_TargetScan_res[[i]]))))
Scan.interp_Mahalanobis_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_ENCORI_res), function(i) 100*nrow(Scan.interp_Mahalanobis_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_ENCORI_res[[i]]))))

Scan.interp_MI_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_MI_NULL_res), function(i) 100*nrow(Scan.interp_MI_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_NULL_res[[i]]))))
Scan.interp_MI_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_MI_TargetScan_res), function(i) 100*nrow(Scan.interp_MI_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_TargetScan_res[[i]]))))
Scan.interp_MI_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_MI_ENCORI_res), function(i) 100*nrow(Scan.interp_MI_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_ENCORI_res[[i]]))))

Scan.interp_MIC_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_MIC_NULL_res), function(i) 100*nrow(Scan.interp_MIC_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_NULL_res[[i]]))))
Scan.interp_MIC_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_MIC_TargetScan_res), function(i) 100*nrow(Scan.interp_MIC_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_TargetScan_res[[i]]))))
Scan.interp_MIC_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_MIC_ENCORI_res), function(i) 100*nrow(Scan.interp_MIC_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_ENCORI_res[[i]]))))

Scan.interp_Lasso_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_NULL_res), function(i) 100*nrow(Scan.interp_Lasso_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_NULL_res[[i]]))))
Scan.interp_Lasso_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_TargetScan_res), function(i) 100*nrow(Scan.interp_Lasso_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_TargetScan_res[[i]]))))
Scan.interp_Lasso_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_ENCORI_res), function(i) 100*nrow(Scan.interp_Lasso_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_ENCORI_res[[i]]))))

Scan.interp_Elastic_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_NULL_res), function(i) 100*nrow(Scan.interp_Elastic_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_NULL_res[[i]]))))
Scan.interp_Elastic_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_TargetScan_res), function(i) 100*nrow(Scan.interp_Elastic_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_TargetScan_res[[i]]))))
Scan.interp_Elastic_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_ENCORI_res), function(i) 100*nrow(Scan.interp_Elastic_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_ENCORI_res[[i]]))))

Scan.interp_Ridge_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_NULL_res), function(i) 100*nrow(Scan.interp_Ridge_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_NULL_res[[i]]))))
Scan.interp_Ridge_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_TargetScan_res), function(i) 100*nrow(Scan.interp_Ridge_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_TargetScan_res[[i]]))))
Scan.interp_Ridge_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_ENCORI_res), function(i) 100*nrow(Scan.interp_Ridge_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_ENCORI_res[[i]]))))

Scan.interp_Phit_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Phit_NULL_res), function(i) 100*nrow(Scan.interp_Phit_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_NULL_res[[i]]))))
Scan.interp_Phit_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Phit_TargetScan_res), function(i) 100*nrow(Scan.interp_Phit_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_TargetScan_res[[i]]))))
Scan.interp_Phit_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Phit_ENCORI_res), function(i) 100*nrow(Scan.interp_Phit_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_ENCORI_res[[i]]))))

Scan.interp_Phis_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Phis_NULL_res), function(i) 100*nrow(Scan.interp_Phis_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_NULL_res[[i]]))))
Scan.interp_Phis_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Phis_TargetScan_res), function(i) 100*nrow(Scan.interp_Phis_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_TargetScan_res[[i]]))))
Scan.interp_Phis_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Phis_ENCORI_res), function(i) 100*nrow(Scan.interp_Phis_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_ENCORI_res[[i]]))))

Scan.interp_Rhop_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_NULL_res), function(i) 100*nrow(Scan.interp_Rhop_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_NULL_res[[i]]))))
Scan.interp_Rhop_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_TargetScan_res), function(i) 100*nrow(Scan.interp_Rhop_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_TargetScan_res[[i]]))))
Scan.interp_Rhop_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_ENCORI_res), function(i) 100*nrow(Scan.interp_Rhop_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_ENCORI_res[[i]]))))

Scan.interp_IDA_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_IDA_NULL_res), function(i) 100*nrow(Scan.interp_IDA_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_NULL_res[[i]]))))
Scan.interp_IDA_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_IDA_TargetScan_res), function(i) 100*nrow(Scan.interp_IDA_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_TargetScan_res[[i]]))))
Scan.interp_IDA_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_IDA_ENCORI_res), function(i) 100*nrow(Scan.interp_IDA_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_ENCORI_res[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.interp_Pearson_NULL_res_Sim <- Sim.network(Scan.interp_Pearson_NULL_res)
Scan.interp_Pearson_TargetScan_res_Sim <- Sim.network(Scan.interp_Pearson_TargetScan_res)
Scan.interp_Pearson_ENCORI_res_Sim <- Sim.network(Scan.interp_Pearson_ENCORI_res)

Scan.interp_Spearman_NULL_res_Sim <- Sim.network(Scan.interp_Spearman_NULL_res)
Scan.interp_Spearman_TargetScan_res_Sim <- Sim.network(Scan.interp_Spearman_TargetScan_res)
Scan.interp_Spearman_ENCORI_res_Sim <- Sim.network(Scan.interp_Spearman_ENCORI_res)

Scan.interp_Kendall_NULL_res_Sim <- Sim.network(Scan.interp_Kendall_NULL_res)
Scan.interp_Kendall_TargetScan_res_Sim <- Sim.network(Scan.interp_Kendall_TargetScan_res)
Scan.interp_Kendall_ENCORI_res_Sim <- Sim.network(Scan.interp_Kendall_ENCORI_res)

Scan.interp_Dcor_NULL_res_Sim <- Sim.network(Scan.interp_Dcor_NULL_res)
Scan.interp_Dcor_TargetScan_res_Sim <- Sim.network(Scan.interp_Dcor_TargetScan_res)
Scan.interp_Dcor_ENCORI_res_Sim <- Sim.network(Scan.interp_Dcor_ENCORI_res)

Scan.interp_RDC_NULL_res_Sim <- Sim.network(Scan.interp_RDC_NULL_res)
Scan.interp_RDC_TargetScan_res_Sim <- Sim.network(Scan.interp_RDC_TargetScan_res)
Scan.interp_RDC_ENCORI_res_Sim <- Sim.network(Scan.interp_RDC_ENCORI_res)

Scan.interp_Hoeffd_NULL_res_Sim <- Sim.network(Scan.interp_Hoeffd_NULL_res)
Scan.interp_Hoeffd_TargetScan_res_Sim <- Sim.network(Scan.interp_Hoeffd_TargetScan_res)
Scan.interp_Hoeffd_ENCORI_res_Sim <- Sim.network(Scan.interp_Hoeffd_ENCORI_res)

Scan.interp_Zscore_NULL_res_Sim <- Sim.network(Scan.interp_Zscore_NULL_res)
Scan.interp_Zscore_TargetScan_res_Sim <- Sim.network(Scan.interp_Zscore_TargetScan_res)
Scan.interp_Zscore_ENCORI_res_Sim <- Sim.network(Scan.interp_Zscore_ENCORI_res)

Scan.interp_Biweight_NULL_res_Sim <- Sim.network(Scan.interp_Biweight_NULL_res)
Scan.interp_Biweight_TargetScan_res_Sim <- Sim.network(Scan.interp_Biweight_TargetScan_res)
Scan.interp_Biweight_ENCORI_res_Sim <- Sim.network(Scan.interp_Biweight_ENCORI_res)

Scan.interp_Weighted_rank_NULL_res_Sim <- Sim.network(Scan.interp_Weighted_rank_NULL_res)
Scan.interp_Weighted_rank_TargetScan_res_Sim <- Sim.network(Scan.interp_Weighted_rank_TargetScan_res)
Scan.interp_Weighted_rank_ENCORI_res_Sim <- Sim.network(Scan.interp_Weighted_rank_ENCORI_res)

Scan.interp_Cosine_NULL_res_Sim <- Sim.network(Scan.interp_Cosine_NULL_res)
Scan.interp_Cosine_TargetScan_res_Sim <- Sim.network(Scan.interp_Cosine_TargetScan_res)
Scan.interp_Cosine_ENCORI_res_Sim <- Sim.network(Scan.interp_Cosine_ENCORI_res)

Scan.interp_Euclidean_NULL_res_Sim <- Sim.network(Scan.interp_Euclidean_NULL_res)
Scan.interp_Euclidean_TargetScan_res_Sim <- Sim.network(Scan.interp_Euclidean_TargetScan_res)
Scan.interp_Euclidean_ENCORI_res_Sim <- Sim.network(Scan.interp_Euclidean_ENCORI_res)

Scan.interp_Manhattan_NULL_res_Sim <- Sim.network(Scan.interp_Manhattan_NULL_res)
Scan.interp_Manhattan_TargetScan_res_Sim <- Sim.network(Scan.interp_Manhattan_TargetScan_res)
Scan.interp_Manhattan_ENCORI_res_Sim <- Sim.network(Scan.interp_Manhattan_ENCORI_res)

Scan.interp_Canberra_NULL_res_Sim <- Sim.network(Scan.interp_Canberra_NULL_res)
Scan.interp_Canberra_TargetScan_res_Sim <- Sim.network(Scan.interp_Canberra_TargetScan_res)
Scan.interp_Canberra_ENCORI_res_Sim <- Sim.network(Scan.interp_Canberra_ENCORI_res)

Scan.interp_Chebyshev_NULL_res_Sim <- Sim.network(Scan.interp_Chebyshev_NULL_res)
Scan.interp_Chebyshev_TargetScan_res_Sim <- Sim.network(Scan.interp_Chebyshev_TargetScan_res)
Scan.interp_Chebyshev_ENCORI_res_Sim <- Sim.network(Scan.interp_Chebyshev_ENCORI_res)

Scan.interp_Dice_NULL_res_Sim <- Sim.network(Scan.interp_Dice_NULL_res)
Scan.interp_Dice_TargetScan_res_Sim <- Sim.network(Scan.interp_Dice_TargetScan_res)
Scan.interp_Dice_ENCORI_res_Sim <- Sim.network(Scan.interp_Dice_ENCORI_res)

Scan.interp_Jaccard_NULL_res_Sim <- Sim.network(Scan.interp_Jaccard_NULL_res)
Scan.interp_Jaccard_TargetScan_res_Sim <- Sim.network(Scan.interp_Jaccard_TargetScan_res)
Scan.interp_Jaccard_ENCORI_res_Sim <- Sim.network(Scan.interp_Jaccard_ENCORI_res)

Scan.interp_Mahalanobis_NULL_res_Sim <- Sim.network(Scan.interp_Mahalanobis_NULL_res)
Scan.interp_Mahalanobis_TargetScan_res_Sim <- Sim.network(Scan.interp_Mahalanobis_TargetScan_res)
Scan.interp_Mahalanobis_ENCORI_res_Sim <- Sim.network(Scan.interp_Mahalanobis_ENCORI_res)

Scan.interp_MI_NULL_res_Sim <- Sim.network(Scan.interp_MI_NULL_res)
Scan.interp_MI_TargetScan_res_Sim <- Sim.network(Scan.interp_MI_TargetScan_res)
Scan.interp_MI_ENCORI_res_Sim <- Sim.network(Scan.interp_MI_ENCORI_res)

Scan.interp_MIC_NULL_res_Sim <- Sim.network(Scan.interp_MIC_NULL_res)
Scan.interp_MIC_TargetScan_res_Sim <- Sim.network(Scan.interp_MIC_TargetScan_res)
Scan.interp_MIC_ENCORI_res_Sim <- Sim.network(Scan.interp_MIC_ENCORI_res)

Scan.interp_Lasso_NULL_res_Sim <- Sim.network(Scan.interp_Lasso_NULL_res)
Scan.interp_Lasso_TargetScan_res_Sim <- Sim.network(Scan.interp_Lasso_TargetScan_res)
Scan.interp_Lasso_ENCORI_res_Sim <- Sim.network(Scan.interp_Lasso_ENCORI_res)

Scan.interp_Elastic_NULL_res_Sim <- Sim.network(Scan.interp_Elastic_NULL_res)
Scan.interp_Elastic_TargetScan_res_Sim <- Sim.network(Scan.interp_Elastic_TargetScan_res)
Scan.interp_Elastic_ENCORI_res_Sim <- Sim.network(Scan.interp_Elastic_ENCORI_res)

Scan.interp_Ridge_NULL_res_Sim <- Sim.network(Scan.interp_Ridge_NULL_res)
Scan.interp_Ridge_TargetScan_res_Sim <- Sim.network(Scan.interp_Ridge_TargetScan_res)
Scan.interp_Ridge_ENCORI_res_Sim <- Sim.network(Scan.interp_Ridge_ENCORI_res)

Scan.interp_Phit_NULL_res_Sim <- Sim.network(Scan.interp_Phit_NULL_res)
Scan.interp_Phit_TargetScan_res_Sim <- Sim.network(Scan.interp_Phit_TargetScan_res)
Scan.interp_Phit_ENCORI_res_Sim <- Sim.network(Scan.interp_Phit_ENCORI_res)

Scan.interp_Phis_NULL_res_Sim <- Sim.network(Scan.interp_Phis_NULL_res)
Scan.interp_Phis_TargetScan_res_Sim <- Sim.network(Scan.interp_Phis_TargetScan_res)
Scan.interp_Phis_ENCORI_res_Sim <- Sim.network(Scan.interp_Phis_ENCORI_res)

Scan.interp_Rhop_NULL_res_Sim <- Sim.network(Scan.interp_Rhop_NULL_res)
Scan.interp_Rhop_TargetScan_res_Sim <- Sim.network(Scan.interp_Rhop_TargetScan_res)
Scan.interp_Rhop_ENCORI_res_Sim <- Sim.network(Scan.interp_Rhop_ENCORI_res)

Scan.interp_IDA_NULL_res_Sim <- Sim.network(Scan.interp_IDA_NULL_res)
Scan.interp_IDA_TargetScan_res_Sim <- Sim.network(Scan.interp_IDA_TargetScan_res)
Scan.interp_IDA_ENCORI_res_Sim <- Sim.network(Scan.interp_IDA_ENCORI_res)

save.image("Scan.interp_K562.RData")

################################################################################# ######################## ############################################################################
########################################################################### Scan.interp application in BRCA dataset ###################################################################
################################################################################# ######################## ############################################################################

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

set.seed(123)
ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

# No prior information
Scan.interp_Pearson_timestart_BRCA <- Sys.time()
Scan.interp_Pearson_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Pearson")
Scan.interp_Pearson_timeend_BRCA <- Sys.time()
Scan.interp_Pearson_runningtime_NULL_BRCA <- Scan.interp_Pearson_timeend_BRCA - Scan.interp_Pearson_timestart_BRCA

Scan.interp_Spearman_timestart_BRCA <- Sys.time()
Scan.interp_Spearman_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Spearman")
Scan.interp_Spearman_timeend_BRCA <- Sys.time()
Scan.interp_Spearman_runningtime_NULL_BRCA <- Scan.interp_Spearman_timeend_BRCA - Scan.interp_Spearman_timestart_BRCA

Scan.interp_Kendall_timestart_BRCA <- Sys.time()
Scan.interp_Kendall_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Kendall")
Scan.interp_Kendall_timeend_BRCA <- Sys.time()
Scan.interp_Kendall_runningtime_NULL_BRCA <- Scan.interp_Kendall_timeend_BRCA - Scan.interp_Kendall_timestart_BRCA

Scan.interp_Dcor_timestart_BRCA <- Sys.time()
Scan.interp_Dcor_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Dcor")
Scan.interp_Dcor_timeend_BRCA <- Sys.time()
Scan.interp_Dcor_runningtime_NULL_BRCA <- Scan.interp_Dcor_timeend_BRCA - Scan.interp_Dcor_timestart_BRCA

Scan.interp_RDC_timestart_BRCA <- Sys.time()
Scan.interp_RDC_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "RDC")
Scan.interp_RDC_timeend_BRCA <- Sys.time()
Scan.interp_RDC_runningtime_NULL_BRCA <- Scan.interp_RDC_timeend_BRCA - Scan.interp_RDC_timestart_BRCA

Scan.interp_Hoeffd_timestart_BRCA <- Sys.time()
Scan.interp_Hoeffd_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Hoeffd")
Scan.interp_Hoeffd_timeend_BRCA <- Sys.time()
Scan.interp_Hoeffd_runningtime_NULL_BRCA <- Scan.interp_Hoeffd_timeend_BRCA - Scan.interp_Hoeffd_timestart_BRCA

Scan.interp_Zscore_timestart_BRCA <- Sys.time()
Scan.interp_Zscore_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Zscore")
Scan.interp_Zscore_timeend_BRCA <- Sys.time()
Scan.interp_Zscore_runningtime_NULL_BRCA <- Scan.interp_Zscore_timeend_BRCA - Scan.interp_Zscore_timestart_BRCA

Scan.interp_Biweight_timestart_BRCA <- Sys.time()
Scan.interp_Biweight_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Biweight")
Scan.interp_Biweight_timeend_BRCA <- Sys.time()
Scan.interp_Biweight_runningtime_NULL_BRCA <- Scan.interp_Biweight_timeend_BRCA - Scan.interp_Biweight_timestart_BRCA

Scan.interp_Weighted_rank_timestart_BRCA <- Sys.time()
Scan.interp_Weighted_rank_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Weighted_rank")
Scan.interp_Weighted_rank_timeend_BRCA <- Sys.time()
Scan.interp_Weighted_rank_runningtime_NULL_BRCA <- Scan.interp_Weighted_rank_timeend_BRCA - Scan.interp_Weighted_rank_timestart_BRCA

Scan.interp_Cosine_timestart_BRCA <- Sys.time()
Scan.interp_Cosine_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Cosine")
Scan.interp_Cosine_timeend_BRCA <- Sys.time()
Scan.interp_Cosine_runningtime_NULL_BRCA <- Scan.interp_Cosine_timeend_BRCA - Scan.interp_Cosine_timestart_BRCA

Scan.interp_Euclidean_timestart_BRCA <- Sys.time()
Scan.interp_Euclidean_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Euclidean")
Scan.interp_Euclidean_timeend_BRCA <- Sys.time()
Scan.interp_Euclidean_runningtime_NULL_BRCA <- Scan.interp_Euclidean_timeend_BRCA - Scan.interp_Euclidean_timestart_BRCA

Scan.interp_Manhattan_timestart_BRCA <- Sys.time()
Scan.interp_Manhattan_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Manhattan")
Scan.interp_Manhattan_timeend_BRCA <- Sys.time()
Scan.interp_Manhattan_runningtime_NULL_BRCA <- Scan.interp_Manhattan_timeend_BRCA - Scan.interp_Manhattan_timestart_BRCA

Scan.interp_Canberra_timestart_BRCA <- Sys.time()
Scan.interp_Canberra_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Canberra")
Scan.interp_Canberra_timeend_BRCA <- Sys.time()
Scan.interp_Canberra_runningtime_NULL_BRCA <- Scan.interp_Canberra_timeend_BRCA - Scan.interp_Canberra_timestart_BRCA

Scan.interp_Chebyshev_timestart_BRCA <- Sys.time()
Scan.interp_Chebyshev_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Chebyshev")
Scan.interp_Chebyshev_timeend_BRCA <- Sys.time()
Scan.interp_Chebyshev_runningtime_NULL_BRCA <- Scan.interp_Chebyshev_timeend_BRCA - Scan.interp_Chebyshev_timestart_BRCA

Scan.interp_Dice_timestart_BRCA <- Sys.time()
Scan.interp_Dice_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Dice")
Scan.interp_Dice_timeend_BRCA <- Sys.time()
Scan.interp_Dice_runningtime_NULL_BRCA <- Scan.interp_Dice_timeend_BRCA - Scan.interp_Dice_timestart_BRCA

Scan.interp_Jaccard_timestart_BRCA <- Sys.time()
Scan.interp_Jaccard_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Jaccard")
Scan.interp_Jaccard_timeend_BRCA <- Sys.time()
Scan.interp_Jaccard_runningtime_NULL_BRCA <- Scan.interp_Jaccard_timeend_BRCA - Scan.interp_Jaccard_timestart_BRCA

Scan.interp_Mahalanobis_timestart_BRCA <- Sys.time()
Scan.interp_Mahalanobis_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Mahalanobis")
Scan.interp_Mahalanobis_timeend_BRCA <- Sys.time()
Scan.interp_Mahalanobis_runningtime_NULL_BRCA <- Scan.interp_Mahalanobis_timeend_BRCA - Scan.interp_Mahalanobis_timestart_BRCA

Scan.interp_MI_timestart_BRCA <- Sys.time()
Scan.interp_MI_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "MI")
Scan.interp_MI_timeend_BRCA <- Sys.time()
Scan.interp_MI_runningtime_NULL_BRCA <- Scan.interp_MI_timeend_BRCA - Scan.interp_MI_timestart_BRCA

Scan.interp_MIC_timestart_BRCA <- Sys.time()
Scan.interp_MIC_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "MIC")
Scan.interp_MIC_timeend_BRCA <- Sys.time()
Scan.interp_MIC_runningtime_NULL_BRCA <- Scan.interp_MIC_timeend_BRCA - Scan.interp_MIC_timestart_BRCA

Scan.interp_Lasso_timestart_BRCA <- Sys.time()
Scan.interp_Lasso_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Lasso")
Scan.interp_Lasso_timeend_BRCA <- Sys.time()
Scan.interp_Lasso_runningtime_NULL_BRCA <- Scan.interp_Lasso_timeend_BRCA - Scan.interp_Lasso_timestart_BRCA

Scan.interp_Elastic_timestart_BRCA <- Sys.time()
Scan.interp_Elastic_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Elastic")
Scan.interp_Elastic_timeend_BRCA <- Sys.time()
Scan.interp_Elastic_runningtime_NULL_BRCA <- Scan.interp_Elastic_timeend_BRCA - Scan.interp_Elastic_timestart_BRCA

Scan.interp_Ridge_timestart_BRCA <- Sys.time()
Scan.interp_Ridge_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Ridge")
Scan.interp_Ridge_timeend_BRCA <- Sys.time()
Scan.interp_Ridge_runningtime_NULL_BRCA <- Scan.interp_Ridge_timeend_BRCA - Scan.interp_Ridge_timestart_BRCA

Scan.interp_Phit_timestart_BRCA <- Sys.time()
Scan.interp_Phit_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Phit")
Scan.interp_Phit_timeend_BRCA <- Sys.time()
Scan.interp_Phit_runningtime_NULL_BRCA <- Scan.interp_Phit_timeend_BRCA - Scan.interp_Phit_timestart_BRCA

Scan.interp_Phis_timestart_BRCA <- Sys.time()
Scan.interp_Phis_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Phis")
Scan.interp_Phis_timeend_BRCA <- Sys.time()
Scan.interp_Phis_runningtime_NULL_BRCA <- Scan.interp_Phis_timeend_BRCA - Scan.interp_Phis_timestart_BRCA

Scan.interp_Rhop_timestart_BRCA <- Sys.time()
Scan.interp_Rhop_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Rhop")
Scan.interp_Rhop_timeend_BRCA <- Sys.time()
Scan.interp_Rhop_runningtime_NULL_BRCA <- Scan.interp_Rhop_timeend_BRCA - Scan.interp_Rhop_timestart_BRCA

Scan.interp_IDA_timestart_BRCA <- Sys.time()
Scan.interp_IDA_NULL_res_BRCA <- Scan.interp(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "IDA", pcmethod = "stable")
Scan.interp_IDA_timeend_BRCA <- Sys.time()
Scan.interp_IDA_runningtime_NULL_BRCA <- Scan.interp_IDA_timeend_BRCA - Scan.interp_IDA_timestart_BRCA

# The prior information is TargetScan
Scan.interp_Pearson_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Pearson_NULL_res_BRCA), function(i) Scan.interp_Pearson_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Spearman_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Spearman_NULL_res_BRCA), function(i) Scan.interp_Spearman_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Kendall_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Kendall_NULL_res_BRCA), function(i) Scan.interp_Kendall_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Dcor_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Dcor_NULL_res_BRCA), function(i) Scan.interp_Dcor_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_RDC_TargetScan_res_BRCA <- lapply(seq(Scan.interp_RDC_NULL_res_BRCA), function(i) Scan.interp_RDC_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Hoeffd_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Hoeffd_NULL_res_BRCA), function(i) Scan.interp_Hoeffd_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Zscore_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Zscore_NULL_res_BRCA), function(i) Scan.interp_Zscore_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Biweight_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Biweight_NULL_res_BRCA), function(i) Scan.interp_Biweight_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Weighted_rank_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Weighted_rank_NULL_res_BRCA), function(i) Scan.interp_Weighted_rank_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Cosine_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Cosine_NULL_res_BRCA), function(i) Scan.interp_Cosine_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Euclidean_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Euclidean_NULL_res_BRCA), function(i) Scan.interp_Euclidean_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Manhattan_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Manhattan_NULL_res_BRCA), function(i) Scan.interp_Manhattan_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Canberra_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Canberra_NULL_res_BRCA), function(i) Scan.interp_Canberra_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Chebyshev_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Chebyshev_NULL_res_BRCA), function(i) Scan.interp_Chebyshev_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Dice_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Dice_NULL_res_BRCA), function(i) Scan.interp_Dice_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Jaccard_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Jaccard_NULL_res_BRCA), function(i) Scan.interp_Jaccard_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Mahalanobis_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Mahalanobis_NULL_res_BRCA), function(i) Scan.interp_Mahalanobis_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_MI_TargetScan_res_BRCA <- lapply(seq(Scan.interp_MI_NULL_res_BRCA), function(i) Scan.interp_MI_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_MIC_TargetScan_res_BRCA <- lapply(seq(Scan.interp_MIC_NULL_res_BRCA), function(i) Scan.interp_MIC_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Lasso_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Lasso_NULL_res_BRCA), function(i) Scan.interp_Lasso_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Elastic_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Elastic_NULL_res_BRCA), function(i) Scan.interp_Elastic_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Ridge_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Ridge_NULL_res_BRCA), function(i) Scan.interp_Ridge_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Phit_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Phit_NULL_res_BRCA), function(i) Scan.interp_Phit_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Phis_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Phis_NULL_res_BRCA), function(i) Scan.interp_Phis_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_Rhop_TargetScan_res_BRCA <- lapply(seq(Scan.interp_Rhop_NULL_res_BRCA), function(i) Scan.interp_Rhop_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.interp_IDA_TargetScan_res_BRCA <- lapply(seq(Scan.interp_IDA_NULL_res_BRCA), function(i) Scan.interp_IDA_NULL_res_BRCA[[i]] %s% TargetScan_graph)

# The prior information is ENCORI
Scan.interp_Pearson_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Pearson_NULL_res_BRCA), function(i) Scan.interp_Pearson_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Spearman_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Spearman_NULL_res_BRCA), function(i) Scan.interp_Spearman_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Kendall_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Kendall_NULL_res_BRCA), function(i) Scan.interp_Kendall_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Dcor_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Dcor_NULL_res_BRCA), function(i) Scan.interp_Dcor_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_RDC_ENCORI_res_BRCA <- lapply(seq(Scan.interp_RDC_NULL_res_BRCA), function(i) Scan.interp_RDC_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Hoeffd_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Hoeffd_NULL_res_BRCA), function(i) Scan.interp_Hoeffd_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Zscore_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Zscore_NULL_res_BRCA), function(i) Scan.interp_Zscore_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Biweight_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Biweight_NULL_res_BRCA), function(i) Scan.interp_Biweight_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Weighted_rank_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Weighted_rank_NULL_res_BRCA), function(i) Scan.interp_Weighted_rank_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Cosine_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Cosine_NULL_res_BRCA), function(i) Scan.interp_Cosine_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Euclidean_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Euclidean_NULL_res_BRCA), function(i) Scan.interp_Euclidean_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Manhattan_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Manhattan_NULL_res_BRCA), function(i) Scan.interp_Manhattan_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Canberra_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Canberra_NULL_res_BRCA), function(i) Scan.interp_Canberra_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Chebyshev_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Chebyshev_NULL_res_BRCA), function(i) Scan.interp_Chebyshev_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Dice_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Dice_NULL_res_BRCA), function(i) Scan.interp_Dice_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Jaccard_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Jaccard_NULL_res_BRCA), function(i) Scan.interp_Jaccard_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Mahalanobis_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Mahalanobis_NULL_res_BRCA), function(i) Scan.interp_Mahalanobis_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_MI_ENCORI_res_BRCA <- lapply(seq(Scan.interp_MI_NULL_res_BRCA), function(i) Scan.interp_MI_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_MIC_ENCORI_res_BRCA <- lapply(seq(Scan.interp_MIC_NULL_res_BRCA), function(i) Scan.interp_MIC_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Lasso_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Lasso_NULL_res_BRCA), function(i) Scan.interp_Lasso_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Elastic_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Elastic_NULL_res_BRCA), function(i) Scan.interp_Elastic_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Ridge_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Ridge_NULL_res_BRCA), function(i) Scan.interp_Ridge_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Phit_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Phit_NULL_res_BRCA), function(i) Scan.interp_Phit_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Phis_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Phis_NULL_res_BRCA), function(i) Scan.interp_Phis_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_Rhop_ENCORI_res_BRCA <- lapply(seq(Scan.interp_Rhop_NULL_res_BRCA), function(i) Scan.interp_Rhop_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_IDA_ENCORI_res_BRCA <- lapply(seq(Scan.interp_IDA_NULL_res_BRCA), function(i) Scan.interp_IDA_NULL_res_BRCA[[i]] %s% ENCORI_graph)

# Number of predicted sample-specific interactions
Scan.interp_Pearson_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Pearson_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Pearson_NULL_res_BRCA[[i]] ))))
Scan.interp_Pearson_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Pearson_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Pearson_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Pearson_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Pearson_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Pearson_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Spearman_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Spearman_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Spearman_NULL_res_BRCA[[i]] ))))
Scan.interp_Spearman_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Spearman_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Spearman_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Spearman_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Spearman_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Spearman_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Kendall_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Kendall_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Kendall_NULL_res_BRCA[[i]] ))))
Scan.interp_Kendall_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Kendall_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Kendall_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Kendall_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Kendall_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Kendall_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Dcor_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dcor_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dcor_NULL_res_BRCA[[i]] ))))
Scan.interp_Dcor_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dcor_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dcor_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Dcor_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dcor_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dcor_ENCORI_res_BRCA[[i]] ))))

Scan.interp_RDC_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_RDC_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_RDC_NULL_res_BRCA[[i]] ))))
Scan.interp_RDC_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_RDC_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_RDC_TargetScan_res_BRCA[[i]] ))))
Scan.interp_RDC_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_RDC_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_RDC_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Hoeffd_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Hoeffd_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_NULL_res_BRCA[[i]] ))))
Scan.interp_Hoeffd_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Hoeffd_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Hoeffd_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Hoeffd_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Hoeffd_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Zscore_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Zscore_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Zscore_NULL_res_BRCA[[i]] ))))
Scan.interp_Zscore_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Zscore_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Zscore_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Zscore_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Zscore_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Zscore_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Biweight_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Biweight_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Biweight_NULL_res_BRCA[[i]] ))))
Scan.interp_Biweight_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Biweight_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Biweight_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Biweight_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Biweight_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Biweight_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Weighted_rank_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_NULL_res_BRCA[[i]] ))))
Scan.interp_Weighted_rank_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Weighted_rank_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Weighted_rank_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Weighted_rank_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Cosine_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Cosine_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Cosine_NULL_res_BRCA[[i]] ))))
Scan.interp_Cosine_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Cosine_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Cosine_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Cosine_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Cosine_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Cosine_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Euclidean_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Euclidean_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Euclidean_NULL_res_BRCA[[i]] ))))
Scan.interp_Euclidean_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Euclidean_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Euclidean_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Euclidean_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Euclidean_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Euclidean_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Manhattan_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Manhattan_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Manhattan_NULL_res_BRCA[[i]] ))))
Scan.interp_Manhattan_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Manhattan_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Manhattan_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Manhattan_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Manhattan_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Manhattan_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Canberra_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Canberra_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Canberra_NULL_res_BRCA[[i]] ))))
Scan.interp_Canberra_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Canberra_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Canberra_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Canberra_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Canberra_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Canberra_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Chebyshev_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Chebyshev_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_NULL_res_BRCA[[i]] ))))
Scan.interp_Chebyshev_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Chebyshev_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Chebyshev_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Chebyshev_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Chebyshev_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Dice_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dice_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dice_NULL_res_BRCA[[i]] ))))
Scan.interp_Dice_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dice_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dice_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Dice_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Dice_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Dice_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Jaccard_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Jaccard_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Jaccard_NULL_res_BRCA[[i]] ))))
Scan.interp_Jaccard_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Jaccard_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Jaccard_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Jaccard_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Jaccard_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Jaccard_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Mahalanobis_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_NULL_res_BRCA[[i]] ))))
Scan.interp_Mahalanobis_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Mahalanobis_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Mahalanobis_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Mahalanobis_ENCORI_res_BRCA[[i]] ))))

Scan.interp_MI_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MI_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MI_NULL_res_BRCA[[i]] ))))
Scan.interp_MI_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MI_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MI_TargetScan_res_BRCA[[i]] ))))
Scan.interp_MI_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MI_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MI_ENCORI_res_BRCA[[i]] ))))

Scan.interp_MIC_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MIC_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MIC_NULL_res_BRCA[[i]] ))))
Scan.interp_MIC_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MIC_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MIC_TargetScan_res_BRCA[[i]] ))))
Scan.interp_MIC_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_MIC_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_MIC_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Lasso_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Lasso_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Lasso_NULL_res_BRCA[[i]] ))))
Scan.interp_Lasso_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Lasso_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Lasso_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Lasso_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Lasso_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Lasso_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Elastic_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Elastic_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Elastic_NULL_res_BRCA[[i]] ))))
Scan.interp_Elastic_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Elastic_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Elastic_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Elastic_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Elastic_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Elastic_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Ridge_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Ridge_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Ridge_NULL_res_BRCA[[i]] ))))
Scan.interp_Ridge_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Ridge_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Ridge_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Ridge_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Ridge_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Ridge_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Phit_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phit_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phit_NULL_res_BRCA[[i]] ))))
Scan.interp_Phit_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phit_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phit_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Phit_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phit_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phit_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Phis_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phis_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phis_NULL_res_BRCA[[i]] ))))
Scan.interp_Phis_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phis_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phis_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Phis_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Phis_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Phis_ENCORI_res_BRCA[[i]] ))))

Scan.interp_Rhop_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Rhop_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Rhop_NULL_res_BRCA[[i]] ))))
Scan.interp_Rhop_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Rhop_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Rhop_TargetScan_res_BRCA[[i]] ))))
Scan.interp_Rhop_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_Rhop_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_Rhop_ENCORI_res_BRCA[[i]] ))))

Scan.interp_IDA_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_IDA_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_IDA_NULL_res_BRCA[[i]] ))))
Scan.interp_IDA_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_IDA_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_IDA_TargetScan_res_BRCA[[i]] ))))
Scan.interp_IDA_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_IDA_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_IDA_ENCORI_res_BRCA[[i]] ))))

# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.interp_Pearson_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Pearson_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Pearson_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Pearson_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Pearson_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Pearson_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Pearson_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Pearson_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Pearson_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Spearman_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Spearman_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Spearman_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Spearman_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Spearman_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Spearman_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Spearman_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Spearman_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Spearman_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Kendall_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Kendall_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Kendall_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Kendall_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Kendall_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Kendall_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Kendall_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Kendall_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Kendall_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Dcor_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Dcor_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Dcor_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dcor_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Dcor_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Dcor_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dcor_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Dcor_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Dcor_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_RDC_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_RDC_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_RDC_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_RDC_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_RDC_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_RDC_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_RDC_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_RDC_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_RDC_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Hoeffd_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Hoeffd_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Hoeffd_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Hoeffd_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Hoeffd_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Hoeffd_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Hoeffd_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Hoeffd_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Hoeffd_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Zscore_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Zscore_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Zscore_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Zscore_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Zscore_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Zscore_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Zscore_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Zscore_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Zscore_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Biweight_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Biweight_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Biweight_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Biweight_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Biweight_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Biweight_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Biweight_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Biweight_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Biweight_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Weighted_rank_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Weighted_rank_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Weighted_rank_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Weighted_rank_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Weighted_rank_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Weighted_rank_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Weighted_rank_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Weighted_rank_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Weighted_rank_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Cosine_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Cosine_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Cosine_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Cosine_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Cosine_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Cosine_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Cosine_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Cosine_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Cosine_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Euclidean_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Euclidean_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Euclidean_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Euclidean_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Euclidean_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Euclidean_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Euclidean_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Euclidean_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Euclidean_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Manhattan_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Manhattan_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Manhattan_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Manhattan_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Manhattan_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Manhattan_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Manhattan_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Manhattan_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Manhattan_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Canberra_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Canberra_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Canberra_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Canberra_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Canberra_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Canberra_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Canberra_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Canberra_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Canberra_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Chebyshev_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Chebyshev_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Chebyshev_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Chebyshev_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Chebyshev_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Chebyshev_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Chebyshev_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Chebyshev_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Chebyshev_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Dice_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Dice_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Dice_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dice_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Dice_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Dice_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Dice_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Dice_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Dice_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Jaccard_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Jaccard_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Jaccard_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Jaccard_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Jaccard_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Jaccard_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Jaccard_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Jaccard_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Jaccard_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Mahalanobis_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Mahalanobis_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Mahalanobis_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Mahalanobis_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Mahalanobis_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Mahalanobis_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Mahalanobis_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Mahalanobis_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Mahalanobis_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_MI_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_MI_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_MI_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MI_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_MI_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_MI_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MI_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_MI_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_MI_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_MIC_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_MIC_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_MIC_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MIC_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_MIC_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_MIC_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_MIC_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_MIC_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_MIC_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Lasso_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Lasso_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Lasso_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Lasso_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Lasso_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Lasso_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Lasso_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Lasso_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Lasso_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Elastic_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Elastic_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Elastic_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Elastic_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Elastic_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Elastic_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Elastic_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Elastic_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Elastic_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Ridge_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Ridge_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Ridge_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Ridge_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Ridge_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Ridge_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Ridge_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Ridge_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Ridge_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Phit_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Phit_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Phit_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phit_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Phit_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Phit_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phit_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Phit_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Phit_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Phis_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Phis_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Phis_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phis_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Phis_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Phis_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Phis_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Phis_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Phis_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_Rhop_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_Rhop_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_Rhop_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Rhop_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_Rhop_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_Rhop_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_Rhop_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_Rhop_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_Rhop_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.interp_IDA_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_IDA_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_IDA_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_IDA_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_IDA_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_IDA_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_IDA_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_IDA_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_IDA_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.interp_Pearson_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Pearson_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_NULL_res_BRCA[[i]]))))
Scan.interp_Pearson_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Pearson_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_TargetScan_res_BRCA[[i]]))))
Scan.interp_Pearson_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Pearson_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Pearson_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Pearson_ENCORI_res_BRCA[[i]]))))

Scan.interp_Spearman_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Spearman_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_NULL_res_BRCA[[i]]))))
Scan.interp_Spearman_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Spearman_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_TargetScan_res_BRCA[[i]]))))
Scan.interp_Spearman_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Spearman_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Spearman_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Spearman_ENCORI_res_BRCA[[i]]))))

Scan.interp_Kendall_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Kendall_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_NULL_res_BRCA[[i]]))))
Scan.interp_Kendall_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Kendall_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_TargetScan_res_BRCA[[i]]))))
Scan.interp_Kendall_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Kendall_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Kendall_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Kendall_ENCORI_res_BRCA[[i]]))))

Scan.interp_Dcor_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Dcor_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_NULL_res_BRCA[[i]]))))
Scan.interp_Dcor_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Dcor_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_TargetScan_res_BRCA[[i]]))))
Scan.interp_Dcor_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dcor_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Dcor_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dcor_ENCORI_res_BRCA[[i]]))))

Scan.interp_RDC_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_RDC_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_RDC_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_NULL_res_BRCA[[i]]))))
Scan.interp_RDC_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_RDC_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_RDC_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_TargetScan_res_BRCA[[i]]))))
Scan.interp_RDC_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_RDC_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_RDC_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_RDC_ENCORI_res_BRCA[[i]]))))

Scan.interp_Hoeffd_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Hoeffd_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_NULL_res_BRCA[[i]]))))
Scan.interp_Hoeffd_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Hoeffd_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_TargetScan_res_BRCA[[i]]))))
Scan.interp_Hoeffd_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Hoeffd_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Hoeffd_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Hoeffd_ENCORI_res_BRCA[[i]]))))

Scan.interp_Zscore_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Zscore_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_NULL_res_BRCA[[i]]))))
Scan.interp_Zscore_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Zscore_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_TargetScan_res_BRCA[[i]]))))
Scan.interp_Zscore_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Zscore_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Zscore_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Zscore_ENCORI_res_BRCA[[i]]))))

Scan.interp_Biweight_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Biweight_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_NULL_res_BRCA[[i]]))))
Scan.interp_Biweight_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Biweight_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_TargetScan_res_BRCA[[i]]))))
Scan.interp_Biweight_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Biweight_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Biweight_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Biweight_ENCORI_res_BRCA[[i]]))))

Scan.interp_Weighted_rank_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Weighted_rank_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_NULL_res_BRCA[[i]]))))
Scan.interp_Weighted_rank_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Weighted_rank_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_TargetScan_res_BRCA[[i]]))))
Scan.interp_Weighted_rank_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Weighted_rank_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Weighted_rank_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Weighted_rank_ENCORI_res_BRCA[[i]]))))

Scan.interp_Cosine_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Cosine_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_NULL_res_BRCA[[i]]))))
Scan.interp_Cosine_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Cosine_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_TargetScan_res_BRCA[[i]]))))
Scan.interp_Cosine_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Cosine_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Cosine_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Cosine_ENCORI_res_BRCA[[i]]))))

Scan.interp_Euclidean_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Euclidean_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_NULL_res_BRCA[[i]]))))
Scan.interp_Euclidean_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Euclidean_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_TargetScan_res_BRCA[[i]]))))
Scan.interp_Euclidean_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Euclidean_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Euclidean_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Euclidean_ENCORI_res_BRCA[[i]]))))

Scan.interp_Manhattan_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Manhattan_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_NULL_res_BRCA[[i]]))))
Scan.interp_Manhattan_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Manhattan_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_TargetScan_res_BRCA[[i]]))))
Scan.interp_Manhattan_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Manhattan_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Manhattan_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Manhattan_ENCORI_res_BRCA[[i]]))))

Scan.interp_Canberra_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Canberra_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_NULL_res_BRCA[[i]]))))
Scan.interp_Canberra_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Canberra_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_TargetScan_res_BRCA[[i]]))))
Scan.interp_Canberra_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Canberra_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Canberra_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Canberra_ENCORI_res_BRCA[[i]]))))

Scan.interp_Chebyshev_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Chebyshev_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_NULL_res_BRCA[[i]]))))
Scan.interp_Chebyshev_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Chebyshev_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_TargetScan_res_BRCA[[i]]))))
Scan.interp_Chebyshev_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Chebyshev_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Chebyshev_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Chebyshev_ENCORI_res_BRCA[[i]]))))

Scan.interp_Dice_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dice_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Dice_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_NULL_res_BRCA[[i]]))))
Scan.interp_Dice_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dice_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Dice_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_TargetScan_res_BRCA[[i]]))))
Scan.interp_Dice_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Dice_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Dice_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Dice_ENCORI_res_BRCA[[i]]))))

Scan.interp_Jaccard_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Jaccard_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_NULL_res_BRCA[[i]]))))
Scan.interp_Jaccard_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Jaccard_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_TargetScan_res_BRCA[[i]]))))
Scan.interp_Jaccard_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Jaccard_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Jaccard_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Jaccard_ENCORI_res_BRCA[[i]]))))

Scan.interp_Mahalanobis_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Mahalanobis_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_NULL_res_BRCA[[i]]))))
Scan.interp_Mahalanobis_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Mahalanobis_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_TargetScan_res_BRCA[[i]]))))
Scan.interp_Mahalanobis_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Mahalanobis_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Mahalanobis_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Mahalanobis_ENCORI_res_BRCA[[i]]))))

Scan.interp_MI_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MI_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_MI_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_NULL_res_BRCA[[i]]))))
Scan.interp_MI_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MI_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_MI_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_TargetScan_res_BRCA[[i]]))))
Scan.interp_MI_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MI_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_MI_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MI_ENCORI_res_BRCA[[i]]))))

Scan.interp_MIC_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MIC_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_MIC_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_NULL_res_BRCA[[i]]))))
Scan.interp_MIC_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MIC_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_MIC_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_TargetScan_res_BRCA[[i]]))))
Scan.interp_MIC_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_MIC_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_MIC_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_MIC_ENCORI_res_BRCA[[i]]))))

Scan.interp_Lasso_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Lasso_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_NULL_res_BRCA[[i]]))))
Scan.interp_Lasso_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Lasso_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_TargetScan_res_BRCA[[i]]))))
Scan.interp_Lasso_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Lasso_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Lasso_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Lasso_ENCORI_res_BRCA[[i]]))))

Scan.interp_Elastic_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Elastic_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_NULL_res_BRCA[[i]]))))
Scan.interp_Elastic_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Elastic_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_TargetScan_res_BRCA[[i]]))))
Scan.interp_Elastic_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Elastic_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Elastic_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Elastic_ENCORI_res_BRCA[[i]]))))

Scan.interp_Ridge_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Ridge_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_NULL_res_BRCA[[i]]))))
Scan.interp_Ridge_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Ridge_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_TargetScan_res_BRCA[[i]]))))
Scan.interp_Ridge_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Ridge_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Ridge_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Ridge_ENCORI_res_BRCA[[i]]))))

Scan.interp_Phit_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phit_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Phit_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_NULL_res_BRCA[[i]]))))
Scan.interp_Phit_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phit_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Phit_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_TargetScan_res_BRCA[[i]]))))
Scan.interp_Phit_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phit_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Phit_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phit_ENCORI_res_BRCA[[i]]))))

Scan.interp_Phis_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phis_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Phis_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_NULL_res_BRCA[[i]]))))
Scan.interp_Phis_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phis_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Phis_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_TargetScan_res_BRCA[[i]]))))
Scan.interp_Phis_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Phis_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Phis_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Phis_ENCORI_res_BRCA[[i]]))))

Scan.interp_Rhop_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_Rhop_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_NULL_res_BRCA[[i]]))))
Scan.interp_Rhop_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_Rhop_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_TargetScan_res_BRCA[[i]]))))
Scan.interp_Rhop_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_Rhop_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_Rhop_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_Rhop_ENCORI_res_BRCA[[i]]))))

Scan.interp_IDA_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_IDA_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_IDA_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_NULL_res_BRCA[[i]]))))
Scan.interp_IDA_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_IDA_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_IDA_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_TargetScan_res_BRCA[[i]]))))
Scan.interp_IDA_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_IDA_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_IDA_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_IDA_ENCORI_res_BRCA[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.interp_Pearson_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Pearson_NULL_res_BRCA)
Scan.interp_Pearson_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Pearson_TargetScan_res_BRCA)
Scan.interp_Pearson_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Pearson_ENCORI_res_BRCA)

Scan.interp_Spearman_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Spearman_NULL_res_BRCA)
Scan.interp_Spearman_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Spearman_TargetScan_res_BRCA)
Scan.interp_Spearman_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Spearman_ENCORI_res_BRCA)

Scan.interp_Kendall_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Kendall_NULL_res_BRCA)
Scan.interp_Kendall_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Kendall_TargetScan_res_BRCA)
Scan.interp_Kendall_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Kendall_ENCORI_res_BRCA)

Scan.interp_Dcor_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dcor_NULL_res_BRCA)
Scan.interp_Dcor_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dcor_TargetScan_res_BRCA)
Scan.interp_Dcor_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dcor_ENCORI_res_BRCA)

Scan.interp_RDC_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_RDC_NULL_res_BRCA)
Scan.interp_RDC_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_RDC_TargetScan_res_BRCA)
Scan.interp_RDC_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_RDC_ENCORI_res_BRCA)

Scan.interp_Hoeffd_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Hoeffd_NULL_res_BRCA)
Scan.interp_Hoeffd_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Hoeffd_TargetScan_res_BRCA)
Scan.interp_Hoeffd_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Hoeffd_ENCORI_res_BRCA)

Scan.interp_Zscore_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Zscore_NULL_res_BRCA)
Scan.interp_Zscore_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Zscore_TargetScan_res_BRCA)
Scan.interp_Zscore_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Zscore_ENCORI_res_BRCA)

Scan.interp_Biweight_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Biweight_NULL_res_BRCA)
Scan.interp_Biweight_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Biweight_TargetScan_res_BRCA)
Scan.interp_Biweight_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Biweight_ENCORI_res_BRCA)

Scan.interp_Weighted_rank_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Weighted_rank_NULL_res_BRCA)
Scan.interp_Weighted_rank_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Weighted_rank_TargetScan_res_BRCA)
Scan.interp_Weighted_rank_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Weighted_rank_ENCORI_res_BRCA)

Scan.interp_Cosine_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Cosine_NULL_res_BRCA)
Scan.interp_Cosine_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Cosine_TargetScan_res_BRCA)
Scan.interp_Cosine_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Cosine_ENCORI_res_BRCA)

Scan.interp_Euclidean_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Euclidean_NULL_res_BRCA)
Scan.interp_Euclidean_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Euclidean_TargetScan_res_BRCA)
Scan.interp_Euclidean_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Euclidean_ENCORI_res_BRCA)

Scan.interp_Manhattan_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Manhattan_NULL_res_BRCA)
Scan.interp_Manhattan_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Manhattan_TargetScan_res_BRCA)
Scan.interp_Manhattan_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Manhattan_ENCORI_res_BRCA)

Scan.interp_Canberra_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Canberra_NULL_res_BRCA)
Scan.interp_Canberra_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Canberra_TargetScan_res_BRCA)
Scan.interp_Canberra_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Canberra_ENCORI_res_BRCA)

Scan.interp_Chebyshev_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Chebyshev_NULL_res_BRCA)
Scan.interp_Chebyshev_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Chebyshev_TargetScan_res_BRCA)
Scan.interp_Chebyshev_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Chebyshev_ENCORI_res_BRCA)

Scan.interp_Dice_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dice_NULL_res_BRCA)
Scan.interp_Dice_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dice_TargetScan_res_BRCA)
Scan.interp_Dice_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Dice_ENCORI_res_BRCA)

Scan.interp_Jaccard_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Jaccard_NULL_res_BRCA)
Scan.interp_Jaccard_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Jaccard_TargetScan_res_BRCA)
Scan.interp_Jaccard_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Jaccard_ENCORI_res_BRCA)

Scan.interp_Mahalanobis_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Mahalanobis_NULL_res_BRCA)
Scan.interp_Mahalanobis_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Mahalanobis_TargetScan_res_BRCA)
Scan.interp_Mahalanobis_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Mahalanobis_ENCORI_res_BRCA)

Scan.interp_MI_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MI_NULL_res_BRCA)
Scan.interp_MI_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MI_TargetScan_res_BRCA)
Scan.interp_MI_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MI_ENCORI_res_BRCA)

Scan.interp_MIC_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MIC_NULL_res_BRCA)
Scan.interp_MIC_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MIC_TargetScan_res_BRCA)
Scan.interp_MIC_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_MIC_ENCORI_res_BRCA)

Scan.interp_Lasso_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Lasso_NULL_res_BRCA)
Scan.interp_Lasso_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Lasso_TargetScan_res_BRCA)
Scan.interp_Lasso_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Lasso_ENCORI_res_BRCA)

Scan.interp_Elastic_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Elastic_NULL_res_BRCA)
Scan.interp_Elastic_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Elastic_TargetScan_res_BRCA)
Scan.interp_Elastic_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Elastic_ENCORI_res_BRCA)

Scan.interp_Ridge_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Ridge_NULL_res_BRCA)
Scan.interp_Ridge_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Ridge_TargetScan_res_BRCA)
Scan.interp_Ridge_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Ridge_ENCORI_res_BRCA)

Scan.interp_Phit_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phit_NULL_res_BRCA)
Scan.interp_Phit_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phit_TargetScan_res_BRCA)
Scan.interp_Phit_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phit_ENCORI_res_BRCA)

Scan.interp_Phis_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phis_NULL_res_BRCA)
Scan.interp_Phis_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phis_TargetScan_res_BRCA)
Scan.interp_Phis_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Phis_ENCORI_res_BRCA)

Scan.interp_Rhop_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Rhop_NULL_res_BRCA)
Scan.interp_Rhop_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Rhop_TargetScan_res_BRCA)
Scan.interp_Rhop_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_Rhop_ENCORI_res_BRCA)

Scan.interp_IDA_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_IDA_NULL_res_BRCA)
Scan.interp_IDA_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_IDA_TargetScan_res_BRCA)
Scan.interp_IDA_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_IDA_ENCORI_res_BRCA)

save.image("Scan.interp_BRCA.RData")


#######################################################################################################################################################################################
########################################################################### Scan.perturb application in K562 dataset ##################################################################
#######################################################################################################################################################################################

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

set.seed(123)
ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

# No prior information
Scan.perturb_Pearson_timestart <- Sys.time()
Scan.perturb_Pearson_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Pearson")
Scan.perturb_Pearson_timeend <- Sys.time()
Scan.perturb_Pearson_runningtime_NULL <- Scan.perturb_Pearson_timeend - Scan.perturb_Pearson_timestart

Scan.perturb_Spearman_timestart <- Sys.time()
Scan.perturb_Spearman_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Spearman")
Scan.perturb_Spearman_timeend <- Sys.time()
Scan.perturb_Spearman_runningtime_NULL <- Scan.perturb_Spearman_timeend - Scan.perturb_Spearman_timestart

Scan.perturb_Kendall_timestart <- Sys.time()
Scan.perturb_Kendall_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Kendall")
Scan.perturb_Kendall_timeend <- Sys.time()
Scan.perturb_Kendall_runningtime_NULL <- Scan.perturb_Kendall_timeend - Scan.perturb_Kendall_timestart

Scan.perturb_Dcor_timestart <- Sys.time()
Scan.perturb_Dcor_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Dcor")
Scan.perturb_Dcor_timeend <- Sys.time()
Scan.perturb_Dcor_runningtime_NULL <- Scan.perturb_Dcor_timeend - Scan.perturb_Dcor_timestart

Scan.perturb_RDC_timestart <- Sys.time()
Scan.perturb_RDC_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "RDC")
Scan.perturb_RDC_timeend <- Sys.time()
Scan.perturb_RDC_runningtime_NULL <- Scan.perturb_RDC_timeend - Scan.perturb_RDC_timestart

Scan.perturb_Hoeffd_timestart <- Sys.time()
Scan.perturb_Hoeffd_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Hoeffd")
Scan.perturb_Hoeffd_timeend <- Sys.time()
Scan.perturb_Hoeffd_runningtime_NULL <- Scan.perturb_Hoeffd_timeend - Scan.perturb_Hoeffd_timestart

Scan.perturb_Zscore_timestart <- Sys.time()
Scan.perturb_Zscore_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Zscore")
Scan.perturb_Zscore_timeend <- Sys.time()
Scan.perturb_Zscore_runningtime_NULL <- Scan.perturb_Zscore_timeend - Scan.perturb_Zscore_timestart

Scan.perturb_Biweight_timestart <- Sys.time()
Scan.perturb_Biweight_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Biweight")
Scan.perturb_Biweight_timeend <- Sys.time()
Scan.perturb_Biweight_runningtime_NULL <- Scan.perturb_Biweight_timeend - Scan.perturb_Biweight_timestart

Scan.perturb_Weighted_rank_timestart <- Sys.time()
Scan.perturb_Weighted_rank_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Weighted_rank")
Scan.perturb_Weighted_rank_timeend <- Sys.time()
Scan.perturb_Weighted_rank_runningtime_NULL <- Scan.perturb_Weighted_rank_timeend - Scan.perturb_Weighted_rank_timestart

Scan.perturb_Cosine_timestart <- Sys.time()
Scan.perturb_Cosine_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Cosine")
Scan.perturb_Cosine_timeend <- Sys.time()
Scan.perturb_Cosine_runningtime_NULL <- Scan.perturb_Cosine_timeend - Scan.perturb_Cosine_timestart

Scan.perturb_Euclidean_timestart <- Sys.time()
Scan.perturb_Euclidean_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Euclidean")
Scan.perturb_Euclidean_timeend <- Sys.time()
Scan.perturb_Euclidean_runningtime_NULL <- Scan.perturb_Euclidean_timeend - Scan.perturb_Euclidean_timestart

Scan.perturb_Manhattan_timestart <- Sys.time()
Scan.perturb_Manhattan_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Manhattan")
Scan.perturb_Manhattan_timeend <- Sys.time()
Scan.perturb_Manhattan_runningtime_NULL <- Scan.perturb_Manhattan_timeend - Scan.perturb_Manhattan_timestart

Scan.perturb_Canberra_timestart <- Sys.time()
Scan.perturb_Canberra_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Canberra")
Scan.perturb_Canberra_timeend <- Sys.time()
Scan.perturb_Canberra_runningtime_NULL <- Scan.perturb_Canberra_timeend - Scan.perturb_Canberra_timestart

Scan.perturb_Chebyshev_timestart <- Sys.time()
Scan.perturb_Chebyshev_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Chebyshev")
Scan.perturb_Chebyshev_timeend <- Sys.time()
Scan.perturb_Chebyshev_runningtime_NULL <- Scan.perturb_Chebyshev_timeend - Scan.perturb_Chebyshev_timestart

Scan.perturb_Dice_timestart <- Sys.time()
Scan.perturb_Dice_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Dice")
Scan.perturb_Dice_timeend <- Sys.time()
Scan.perturb_Dice_runningtime_NULL <- Scan.perturb_Dice_timeend - Scan.perturb_Dice_timestart

Scan.perturb_Jaccard_timestart <- Sys.time()
Scan.perturb_Jaccard_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Jaccard")
Scan.perturb_Jaccard_timeend <- Sys.time()
Scan.perturb_Jaccard_runningtime_NULL <- Scan.perturb_Jaccard_timeend - Scan.perturb_Jaccard_timestart

Scan.perturb_Mahalanobis_timestart <- Sys.time()
Scan.perturb_Mahalanobis_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Mahalanobis")
Scan.perturb_Mahalanobis_timeend <- Sys.time()
Scan.perturb_Mahalanobis_runningtime_NULL <- Scan.perturb_Mahalanobis_timeend - Scan.perturb_Mahalanobis_timestart

Scan.perturb_MI_timestart <- Sys.time()
Scan.perturb_MI_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "MI")
Scan.perturb_MI_timeend <- Sys.time()
Scan.perturb_MI_runningtime_NULL <- Scan.perturb_MI_timeend - Scan.perturb_MI_timestart

Scan.perturb_MIC_timestart <- Sys.time()
Scan.perturb_MIC_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "MIC")
Scan.perturb_MIC_timeend <- Sys.time()
Scan.perturb_MIC_runningtime_NULL <- Scan.perturb_MIC_timeend - Scan.perturb_MIC_timestart

Scan.perturb_Lasso_timestart <- Sys.time()
Scan.perturb_Lasso_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Lasso")
Scan.perturb_Lasso_timeend <- Sys.time()
Scan.perturb_Lasso_runningtime_NULL <- Scan.perturb_Lasso_timeend - Scan.perturb_Lasso_timestart

Scan.perturb_Elastic_timestart <- Sys.time()
Scan.perturb_Elastic_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Elastic")
Scan.perturb_Elastic_timeend <- Sys.time()
Scan.perturb_Elastic_runningtime_NULL <- Scan.perturb_Elastic_timeend - Scan.perturb_Elastic_timestart

Scan.perturb_Ridge_timestart <- Sys.time()
Scan.perturb_Ridge_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Ridge")
Scan.perturb_Ridge_timeend <- Sys.time()
Scan.perturb_Ridge_runningtime_NULL <- Scan.perturb_Ridge_timeend - Scan.perturb_Ridge_timestart

Scan.perturb_Phit_timestart <- Sys.time()
Scan.perturb_Phit_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Phit")
Scan.perturb_Phit_timeend <- Sys.time()
Scan.perturb_Phit_runningtime_NULL <- Scan.perturb_Phit_timeend - Scan.perturb_Phit_timestart

Scan.perturb_Phis_timestart <- Sys.time()
Scan.perturb_Phis_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Phis")
Scan.perturb_Phis_timeend <- Sys.time()
Scan.perturb_Phis_runningtime_NULL <- Scan.perturb_Phis_timeend - Scan.perturb_Phis_timestart

Scan.perturb_Rhop_timestart <- Sys.time()
Scan.perturb_Rhop_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "Rhop")
Scan.perturb_Rhop_timeend <- Sys.time()
Scan.perturb_Rhop_runningtime_NULL <- Scan.perturb_Rhop_timeend - Scan.perturb_Rhop_timestart

Scan.perturb_IDA_timestart <- Sys.time()
Scan.perturb_IDA_NULL_res <- Scan.perturb(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, method = "IDA", pcmethod = "stable")
Scan.perturb_IDA_timeend <- Sys.time()
Scan.perturb_IDA_runningtime_NULL <- Scan.perturb_IDA_timeend - Scan.perturb_IDA_timestart

# The prior information is TargetScan
Scan.perturb_Pearson_TargetScan_res <- lapply(seq(Scan.perturb_Pearson_NULL_res), function(i) Scan.perturb_Pearson_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Spearman_TargetScan_res <- lapply(seq(Scan.perturb_Spearman_NULL_res), function(i) Scan.perturb_Spearman_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Kendall_TargetScan_res <- lapply(seq(Scan.perturb_Kendall_NULL_res), function(i) Scan.perturb_Kendall_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Dcor_TargetScan_res <- lapply(seq(Scan.perturb_Dcor_NULL_res), function(i) Scan.perturb_Dcor_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_RDC_TargetScan_res <- lapply(seq(Scan.perturb_RDC_NULL_res), function(i) Scan.perturb_RDC_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Hoeffd_TargetScan_res <- lapply(seq(Scan.perturb_Hoeffd_NULL_res), function(i) Scan.perturb_Hoeffd_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Zscore_TargetScan_res <- lapply(seq(Scan.perturb_Zscore_NULL_res), function(i) Scan.perturb_Zscore_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Biweight_TargetScan_res <- lapply(seq(Scan.perturb_Biweight_NULL_res), function(i) Scan.perturb_Biweight_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Weighted_rank_TargetScan_res <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res), function(i) Scan.perturb_Weighted_rank_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Cosine_TargetScan_res <- lapply(seq(Scan.perturb_Cosine_NULL_res), function(i) Scan.perturb_Cosine_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Euclidean_TargetScan_res <- lapply(seq(Scan.perturb_Euclidean_NULL_res), function(i) Scan.perturb_Euclidean_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Manhattan_TargetScan_res <- lapply(seq(Scan.perturb_Manhattan_NULL_res), function(i) Scan.perturb_Manhattan_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Canberra_TargetScan_res <- lapply(seq(Scan.perturb_Canberra_NULL_res), function(i) Scan.perturb_Canberra_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Chebyshev_TargetScan_res <- lapply(seq(Scan.perturb_Chebyshev_NULL_res), function(i) Scan.perturb_Chebyshev_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Dice_TargetScan_res <- lapply(seq(Scan.perturb_Dice_NULL_res), function(i) Scan.perturb_Dice_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Jaccard_TargetScan_res <- lapply(seq(Scan.perturb_Jaccard_NULL_res), function(i) Scan.perturb_Jaccard_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Mahalanobis_TargetScan_res <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res), function(i) Scan.perturb_Mahalanobis_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_MI_TargetScan_res <- lapply(seq(Scan.perturb_MI_NULL_res), function(i) Scan.perturb_MI_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_MIC_TargetScan_res <- lapply(seq(Scan.perturb_MIC_NULL_res), function(i) Scan.perturb_MIC_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Lasso_TargetScan_res <- lapply(seq(Scan.perturb_Lasso_NULL_res), function(i) Scan.perturb_Lasso_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Elastic_TargetScan_res <- lapply(seq(Scan.perturb_Elastic_NULL_res), function(i) Scan.perturb_Elastic_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Ridge_TargetScan_res <- lapply(seq(Scan.perturb_Ridge_NULL_res), function(i) Scan.perturb_Ridge_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Phit_TargetScan_res <- lapply(seq(Scan.perturb_Phit_NULL_res), function(i) Scan.perturb_Phit_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Phis_TargetScan_res <- lapply(seq(Scan.perturb_Phis_NULL_res), function(i) Scan.perturb_Phis_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_Rhop_TargetScan_res <- lapply(seq(Scan.perturb_Rhop_NULL_res), function(i) Scan.perturb_Rhop_NULL_res[[i]] %s% TargetScan_graph)
Scan.perturb_IDA_TargetScan_res <- lapply(seq(Scan.perturb_IDA_NULL_res), function(i) Scan.perturb_IDA_NULL_res[[i]] %s% TargetScan_graph)

# The prior information is ENCORI
Scan.perturb_Pearson_ENCORI_res <- lapply(seq(Scan.perturb_Pearson_NULL_res), function(i) Scan.perturb_Pearson_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Spearman_ENCORI_res <- lapply(seq(Scan.perturb_Spearman_NULL_res), function(i) Scan.perturb_Spearman_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Kendall_ENCORI_res <- lapply(seq(Scan.perturb_Kendall_NULL_res), function(i) Scan.perturb_Kendall_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Dcor_ENCORI_res <- lapply(seq(Scan.perturb_Dcor_NULL_res), function(i) Scan.perturb_Dcor_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_RDC_ENCORI_res <- lapply(seq(Scan.perturb_RDC_NULL_res), function(i) Scan.perturb_RDC_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Hoeffd_ENCORI_res <- lapply(seq(Scan.perturb_Hoeffd_NULL_res), function(i) Scan.perturb_Hoeffd_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Zscore_ENCORI_res <- lapply(seq(Scan.perturb_Zscore_NULL_res), function(i) Scan.perturb_Zscore_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Biweight_ENCORI_res <- lapply(seq(Scan.perturb_Biweight_NULL_res), function(i) Scan.perturb_Biweight_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Weighted_rank_ENCORI_res <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res), function(i) Scan.perturb_Weighted_rank_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Cosine_ENCORI_res <- lapply(seq(Scan.perturb_Cosine_NULL_res), function(i) Scan.perturb_Cosine_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Euclidean_ENCORI_res <- lapply(seq(Scan.perturb_Euclidean_NULL_res), function(i) Scan.perturb_Euclidean_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Manhattan_ENCORI_res <- lapply(seq(Scan.perturb_Manhattan_NULL_res), function(i) Scan.perturb_Manhattan_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Canberra_ENCORI_res <- lapply(seq(Scan.perturb_Canberra_NULL_res), function(i) Scan.perturb_Canberra_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Chebyshev_ENCORI_res <- lapply(seq(Scan.perturb_Chebyshev_NULL_res), function(i) Scan.perturb_Chebyshev_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Dice_ENCORI_res <- lapply(seq(Scan.perturb_Dice_NULL_res), function(i) Scan.perturb_Dice_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Jaccard_ENCORI_res <- lapply(seq(Scan.perturb_Jaccard_NULL_res), function(i) Scan.perturb_Jaccard_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Mahalanobis_ENCORI_res <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res), function(i) Scan.perturb_Mahalanobis_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_MI_ENCORI_res <- lapply(seq(Scan.perturb_MI_NULL_res), function(i) Scan.perturb_MI_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_MIC_ENCORI_res <- lapply(seq(Scan.perturb_MIC_NULL_res), function(i) Scan.perturb_MIC_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Lasso_ENCORI_res <- lapply(seq(Scan.perturb_Lasso_NULL_res), function(i) Scan.perturb_Lasso_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Elastic_ENCORI_res <- lapply(seq(Scan.perturb_Elastic_NULL_res), function(i) Scan.perturb_Elastic_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Ridge_ENCORI_res <- lapply(seq(Scan.perturb_Ridge_NULL_res), function(i) Scan.perturb_Ridge_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Phit_ENCORI_res <- lapply(seq(Scan.perturb_Phit_NULL_res), function(i) Scan.perturb_Phit_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Phis_ENCORI_res <- lapply(seq(Scan.perturb_Phis_NULL_res), function(i) Scan.perturb_Phis_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_Rhop_ENCORI_res <- lapply(seq(Scan.perturb_Rhop_NULL_res), function(i) Scan.perturb_Rhop_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_IDA_ENCORI_res <- lapply(seq(Scan.perturb_IDA_NULL_res), function(i) Scan.perturb_IDA_NULL_res[[i]] %s% ENCORI_graph)


# Number of predicted sample-specific interactions
Scan.perturb_Pearson_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Pearson_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Pearson_NULL_res[[i]] ))))
Scan.perturb_Pearson_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Pearson_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Pearson_TargetScan_res[[i]] ))))
Scan.perturb_Pearson_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Pearson_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Pearson_ENCORI_res[[i]] ))))

Scan.perturb_Spearman_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Spearman_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Spearman_NULL_res[[i]] ))))
Scan.perturb_Spearman_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Spearman_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Spearman_TargetScan_res[[i]] ))))
Scan.perturb_Spearman_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Spearman_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Spearman_ENCORI_res[[i]] ))))

Scan.perturb_Kendall_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Kendall_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Kendall_NULL_res[[i]] ))))
Scan.perturb_Kendall_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Kendall_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Kendall_TargetScan_res[[i]] ))))
Scan.perturb_Kendall_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Kendall_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Kendall_ENCORI_res[[i]] ))))

Scan.perturb_Dcor_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Dcor_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Dcor_NULL_res[[i]] ))))
Scan.perturb_Dcor_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Dcor_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Dcor_TargetScan_res[[i]] ))))
Scan.perturb_Dcor_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Dcor_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Dcor_ENCORI_res[[i]] ))))

Scan.perturb_RDC_NULL_res_num <- unlist(lapply(seq(Scan.perturb_RDC_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_RDC_NULL_res[[i]] ))))
Scan.perturb_RDC_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_RDC_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_RDC_TargetScan_res[[i]] ))))
Scan.perturb_RDC_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_RDC_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_RDC_ENCORI_res[[i]] ))))

Scan.perturb_Hoeffd_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_NULL_res[[i]] ))))
Scan.perturb_Hoeffd_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_TargetScan_res[[i]] ))))
Scan.perturb_Hoeffd_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_ENCORI_res[[i]] ))))

Scan.perturb_Zscore_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Zscore_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Zscore_NULL_res[[i]] ))))
Scan.perturb_Zscore_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Zscore_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Zscore_TargetScan_res[[i]] ))))
Scan.perturb_Zscore_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Zscore_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Zscore_ENCORI_res[[i]] ))))

Scan.perturb_Biweight_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Biweight_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Biweight_NULL_res[[i]] ))))
Scan.perturb_Biweight_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Biweight_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Biweight_TargetScan_res[[i]] ))))
Scan.perturb_Biweight_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Biweight_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Biweight_ENCORI_res[[i]] ))))

Scan.perturb_Weighted_rank_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_NULL_res[[i]] ))))
Scan.perturb_Weighted_rank_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res[[i]] ))))
Scan.perturb_Weighted_rank_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res[[i]] ))))

Scan.perturb_Cosine_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Cosine_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Cosine_NULL_res[[i]] ))))
Scan.perturb_Cosine_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Cosine_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Cosine_TargetScan_res[[i]] ))))
Scan.perturb_Cosine_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Cosine_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Cosine_ENCORI_res[[i]] ))))

Scan.perturb_Euclidean_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Euclidean_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_NULL_res[[i]] ))))
Scan.perturb_Euclidean_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Euclidean_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_TargetScan_res[[i]] ))))
Scan.perturb_Euclidean_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Euclidean_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_ENCORI_res[[i]] ))))

Scan.perturb_Manhattan_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Manhattan_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_NULL_res[[i]] ))))
Scan.perturb_Manhattan_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Manhattan_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_TargetScan_res[[i]] ))))
Scan.perturb_Manhattan_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Manhattan_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_ENCORI_res[[i]] ))))

Scan.perturb_Canberra_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Canberra_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Canberra_NULL_res[[i]] ))))
Scan.perturb_Canberra_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Canberra_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Canberra_TargetScan_res[[i]] ))))
Scan.perturb_Canberra_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Canberra_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Canberra_ENCORI_res[[i]] ))))

Scan.perturb_Chebyshev_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_NULL_res[[i]] ))))
Scan.perturb_Chebyshev_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_TargetScan_res[[i]] ))))
Scan.perturb_Chebyshev_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_ENCORI_res[[i]] ))))

Scan.perturb_Dice_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Dice_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Dice_NULL_res[[i]] ))))
Scan.perturb_Dice_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Dice_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Dice_TargetScan_res[[i]] ))))
Scan.perturb_Dice_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Dice_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Dice_ENCORI_res[[i]] ))))

Scan.perturb_Jaccard_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Jaccard_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_NULL_res[[i]] ))))
Scan.perturb_Jaccard_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Jaccard_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_TargetScan_res[[i]] ))))
Scan.perturb_Jaccard_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Jaccard_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_ENCORI_res[[i]] ))))

Scan.perturb_Mahalanobis_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_NULL_res[[i]] ))))
Scan.perturb_Mahalanobis_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res[[i]] ))))
Scan.perturb_Mahalanobis_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res[[i]] ))))

Scan.perturb_MI_NULL_res_num <- unlist(lapply(seq(Scan.perturb_MI_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_MI_NULL_res[[i]] ))))
Scan.perturb_MI_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_MI_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_MI_TargetScan_res[[i]] ))))
Scan.perturb_MI_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_MI_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_MI_ENCORI_res[[i]] ))))

Scan.perturb_MIC_NULL_res_num <- unlist(lapply(seq(Scan.perturb_MIC_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_MIC_NULL_res[[i]] ))))
Scan.perturb_MIC_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_MIC_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_MIC_TargetScan_res[[i]] ))))
Scan.perturb_MIC_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_MIC_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_MIC_ENCORI_res[[i]] ))))

Scan.perturb_Lasso_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Lasso_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Lasso_NULL_res[[i]] ))))
Scan.perturb_Lasso_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Lasso_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Lasso_TargetScan_res[[i]] ))))
Scan.perturb_Lasso_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Lasso_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Lasso_ENCORI_res[[i]] ))))

Scan.perturb_Elastic_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Elastic_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Elastic_NULL_res[[i]] ))))
Scan.perturb_Elastic_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Elastic_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Elastic_TargetScan_res[[i]] ))))
Scan.perturb_Elastic_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Elastic_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Elastic_ENCORI_res[[i]] ))))

Scan.perturb_Ridge_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Ridge_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Ridge_NULL_res[[i]] ))))
Scan.perturb_Ridge_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Ridge_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Ridge_TargetScan_res[[i]] ))))
Scan.perturb_Ridge_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Ridge_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Ridge_ENCORI_res[[i]] ))))

Scan.perturb_Phit_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Phit_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Phit_NULL_res[[i]] ))))
Scan.perturb_Phit_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Phit_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Phit_TargetScan_res[[i]] ))))
Scan.perturb_Phit_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Phit_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Phit_ENCORI_res[[i]] ))))

Scan.perturb_Phis_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Phis_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Phis_NULL_res[[i]] ))))
Scan.perturb_Phis_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Phis_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Phis_TargetScan_res[[i]] ))))
Scan.perturb_Phis_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Phis_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Phis_ENCORI_res[[i]] ))))

Scan.perturb_Rhop_NULL_res_num <- unlist(lapply(seq(Scan.perturb_Rhop_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_Rhop_NULL_res[[i]] ))))
Scan.perturb_Rhop_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_Rhop_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_Rhop_TargetScan_res[[i]] ))))
Scan.perturb_Rhop_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_Rhop_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_Rhop_ENCORI_res[[i]] ))))

Scan.perturb_IDA_NULL_res_num <- unlist(lapply(seq(Scan.perturb_IDA_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_IDA_NULL_res[[i]] ))))
Scan.perturb_IDA_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_IDA_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_IDA_TargetScan_res[[i]] ))))
Scan.perturb_IDA_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_IDA_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_IDA_ENCORI_res[[i]] ))))


# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.perturb_Pearson_NULL_res_validated <- lapply(seq(Scan.perturb_Pearson_NULL_res), function(i) as_data_frame(Scan.perturb_Pearson_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Pearson_TargetScan_res_validated <- lapply(seq(Scan.perturb_Pearson_TargetScan_res), function(i) as_data_frame(Scan.perturb_Pearson_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Pearson_ENCORI_res_validated <- lapply(seq(Scan.perturb_Pearson_ENCORI_res), function(i) as_data_frame(Scan.perturb_Pearson_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Spearman_NULL_res_validated <- lapply(seq(Scan.perturb_Spearman_NULL_res), function(i) as_data_frame(Scan.perturb_Spearman_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Spearman_TargetScan_res_validated <- lapply(seq(Scan.perturb_Spearman_TargetScan_res), function(i) as_data_frame(Scan.perturb_Spearman_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Spearman_ENCORI_res_validated <- lapply(seq(Scan.perturb_Spearman_ENCORI_res), function(i) as_data_frame(Scan.perturb_Spearman_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Kendall_NULL_res_validated <- lapply(seq(Scan.perturb_Kendall_NULL_res), function(i) as_data_frame(Scan.perturb_Kendall_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Kendall_TargetScan_res_validated <- lapply(seq(Scan.perturb_Kendall_TargetScan_res), function(i) as_data_frame(Scan.perturb_Kendall_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Kendall_ENCORI_res_validated <- lapply(seq(Scan.perturb_Kendall_ENCORI_res), function(i) as_data_frame(Scan.perturb_Kendall_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Dcor_NULL_res_validated <- lapply(seq(Scan.perturb_Dcor_NULL_res), function(i) as_data_frame(Scan.perturb_Dcor_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dcor_TargetScan_res_validated <- lapply(seq(Scan.perturb_Dcor_TargetScan_res), function(i) as_data_frame(Scan.perturb_Dcor_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dcor_ENCORI_res_validated <- lapply(seq(Scan.perturb_Dcor_ENCORI_res), function(i) as_data_frame(Scan.perturb_Dcor_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_RDC_NULL_res_validated <- lapply(seq(Scan.perturb_RDC_NULL_res), function(i) as_data_frame(Scan.perturb_RDC_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_RDC_TargetScan_res_validated <- lapply(seq(Scan.perturb_RDC_TargetScan_res), function(i) as_data_frame(Scan.perturb_RDC_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_RDC_ENCORI_res_validated <- lapply(seq(Scan.perturb_RDC_ENCORI_res), function(i) as_data_frame(Scan.perturb_RDC_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Hoeffd_NULL_res_validated <- lapply(seq(Scan.perturb_Hoeffd_NULL_res), function(i) as_data_frame(Scan.perturb_Hoeffd_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Hoeffd_TargetScan_res_validated <- lapply(seq(Scan.perturb_Hoeffd_TargetScan_res), function(i) as_data_frame(Scan.perturb_Hoeffd_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Hoeffd_ENCORI_res_validated <- lapply(seq(Scan.perturb_Hoeffd_ENCORI_res), function(i) as_data_frame(Scan.perturb_Hoeffd_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Zscore_NULL_res_validated <- lapply(seq(Scan.perturb_Zscore_NULL_res), function(i) as_data_frame(Scan.perturb_Zscore_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Zscore_TargetScan_res_validated <- lapply(seq(Scan.perturb_Zscore_TargetScan_res), function(i) as_data_frame(Scan.perturb_Zscore_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Zscore_ENCORI_res_validated <- lapply(seq(Scan.perturb_Zscore_ENCORI_res), function(i) as_data_frame(Scan.perturb_Zscore_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Biweight_NULL_res_validated <- lapply(seq(Scan.perturb_Biweight_NULL_res), function(i) as_data_frame(Scan.perturb_Biweight_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Biweight_TargetScan_res_validated <- lapply(seq(Scan.perturb_Biweight_TargetScan_res), function(i) as_data_frame(Scan.perturb_Biweight_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Biweight_ENCORI_res_validated <- lapply(seq(Scan.perturb_Biweight_ENCORI_res), function(i) as_data_frame(Scan.perturb_Biweight_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Weighted_rank_NULL_res_validated <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res), function(i) as_data_frame(Scan.perturb_Weighted_rank_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Weighted_rank_TargetScan_res_validated <- lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res), function(i) as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Weighted_rank_ENCORI_res_validated <- lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res), function(i) as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Cosine_NULL_res_validated <- lapply(seq(Scan.perturb_Cosine_NULL_res), function(i) as_data_frame(Scan.perturb_Cosine_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Cosine_TargetScan_res_validated <- lapply(seq(Scan.perturb_Cosine_TargetScan_res), function(i) as_data_frame(Scan.perturb_Cosine_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Cosine_ENCORI_res_validated <- lapply(seq(Scan.perturb_Cosine_ENCORI_res), function(i) as_data_frame(Scan.perturb_Cosine_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Euclidean_NULL_res_validated <- lapply(seq(Scan.perturb_Euclidean_NULL_res), function(i) as_data_frame(Scan.perturb_Euclidean_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Euclidean_TargetScan_res_validated <- lapply(seq(Scan.perturb_Euclidean_TargetScan_res), function(i) as_data_frame(Scan.perturb_Euclidean_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Euclidean_ENCORI_res_validated <- lapply(seq(Scan.perturb_Euclidean_ENCORI_res), function(i) as_data_frame(Scan.perturb_Euclidean_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Manhattan_NULL_res_validated <- lapply(seq(Scan.perturb_Manhattan_NULL_res), function(i) as_data_frame(Scan.perturb_Manhattan_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Manhattan_TargetScan_res_validated <- lapply(seq(Scan.perturb_Manhattan_TargetScan_res), function(i) as_data_frame(Scan.perturb_Manhattan_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Manhattan_ENCORI_res_validated <- lapply(seq(Scan.perturb_Manhattan_ENCORI_res), function(i) as_data_frame(Scan.perturb_Manhattan_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Canberra_NULL_res_validated <- lapply(seq(Scan.perturb_Canberra_NULL_res), function(i) as_data_frame(Scan.perturb_Canberra_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Canberra_TargetScan_res_validated <- lapply(seq(Scan.perturb_Canberra_TargetScan_res), function(i) as_data_frame(Scan.perturb_Canberra_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Canberra_ENCORI_res_validated <- lapply(seq(Scan.perturb_Canberra_ENCORI_res), function(i) as_data_frame(Scan.perturb_Canberra_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Chebyshev_NULL_res_validated <- lapply(seq(Scan.perturb_Chebyshev_NULL_res), function(i) as_data_frame(Scan.perturb_Chebyshev_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Chebyshev_TargetScan_res_validated <- lapply(seq(Scan.perturb_Chebyshev_TargetScan_res), function(i) as_data_frame(Scan.perturb_Chebyshev_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Chebyshev_ENCORI_res_validated <- lapply(seq(Scan.perturb_Chebyshev_ENCORI_res), function(i) as_data_frame(Scan.perturb_Chebyshev_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Dice_NULL_res_validated <- lapply(seq(Scan.perturb_Dice_NULL_res), function(i) as_data_frame(Scan.perturb_Dice_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dice_TargetScan_res_validated <- lapply(seq(Scan.perturb_Dice_TargetScan_res), function(i) as_data_frame(Scan.perturb_Dice_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dice_ENCORI_res_validated <- lapply(seq(Scan.perturb_Dice_ENCORI_res), function(i) as_data_frame(Scan.perturb_Dice_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Jaccard_NULL_res_validated <- lapply(seq(Scan.perturb_Jaccard_NULL_res), function(i) as_data_frame(Scan.perturb_Jaccard_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Jaccard_TargetScan_res_validated <- lapply(seq(Scan.perturb_Jaccard_TargetScan_res), function(i) as_data_frame(Scan.perturb_Jaccard_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Jaccard_ENCORI_res_validated <- lapply(seq(Scan.perturb_Jaccard_ENCORI_res), function(i) as_data_frame(Scan.perturb_Jaccard_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Mahalanobis_NULL_res_validated <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res), function(i) as_data_frame(Scan.perturb_Mahalanobis_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Mahalanobis_TargetScan_res_validated <- lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res), function(i) as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Mahalanobis_ENCORI_res_validated <- lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res), function(i) as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_MI_NULL_res_validated <- lapply(seq(Scan.perturb_MI_NULL_res), function(i) as_data_frame(Scan.perturb_MI_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MI_TargetScan_res_validated <- lapply(seq(Scan.perturb_MI_TargetScan_res), function(i) as_data_frame(Scan.perturb_MI_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MI_ENCORI_res_validated <- lapply(seq(Scan.perturb_MI_ENCORI_res), function(i) as_data_frame(Scan.perturb_MI_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_MIC_NULL_res_validated <- lapply(seq(Scan.perturb_MIC_NULL_res), function(i) as_data_frame(Scan.perturb_MIC_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MIC_TargetScan_res_validated <- lapply(seq(Scan.perturb_MIC_TargetScan_res), function(i) as_data_frame(Scan.perturb_MIC_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MIC_ENCORI_res_validated <- lapply(seq(Scan.perturb_MIC_ENCORI_res), function(i) as_data_frame(Scan.perturb_MIC_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Lasso_NULL_res_validated <- lapply(seq(Scan.perturb_Lasso_NULL_res), function(i) as_data_frame(Scan.perturb_Lasso_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Lasso_TargetScan_res_validated <- lapply(seq(Scan.perturb_Lasso_TargetScan_res), function(i) as_data_frame(Scan.perturb_Lasso_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Lasso_ENCORI_res_validated <- lapply(seq(Scan.perturb_Lasso_ENCORI_res), function(i) as_data_frame(Scan.perturb_Lasso_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Elastic_NULL_res_validated <- lapply(seq(Scan.perturb_Elastic_NULL_res), function(i) as_data_frame(Scan.perturb_Elastic_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Elastic_TargetScan_res_validated <- lapply(seq(Scan.perturb_Elastic_TargetScan_res), function(i) as_data_frame(Scan.perturb_Elastic_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Elastic_ENCORI_res_validated <- lapply(seq(Scan.perturb_Elastic_ENCORI_res), function(i) as_data_frame(Scan.perturb_Elastic_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Ridge_NULL_res_validated <- lapply(seq(Scan.perturb_Ridge_NULL_res), function(i) as_data_frame(Scan.perturb_Ridge_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Ridge_TargetScan_res_validated <- lapply(seq(Scan.perturb_Ridge_TargetScan_res), function(i) as_data_frame(Scan.perturb_Ridge_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Ridge_ENCORI_res_validated <- lapply(seq(Scan.perturb_Ridge_ENCORI_res), function(i) as_data_frame(Scan.perturb_Ridge_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Phit_NULL_res_validated <- lapply(seq(Scan.perturb_Phit_NULL_res), function(i) as_data_frame(Scan.perturb_Phit_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phit_TargetScan_res_validated <- lapply(seq(Scan.perturb_Phit_TargetScan_res), function(i) as_data_frame(Scan.perturb_Phit_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phit_ENCORI_res_validated <- lapply(seq(Scan.perturb_Phit_ENCORI_res), function(i) as_data_frame(Scan.perturb_Phit_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Phis_NULL_res_validated <- lapply(seq(Scan.perturb_Phis_NULL_res), function(i) as_data_frame(Scan.perturb_Phis_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phis_TargetScan_res_validated <- lapply(seq(Scan.perturb_Phis_TargetScan_res), function(i) as_data_frame(Scan.perturb_Phis_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phis_ENCORI_res_validated <- lapply(seq(Scan.perturb_Phis_ENCORI_res), function(i) as_data_frame(Scan.perturb_Phis_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Rhop_NULL_res_validated <- lapply(seq(Scan.perturb_Rhop_NULL_res), function(i) as_data_frame(Scan.perturb_Rhop_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Rhop_TargetScan_res_validated <- lapply(seq(Scan.perturb_Rhop_TargetScan_res), function(i) as_data_frame(Scan.perturb_Rhop_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Rhop_ENCORI_res_validated <- lapply(seq(Scan.perturb_Rhop_ENCORI_res), function(i) as_data_frame(Scan.perturb_Rhop_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_IDA_NULL_res_validated <- lapply(seq(Scan.perturb_IDA_NULL_res), function(i) as_data_frame(Scan.perturb_IDA_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_IDA_TargetScan_res_validated <- lapply(seq(Scan.perturb_IDA_TargetScan_res), function(i) as_data_frame(Scan.perturb_IDA_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_IDA_ENCORI_res_validated <- lapply(seq(Scan.perturb_IDA_ENCORI_res), function(i) as_data_frame(Scan.perturb_IDA_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.perturb_Pearson_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_NULL_res), function(i) 100*nrow(Scan.perturb_Pearson_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_NULL_res[[i]]))))
Scan.perturb_Pearson_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_TargetScan_res), function(i) 100*nrow(Scan.perturb_Pearson_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_TargetScan_res[[i]]))))
Scan.perturb_Pearson_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_ENCORI_res), function(i) 100*nrow(Scan.perturb_Pearson_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_ENCORI_res[[i]]))))

Scan.perturb_Spearman_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_NULL_res), function(i) 100*nrow(Scan.perturb_Spearman_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_NULL_res[[i]]))))
Scan.perturb_Spearman_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_TargetScan_res), function(i) 100*nrow(Scan.perturb_Spearman_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_TargetScan_res[[i]]))))
Scan.perturb_Spearman_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_ENCORI_res), function(i) 100*nrow(Scan.perturb_Spearman_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_ENCORI_res[[i]]))))

Scan.perturb_Kendall_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_NULL_res), function(i) 100*nrow(Scan.perturb_Kendall_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_NULL_res[[i]]))))
Scan.perturb_Kendall_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_TargetScan_res), function(i) 100*nrow(Scan.perturb_Kendall_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_TargetScan_res[[i]]))))
Scan.perturb_Kendall_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_ENCORI_res), function(i) 100*nrow(Scan.perturb_Kendall_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_ENCORI_res[[i]]))))

Scan.perturb_Dcor_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_NULL_res), function(i) 100*nrow(Scan.perturb_Dcor_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_NULL_res[[i]]))))
Scan.perturb_Dcor_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_TargetScan_res), function(i) 100*nrow(Scan.perturb_Dcor_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_TargetScan_res[[i]]))))
Scan.perturb_Dcor_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_ENCORI_res), function(i) 100*nrow(Scan.perturb_Dcor_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_ENCORI_res[[i]]))))

Scan.perturb_RDC_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_NULL_res), function(i) 100*nrow(Scan.perturb_RDC_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_NULL_res[[i]]))))
Scan.perturb_RDC_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_TargetScan_res), function(i) 100*nrow(Scan.perturb_RDC_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_TargetScan_res[[i]]))))
Scan.perturb_RDC_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_ENCORI_res), function(i) 100*nrow(Scan.perturb_RDC_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_ENCORI_res[[i]]))))

Scan.perturb_Hoeffd_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_NULL_res), function(i) 100*nrow(Scan.perturb_Hoeffd_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_NULL_res[[i]]))))
Scan.perturb_Hoeffd_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_TargetScan_res), function(i) 100*nrow(Scan.perturb_Hoeffd_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_TargetScan_res[[i]]))))
Scan.perturb_Hoeffd_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_ENCORI_res), function(i) 100*nrow(Scan.perturb_Hoeffd_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_ENCORI_res[[i]]))))

Scan.perturb_Zscore_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_NULL_res), function(i) 100*nrow(Scan.perturb_Zscore_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_NULL_res[[i]]))))
Scan.perturb_Zscore_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_TargetScan_res), function(i) 100*nrow(Scan.perturb_Zscore_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_TargetScan_res[[i]]))))
Scan.perturb_Zscore_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_ENCORI_res), function(i) 100*nrow(Scan.perturb_Zscore_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_ENCORI_res[[i]]))))

Scan.perturb_Biweight_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_NULL_res), function(i) 100*nrow(Scan.perturb_Biweight_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_NULL_res[[i]]))))
Scan.perturb_Biweight_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_TargetScan_res), function(i) 100*nrow(Scan.perturb_Biweight_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_TargetScan_res[[i]]))))
Scan.perturb_Biweight_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_ENCORI_res), function(i) 100*nrow(Scan.perturb_Biweight_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_ENCORI_res[[i]]))))

Scan.perturb_Weighted_rank_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_NULL_res), function(i) 100*nrow(Scan.perturb_Weighted_rank_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_NULL_res[[i]]))))
Scan.perturb_Weighted_rank_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res), function(i) 100*nrow(Scan.perturb_Weighted_rank_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res[[i]]))))
Scan.perturb_Weighted_rank_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res), function(i) 100*nrow(Scan.perturb_Weighted_rank_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res[[i]]))))

Scan.perturb_Cosine_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_NULL_res), function(i) 100*nrow(Scan.perturb_Cosine_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_NULL_res[[i]]))))
Scan.perturb_Cosine_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_TargetScan_res), function(i) 100*nrow(Scan.perturb_Cosine_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_TargetScan_res[[i]]))))
Scan.perturb_Cosine_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_ENCORI_res), function(i) 100*nrow(Scan.perturb_Cosine_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_ENCORI_res[[i]]))))

Scan.perturb_Euclidean_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_NULL_res), function(i) 100*nrow(Scan.perturb_Euclidean_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_NULL_res[[i]]))))
Scan.perturb_Euclidean_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_TargetScan_res), function(i) 100*nrow(Scan.perturb_Euclidean_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_TargetScan_res[[i]]))))
Scan.perturb_Euclidean_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_ENCORI_res), function(i) 100*nrow(Scan.perturb_Euclidean_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_ENCORI_res[[i]]))))

Scan.perturb_Manhattan_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_NULL_res), function(i) 100*nrow(Scan.perturb_Manhattan_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_NULL_res[[i]]))))
Scan.perturb_Manhattan_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_TargetScan_res), function(i) 100*nrow(Scan.perturb_Manhattan_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_TargetScan_res[[i]]))))
Scan.perturb_Manhattan_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_ENCORI_res), function(i) 100*nrow(Scan.perturb_Manhattan_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_ENCORI_res[[i]]))))

Scan.perturb_Canberra_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_NULL_res), function(i) 100*nrow(Scan.perturb_Canberra_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_NULL_res[[i]]))))
Scan.perturb_Canberra_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_TargetScan_res), function(i) 100*nrow(Scan.perturb_Canberra_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_TargetScan_res[[i]]))))
Scan.perturb_Canberra_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_ENCORI_res), function(i) 100*nrow(Scan.perturb_Canberra_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_ENCORI_res[[i]]))))

Scan.perturb_Chebyshev_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_NULL_res), function(i) 100*nrow(Scan.perturb_Chebyshev_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_NULL_res[[i]]))))
Scan.perturb_Chebyshev_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_TargetScan_res), function(i) 100*nrow(Scan.perturb_Chebyshev_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_TargetScan_res[[i]]))))
Scan.perturb_Chebyshev_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_ENCORI_res), function(i) 100*nrow(Scan.perturb_Chebyshev_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_ENCORI_res[[i]]))))

Scan.perturb_Dice_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_NULL_res), function(i) 100*nrow(Scan.perturb_Dice_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_NULL_res[[i]]))))
Scan.perturb_Dice_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_TargetScan_res), function(i) 100*nrow(Scan.perturb_Dice_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_TargetScan_res[[i]]))))
Scan.perturb_Dice_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_ENCORI_res), function(i) 100*nrow(Scan.perturb_Dice_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_ENCORI_res[[i]]))))

Scan.perturb_Jaccard_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_NULL_res), function(i) 100*nrow(Scan.perturb_Jaccard_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_NULL_res[[i]]))))
Scan.perturb_Jaccard_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_TargetScan_res), function(i) 100*nrow(Scan.perturb_Jaccard_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_TargetScan_res[[i]]))))
Scan.perturb_Jaccard_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_ENCORI_res), function(i) 100*nrow(Scan.perturb_Jaccard_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_ENCORI_res[[i]]))))

Scan.perturb_Mahalanobis_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_NULL_res), function(i) 100*nrow(Scan.perturb_Mahalanobis_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_NULL_res[[i]]))))
Scan.perturb_Mahalanobis_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res), function(i) 100*nrow(Scan.perturb_Mahalanobis_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res[[i]]))))
Scan.perturb_Mahalanobis_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res), function(i) 100*nrow(Scan.perturb_Mahalanobis_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res[[i]]))))

Scan.perturb_MI_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_MI_NULL_res), function(i) 100*nrow(Scan.perturb_MI_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_NULL_res[[i]]))))
Scan.perturb_MI_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_MI_TargetScan_res), function(i) 100*nrow(Scan.perturb_MI_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_TargetScan_res[[i]]))))
Scan.perturb_MI_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_MI_ENCORI_res), function(i) 100*nrow(Scan.perturb_MI_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_ENCORI_res[[i]]))))

Scan.perturb_MIC_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_NULL_res), function(i) 100*nrow(Scan.perturb_MIC_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_NULL_res[[i]]))))
Scan.perturb_MIC_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_TargetScan_res), function(i) 100*nrow(Scan.perturb_MIC_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_TargetScan_res[[i]]))))
Scan.perturb_MIC_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_ENCORI_res), function(i) 100*nrow(Scan.perturb_MIC_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_ENCORI_res[[i]]))))

Scan.perturb_Lasso_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_NULL_res), function(i) 100*nrow(Scan.perturb_Lasso_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_NULL_res[[i]]))))
Scan.perturb_Lasso_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_TargetScan_res), function(i) 100*nrow(Scan.perturb_Lasso_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_TargetScan_res[[i]]))))
Scan.perturb_Lasso_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_ENCORI_res), function(i) 100*nrow(Scan.perturb_Lasso_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_ENCORI_res[[i]]))))

Scan.perturb_Elastic_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_NULL_res), function(i) 100*nrow(Scan.perturb_Elastic_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_NULL_res[[i]]))))
Scan.perturb_Elastic_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_TargetScan_res), function(i) 100*nrow(Scan.perturb_Elastic_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_TargetScan_res[[i]]))))
Scan.perturb_Elastic_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_ENCORI_res), function(i) 100*nrow(Scan.perturb_Elastic_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_ENCORI_res[[i]]))))

Scan.perturb_Ridge_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_NULL_res), function(i) 100*nrow(Scan.perturb_Ridge_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_NULL_res[[i]]))))
Scan.perturb_Ridge_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_TargetScan_res), function(i) 100*nrow(Scan.perturb_Ridge_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_TargetScan_res[[i]]))))
Scan.perturb_Ridge_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_ENCORI_res), function(i) 100*nrow(Scan.perturb_Ridge_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_ENCORI_res[[i]]))))

Scan.perturb_Phit_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_NULL_res), function(i) 100*nrow(Scan.perturb_Phit_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_NULL_res[[i]]))))
Scan.perturb_Phit_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_TargetScan_res), function(i) 100*nrow(Scan.perturb_Phit_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_TargetScan_res[[i]]))))
Scan.perturb_Phit_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_ENCORI_res), function(i) 100*nrow(Scan.perturb_Phit_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_ENCORI_res[[i]]))))

Scan.perturb_Phis_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_NULL_res), function(i) 100*nrow(Scan.perturb_Phis_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_NULL_res[[i]]))))
Scan.perturb_Phis_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_TargetScan_res), function(i) 100*nrow(Scan.perturb_Phis_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_TargetScan_res[[i]]))))
Scan.perturb_Phis_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_ENCORI_res), function(i) 100*nrow(Scan.perturb_Phis_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_ENCORI_res[[i]]))))

Scan.perturb_Rhop_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_NULL_res), function(i) 100*nrow(Scan.perturb_Rhop_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_NULL_res[[i]]))))
Scan.perturb_Rhop_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_TargetScan_res), function(i) 100*nrow(Scan.perturb_Rhop_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_TargetScan_res[[i]]))))
Scan.perturb_Rhop_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_ENCORI_res), function(i) 100*nrow(Scan.perturb_Rhop_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_ENCORI_res[[i]]))))

Scan.perturb_IDA_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_NULL_res), function(i) 100*nrow(Scan.perturb_IDA_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_NULL_res[[i]]))))
Scan.perturb_IDA_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_TargetScan_res), function(i) 100*nrow(Scan.perturb_IDA_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_TargetScan_res[[i]]))))
Scan.perturb_IDA_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_ENCORI_res), function(i) 100*nrow(Scan.perturb_IDA_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_ENCORI_res[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.perturb_Pearson_NULL_res_Sim <- Sim.network(Scan.perturb_Pearson_NULL_res)
Scan.perturb_Pearson_TargetScan_res_Sim <- Sim.network(Scan.perturb_Pearson_TargetScan_res)
Scan.perturb_Pearson_ENCORI_res_Sim <- Sim.network(Scan.perturb_Pearson_ENCORI_res)

Scan.perturb_Spearman_NULL_res_Sim <- Sim.network(Scan.perturb_Spearman_NULL_res)
Scan.perturb_Spearman_TargetScan_res_Sim <- Sim.network(Scan.perturb_Spearman_TargetScan_res)
Scan.perturb_Spearman_ENCORI_res_Sim <- Sim.network(Scan.perturb_Spearman_ENCORI_res)

Scan.perturb_Kendall_NULL_res_Sim <- Sim.network(Scan.perturb_Kendall_NULL_res)
Scan.perturb_Kendall_TargetScan_res_Sim <- Sim.network(Scan.perturb_Kendall_TargetScan_res)
Scan.perturb_Kendall_ENCORI_res_Sim <- Sim.network(Scan.perturb_Kendall_ENCORI_res)

Scan.perturb_Dcor_NULL_res_Sim <- Sim.network(Scan.perturb_Dcor_NULL_res)
Scan.perturb_Dcor_TargetScan_res_Sim <- Sim.network(Scan.perturb_Dcor_TargetScan_res)
Scan.perturb_Dcor_ENCORI_res_Sim <- Sim.network(Scan.perturb_Dcor_ENCORI_res)

Scan.perturb_RDC_NULL_res_Sim <- Sim.network(Scan.perturb_RDC_NULL_res)
Scan.perturb_RDC_TargetScan_res_Sim <- Sim.network(Scan.perturb_RDC_TargetScan_res)
Scan.perturb_RDC_ENCORI_res_Sim <- Sim.network(Scan.perturb_RDC_ENCORI_res)

Scan.perturb_Hoeffd_NULL_res_Sim <- Sim.network(Scan.perturb_Hoeffd_NULL_res)
Scan.perturb_Hoeffd_TargetScan_res_Sim <- Sim.network(Scan.perturb_Hoeffd_TargetScan_res)
Scan.perturb_Hoeffd_ENCORI_res_Sim <- Sim.network(Scan.perturb_Hoeffd_ENCORI_res)

Scan.perturb_Zscore_NULL_res_Sim <- Sim.network(Scan.perturb_Zscore_NULL_res)
Scan.perturb_Zscore_TargetScan_res_Sim <- Sim.network(Scan.perturb_Zscore_TargetScan_res)
Scan.perturb_Zscore_ENCORI_res_Sim <- Sim.network(Scan.perturb_Zscore_ENCORI_res)

Scan.perturb_Biweight_NULL_res_Sim <- Sim.network(Scan.perturb_Biweight_NULL_res)
Scan.perturb_Biweight_TargetScan_res_Sim <- Sim.network(Scan.perturb_Biweight_TargetScan_res)
Scan.perturb_Biweight_ENCORI_res_Sim <- Sim.network(Scan.perturb_Biweight_ENCORI_res)

Scan.perturb_Weighted_rank_NULL_res_Sim <- Sim.network(Scan.perturb_Weighted_rank_NULL_res)
Scan.perturb_Weighted_rank_TargetScan_res_Sim <- Sim.network(Scan.perturb_Weighted_rank_TargetScan_res)
Scan.perturb_Weighted_rank_ENCORI_res_Sim <- Sim.network(Scan.perturb_Weighted_rank_ENCORI_res)

Scan.perturb_Cosine_NULL_res_Sim <- Sim.network(Scan.perturb_Cosine_NULL_res)
Scan.perturb_Cosine_TargetScan_res_Sim <- Sim.network(Scan.perturb_Cosine_TargetScan_res)
Scan.perturb_Cosine_ENCORI_res_Sim <- Sim.network(Scan.perturb_Cosine_ENCORI_res)

Scan.perturb_Euclidean_NULL_res_Sim <- Sim.network(Scan.perturb_Euclidean_NULL_res)
Scan.perturb_Euclidean_TargetScan_res_Sim <- Sim.network(Scan.perturb_Euclidean_TargetScan_res)
Scan.perturb_Euclidean_ENCORI_res_Sim <- Sim.network(Scan.perturb_Euclidean_ENCORI_res)

Scan.perturb_Manhattan_NULL_res_Sim <- Sim.network(Scan.perturb_Manhattan_NULL_res)
Scan.perturb_Manhattan_TargetScan_res_Sim <- Sim.network(Scan.perturb_Manhattan_TargetScan_res)
Scan.perturb_Manhattan_ENCORI_res_Sim <- Sim.network(Scan.perturb_Manhattan_ENCORI_res)

Scan.perturb_Canberra_NULL_res_Sim <- Sim.network(Scan.perturb_Canberra_NULL_res)
Scan.perturb_Canberra_TargetScan_res_Sim <- Sim.network(Scan.perturb_Canberra_TargetScan_res)
Scan.perturb_Canberra_ENCORI_res_Sim <- Sim.network(Scan.perturb_Canberra_ENCORI_res)

Scan.perturb_Chebyshev_NULL_res_Sim <- Sim.network(Scan.perturb_Chebyshev_NULL_res)
Scan.perturb_Chebyshev_TargetScan_res_Sim <- Sim.network(Scan.perturb_Chebyshev_TargetScan_res)
Scan.perturb_Chebyshev_ENCORI_res_Sim <- Sim.network(Scan.perturb_Chebyshev_ENCORI_res)

Scan.perturb_Dice_NULL_res_Sim <- Sim.network(Scan.perturb_Dice_NULL_res)
Scan.perturb_Dice_TargetScan_res_Sim <- Sim.network(Scan.perturb_Dice_TargetScan_res)
Scan.perturb_Dice_ENCORI_res_Sim <- Sim.network(Scan.perturb_Dice_ENCORI_res)

Scan.perturb_Jaccard_NULL_res_Sim <- Sim.network(Scan.perturb_Jaccard_NULL_res)
Scan.perturb_Jaccard_TargetScan_res_Sim <- Sim.network(Scan.perturb_Jaccard_TargetScan_res)
Scan.perturb_Jaccard_ENCORI_res_Sim <- Sim.network(Scan.perturb_Jaccard_ENCORI_res)

Scan.perturb_Mahalanobis_NULL_res_Sim <- Sim.network(Scan.perturb_Mahalanobis_NULL_res)
Scan.perturb_Mahalanobis_TargetScan_res_Sim <- Sim.network(Scan.perturb_Mahalanobis_TargetScan_res)
Scan.perturb_Mahalanobis_ENCORI_res_Sim <- Sim.network(Scan.perturb_Mahalanobis_ENCORI_res)

Scan.perturb_MI_NULL_res_Sim <- Sim.network(Scan.perturb_MI_NULL_res)
Scan.perturb_MI_TargetScan_res_Sim <- Sim.network(Scan.perturb_MI_TargetScan_res)
Scan.perturb_MI_ENCORI_res_Sim <- Sim.network(Scan.perturb_MI_ENCORI_res)

Scan.perturb_MIC_NULL_res_Sim <- Sim.network(Scan.perturb_MIC_NULL_res)
Scan.perturb_MIC_TargetScan_res_Sim <- Sim.network(Scan.perturb_MIC_TargetScan_res)
Scan.perturb_MIC_ENCORI_res_Sim <- Sim.network(Scan.perturb_MIC_ENCORI_res)

Scan.perturb_Lasso_NULL_res_Sim <- Sim.network(Scan.perturb_Lasso_NULL_res)
Scan.perturb_Lasso_TargetScan_res_Sim <- Sim.network(Scan.perturb_Lasso_TargetScan_res)
Scan.perturb_Lasso_ENCORI_res_Sim <- Sim.network(Scan.perturb_Lasso_ENCORI_res)

Scan.perturb_Elastic_NULL_res_Sim <- Sim.network(Scan.perturb_Elastic_NULL_res)
Scan.perturb_Elastic_TargetScan_res_Sim <- Sim.network(Scan.perturb_Elastic_TargetScan_res)
Scan.perturb_Elastic_ENCORI_res_Sim <- Sim.network(Scan.perturb_Elastic_ENCORI_res)

Scan.perturb_Ridge_NULL_res_Sim <- Sim.network(Scan.perturb_Ridge_NULL_res)
Scan.perturb_Ridge_TargetScan_res_Sim <- Sim.network(Scan.perturb_Ridge_TargetScan_res)
Scan.perturb_Ridge_ENCORI_res_Sim <- Sim.network(Scan.perturb_Ridge_ENCORI_res)

Scan.perturb_Phit_NULL_res_Sim <- Sim.network(Scan.perturb_Phit_NULL_res)
Scan.perturb_Phit_TargetScan_res_Sim <- Sim.network(Scan.perturb_Phit_TargetScan_res)
Scan.perturb_Phit_ENCORI_res_Sim <- Sim.network(Scan.perturb_Phit_ENCORI_res)

Scan.perturb_Phis_NULL_res_Sim <- Sim.network(Scan.perturb_Phis_NULL_res)
Scan.perturb_Phis_TargetScan_res_Sim <- Sim.network(Scan.perturb_Phis_TargetScan_res)
Scan.perturb_Phis_ENCORI_res_Sim <- Sim.network(Scan.perturb_Phis_ENCORI_res)

Scan.perturb_Rhop_NULL_res_Sim <- Sim.network(Scan.perturb_Rhop_NULL_res)
Scan.perturb_Rhop_TargetScan_res_Sim <- Sim.network(Scan.perturb_Rhop_TargetScan_res)
Scan.perturb_Rhop_ENCORI_res_Sim <- Sim.network(Scan.perturb_Rhop_ENCORI_res)

Scan.perturb_IDA_NULL_res_Sim <- Sim.network(Scan.perturb_IDA_NULL_res)
Scan.perturb_IDA_TargetScan_res_Sim <- Sim.network(Scan.perturb_IDA_TargetScan_res)
Scan.perturb_IDA_ENCORI_res_Sim <- Sim.network(Scan.perturb_IDA_ENCORI_res)

save.image("Scan.perturb_K562.RData")

#######################################################################################################################################################################################
############################################################################ Scan.perturb application in BRCA dataset #################################################################
#######################################################################################################################################################################################

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

set.seed(123)
ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

# No prior information
Scan.perturb_Pearson_timestart_BRCA <- Sys.time()
Scan.perturb_Pearson_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Pearson")
Scan.perturb_Pearson_timeend_BRCA <- Sys.time()
Scan.perturb_Pearson_runningtime_NULL_BRCA <- Scan.perturb_Pearson_timeend_BRCA - Scan.perturb_Pearson_timestart_BRCA

Scan.perturb_Spearman_timestart_BRCA <- Sys.time()
Scan.perturb_Spearman_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Spearman")
Scan.perturb_Spearman_timeend_BRCA <- Sys.time()
Scan.perturb_Spearman_runningtime_NULL_BRCA <- Scan.perturb_Spearman_timeend_BRCA - Scan.perturb_Spearman_timestart_BRCA

Scan.perturb_Kendall_timestart_BRCA <- Sys.time()
Scan.perturb_Kendall_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Kendall")
Scan.perturb_Kendall_timeend_BRCA <- Sys.time()
Scan.perturb_Kendall_runningtime_NULL_BRCA <- Scan.perturb_Kendall_timeend_BRCA - Scan.perturb_Kendall_timestart_BRCA

Scan.perturb_Dcor_timestart_BRCA <- Sys.time()
Scan.perturb_Dcor_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Dcor")
Scan.perturb_Dcor_timeend_BRCA <- Sys.time()
Scan.perturb_Dcor_runningtime_NULL_BRCA <- Scan.perturb_Dcor_timeend_BRCA - Scan.perturb_Dcor_timestart_BRCA

Scan.perturb_RDC_timestart_BRCA <- Sys.time()
Scan.perturb_RDC_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "RDC")
Scan.perturb_RDC_timeend_BRCA <- Sys.time()
Scan.perturb_RDC_runningtime_NULL_BRCA <- Scan.perturb_RDC_timeend_BRCA - Scan.perturb_RDC_timestart_BRCA

Scan.perturb_Hoeffd_timestart_BRCA <- Sys.time()
Scan.perturb_Hoeffd_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Hoeffd")
Scan.perturb_Hoeffd_timeend_BRCA <- Sys.time()
Scan.perturb_Hoeffd_runningtime_NULL_BRCA <- Scan.perturb_Hoeffd_timeend_BRCA - Scan.perturb_Hoeffd_timestart_BRCA

Scan.perturb_Zscore_timestart_BRCA <- Sys.time()
Scan.perturb_Zscore_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Zscore")
Scan.perturb_Zscore_timeend_BRCA <- Sys.time()
Scan.perturb_Zscore_runningtime_NULL_BRCA <- Scan.perturb_Zscore_timeend_BRCA - Scan.perturb_Zscore_timestart_BRCA

Scan.perturb_Biweight_timestart_BRCA <- Sys.time()
Scan.perturb_Biweight_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Biweight")
Scan.perturb_Biweight_timeend_BRCA <- Sys.time()
Scan.perturb_Biweight_runningtime_NULL_BRCA <- Scan.perturb_Biweight_timeend_BRCA - Scan.perturb_Biweight_timestart_BRCA

Scan.perturb_Weighted_rank_timestart_BRCA <- Sys.time()
Scan.perturb_Weighted_rank_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Weighted_rank")
Scan.perturb_Weighted_rank_timeend_BRCA <- Sys.time()
Scan.perturb_Weighted_rank_runningtime_NULL_BRCA <- Scan.perturb_Weighted_rank_timeend_BRCA - Scan.perturb_Weighted_rank_timestart_BRCA

Scan.perturb_Cosine_timestart_BRCA <- Sys.time()
Scan.perturb_Cosine_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Cosine")
Scan.perturb_Cosine_timeend_BRCA <- Sys.time()
Scan.perturb_Cosine_runningtime_NULL_BRCA <- Scan.perturb_Cosine_timeend_BRCA - Scan.perturb_Cosine_timestart_BRCA

Scan.perturb_Euclidean_timestart_BRCA <- Sys.time()
Scan.perturb_Euclidean_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Euclidean")
Scan.perturb_Euclidean_timeend_BRCA <- Sys.time()
Scan.perturb_Euclidean_runningtime_NULL_BRCA <- Scan.perturb_Euclidean_timeend_BRCA - Scan.perturb_Euclidean_timestart_BRCA

Scan.perturb_Manhattan_timestart_BRCA <- Sys.time()
Scan.perturb_Manhattan_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Manhattan")
Scan.perturb_Manhattan_timeend_BRCA <- Sys.time()
Scan.perturb_Manhattan_runningtime_NULL_BRCA <- Scan.perturb_Manhattan_timeend_BRCA - Scan.perturb_Manhattan_timestart_BRCA

Scan.perturb_Canberra_timestart_BRCA <- Sys.time()
Scan.perturb_Canberra_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Canberra")
Scan.perturb_Canberra_timeend_BRCA <- Sys.time()
Scan.perturb_Canberra_runningtime_NULL_BRCA <- Scan.perturb_Canberra_timeend_BRCA - Scan.perturb_Canberra_timestart_BRCA

Scan.perturb_Chebyshev_timestart_BRCA <- Sys.time()
Scan.perturb_Chebyshev_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Chebyshev")
Scan.perturb_Chebyshev_timeend_BRCA <- Sys.time()
Scan.perturb_Chebyshev_runningtime_NULL_BRCA <- Scan.perturb_Chebyshev_timeend_BRCA - Scan.perturb_Chebyshev_timestart_BRCA

Scan.perturb_Dice_timestart_BRCA <- Sys.time()
Scan.perturb_Dice_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Dice")
Scan.perturb_Dice_timeend_BRCA <- Sys.time()
Scan.perturb_Dice_runningtime_NULL_BRCA <- Scan.perturb_Dice_timeend_BRCA - Scan.perturb_Dice_timestart_BRCA

Scan.perturb_Jaccard_timestart_BRCA <- Sys.time()
Scan.perturb_Jaccard_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Jaccard")
Scan.perturb_Jaccard_timeend_BRCA <- Sys.time()
Scan.perturb_Jaccard_runningtime_NULL_BRCA <- Scan.perturb_Jaccard_timeend_BRCA - Scan.perturb_Jaccard_timestart_BRCA

Scan.perturb_Mahalanobis_timestart_BRCA <- Sys.time()
Scan.perturb_Mahalanobis_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Mahalanobis")
Scan.perturb_Mahalanobis_timeend_BRCA <- Sys.time()
Scan.perturb_Mahalanobis_runningtime_NULL_BRCA <- Scan.perturb_Mahalanobis_timeend_BRCA - Scan.perturb_Mahalanobis_timestart_BRCA

Scan.perturb_MI_timestart_BRCA <- Sys.time()
Scan.perturb_MI_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "MI")
Scan.perturb_MI_timeend_BRCA <- Sys.time()
Scan.perturb_MI_runningtime_NULL_BRCA <- Scan.perturb_MI_timeend_BRCA - Scan.perturb_MI_timestart_BRCA

Scan.perturb_MIC_timestart_BRCA <- Sys.time()
Scan.perturb_MIC_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "MIC")
Scan.perturb_MIC_timeend_BRCA <- Sys.time()
Scan.perturb_MIC_runningtime_NULL_BRCA <- Scan.perturb_MIC_timeend_BRCA - Scan.perturb_MIC_timestart_BRCA

Scan.perturb_Lasso_timestart_BRCA <- Sys.time()
Scan.perturb_Lasso_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Lasso")
Scan.perturb_Lasso_timeend_BRCA <- Sys.time()
Scan.perturb_Lasso_runningtime_NULL_BRCA <- Scan.perturb_Lasso_timeend_BRCA - Scan.perturb_Lasso_timestart_BRCA

Scan.perturb_Elastic_timestart_BRCA <- Sys.time()
Scan.perturb_Elastic_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Elastic")
Scan.perturb_Elastic_timeend_BRCA <- Sys.time()
Scan.perturb_Elastic_runningtime_NULL_BRCA <- Scan.perturb_Elastic_timeend_BRCA - Scan.perturb_Elastic_timestart_BRCA

Scan.perturb_Ridge_timestart_BRCA <- Sys.time()
Scan.perturb_Ridge_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Ridge")
Scan.perturb_Ridge_timeend_BRCA <- Sys.time()
Scan.perturb_Ridge_runningtime_NULL_BRCA <- Scan.perturb_Ridge_timeend_BRCA - Scan.perturb_Ridge_timestart_BRCA

Scan.perturb_Phit_timestart_BRCA <- Sys.time()
Scan.perturb_Phit_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Phit")
Scan.perturb_Phit_timeend_BRCA <- Sys.time()
Scan.perturb_Phit_runningtime_NULL_BRCA <- Scan.perturb_Phit_timeend_BRCA - Scan.perturb_Phit_timestart_BRCA

Scan.perturb_Phis_timestart_BRCA <- Sys.time()
Scan.perturb_Phis_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Phis")
Scan.perturb_Phis_timeend_BRCA <- Sys.time()
Scan.perturb_Phis_runningtime_NULL_BRCA <- Scan.perturb_Phis_timeend_BRCA - Scan.perturb_Phis_timestart_BRCA

Scan.perturb_Rhop_timestart_BRCA <- Sys.time()
Scan.perturb_Rhop_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "Rhop")
Scan.perturb_Rhop_timeend_BRCA <- Sys.time()
Scan.perturb_Rhop_runningtime_NULL_BRCA <- Scan.perturb_Rhop_timeend_BRCA - Scan.perturb_Rhop_timestart_BRCA

Scan.perturb_IDA_timestart_BRCA <- Sys.time()
Scan.perturb_IDA_NULL_res_BRCA <- Scan.perturb(BRCA_miRNA_Exp_DEG[[2]], BRCA_mRNA_Exp_DEG[[2]], method = "IDA", pcmethod = "stable)
Scan.perturb_IDA_timeend_BRCA <- Sys.time()
Scan.perturb_IDA_runningtime_NULL_BRCA <- Scan.perturb_IDAp_timeend_BRCA - Scan.perturb_IDA_timestart_BRCA

# The prior information is TargetScan
Scan.perturb_Pearson_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Pearson_NULL_res_BRCA), function(i) Scan.perturb_Pearson_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Spearman_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Spearman_NULL_res_BRCA), function(i) Scan.perturb_Spearman_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Kendall_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Kendall_NULL_res_BRCA), function(i) Scan.perturb_Kendall_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Dcor_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Dcor_NULL_res_BRCA), function(i) Scan.perturb_Dcor_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_RDC_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_RDC_NULL_res_BRCA), function(i) Scan.perturb_RDC_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Hoeffd_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Hoeffd_NULL_res_BRCA), function(i) Scan.perturb_Hoeffd_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Zscore_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Zscore_NULL_res_BRCA), function(i) Scan.perturb_Zscore_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Biweight_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Biweight_NULL_res_BRCA), function(i) Scan.perturb_Biweight_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Weighted_rank_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res_BRCA), function(i) Scan.perturb_Weighted_rank_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Cosine_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Cosine_NULL_res_BRCA), function(i) Scan.perturb_Cosine_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Euclidean_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Euclidean_NULL_res_BRCA), function(i) Scan.perturb_Euclidean_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Manhattan_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Manhattan_NULL_res_BRCA), function(i) Scan.perturb_Manhattan_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Canberra_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Canberra_NULL_res_BRCA), function(i) Scan.perturb_Canberra_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Chebyshev_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Chebyshev_NULL_res_BRCA), function(i) Scan.perturb_Chebyshev_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Dice_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Dice_NULL_res_BRCA), function(i) Scan.perturb_Dice_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Jaccard_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Jaccard_NULL_res_BRCA), function(i) Scan.perturb_Jaccard_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Mahalanobis_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res_BRCA), function(i) Scan.perturb_Mahalanobis_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_MI_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_MI_NULL_res_BRCA), function(i) Scan.perturb_MI_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_MIC_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_MIC_NULL_res_BRCA), function(i) Scan.perturb_MIC_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Lasso_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Lasso_NULL_res_BRCA), function(i) Scan.perturb_Lasso_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Elastic_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Elastic_NULL_res_BRCA), function(i) Scan.perturb_Elastic_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Ridge_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Ridge_NULL_res_BRCA), function(i) Scan.perturb_Ridge_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Phit_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Phit_NULL_res_BRCA), function(i) Scan.perturb_Phit_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Phis_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Phis_NULL_res_BRCA), function(i) Scan.perturb_Phis_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_Rhop_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_Rhop_NULL_res_BRCA), function(i) Scan.perturb_Rhop_NULL_res_BRCA[[i]] %s% TargetScan_graph)
Scan.perturb_IDA_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_IDA_NULL_res_BRCA), function(i) Scan.perturb_IDA_NULL_res_BRCA[[i]] %s% TargetScan_graph)

# The prior information is ENCORI
Scan.perturb_Pearson_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Pearson_NULL_res_BRCA), function(i) Scan.perturb_Pearson_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Spearman_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Spearman_NULL_res_BRCA), function(i) Scan.perturb_Spearman_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Kendall_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Kendall_NULL_res_BRCA), function(i) Scan.perturb_Kendall_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Dcor_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Dcor_NULL_res_BRCA), function(i) Scan.perturb_Dcor_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_RDC_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_RDC_NULL_res_BRCA), function(i) Scan.perturb_RDC_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Hoeffd_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Hoeffd_NULL_res_BRCA), function(i) Scan.perturb_Hoeffd_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Zscore_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Zscore_NULL_res_BRCA), function(i) Scan.perturb_Zscore_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Biweight_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Biweight_NULL_res_BRCA), function(i) Scan.perturb_Biweight_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Weighted_rank_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res_BRCA), function(i) Scan.perturb_Weighted_rank_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Cosine_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Cosine_NULL_res_BRCA), function(i) Scan.perturb_Cosine_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Euclidean_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Euclidean_NULL_res_BRCA), function(i) Scan.perturb_Euclidean_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Manhattan_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Manhattan_NULL_res_BRCA), function(i) Scan.perturb_Manhattan_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Canberra_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Canberra_NULL_res_BRCA), function(i) Scan.perturb_Canberra_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Chebyshev_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Chebyshev_NULL_res_BRCA), function(i) Scan.perturb_Chebyshev_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Dice_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Dice_NULL_res_BRCA), function(i) Scan.perturb_Dice_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Jaccard_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Jaccard_NULL_res_BRCA), function(i) Scan.perturb_Jaccard_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Mahalanobis_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res_BRCA), function(i) Scan.perturb_Mahalanobis_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_MI_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_MI_NULL_res_BRCA), function(i) Scan.perturb_MI_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_MIC_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_MIC_NULL_res_BRCA), function(i) Scan.perturb_MIC_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Lasso_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Lasso_NULL_res_BRCA), function(i) Scan.perturb_Lasso_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Elastic_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Elastic_NULL_res_BRCA), function(i) Scan.perturb_Elastic_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Ridge_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Ridge_NULL_res_BRCA), function(i) Scan.perturb_Ridge_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Phit_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Phit_NULL_res_BRCA), function(i) Scan.perturb_Phit_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Phis_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Phis_NULL_res_BRCA), function(i) Scan.perturb_Phis_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_Rhop_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_Rhop_NULL_res_BRCA), function(i) Scan.perturb_Rhop_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_IDA_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_IDA_NULL_res_BRCA), function(i) Scan.perturb_IDA_NULL_res_BRCA[[i]] %s% ENCORI_graph)

# Number of predicted sample-specific interactions
Scan.perturb_Pearson_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Pearson_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Pearson_NULL_res_BRCA[[i]] ))))
Scan.perturb_Pearson_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Pearson_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Pearson_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Pearson_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Pearson_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Pearson_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Spearman_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Spearman_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Spearman_NULL_res_BRCA[[i]] ))))
Scan.perturb_Spearman_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Spearman_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Spearman_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Spearman_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Spearman_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Spearman_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Kendall_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Kendall_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Kendall_NULL_res_BRCA[[i]] ))))
Scan.perturb_Kendall_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Kendall_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Kendall_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Kendall_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Kendall_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Kendall_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Dcor_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dcor_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dcor_NULL_res_BRCA[[i]] ))))
Scan.perturb_Dcor_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dcor_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dcor_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Dcor_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dcor_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dcor_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_RDC_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_RDC_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_RDC_NULL_res_BRCA[[i]] ))))
Scan.perturb_RDC_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_RDC_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_RDC_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_RDC_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_RDC_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_RDC_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Hoeffd_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_NULL_res_BRCA[[i]] ))))
Scan.perturb_Hoeffd_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Hoeffd_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Hoeffd_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Hoeffd_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Zscore_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Zscore_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Zscore_NULL_res_BRCA[[i]] ))))
Scan.perturb_Zscore_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Zscore_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Zscore_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Zscore_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Zscore_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Zscore_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Biweight_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Biweight_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Biweight_NULL_res_BRCA[[i]] ))))
Scan.perturb_Biweight_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Biweight_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Biweight_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Biweight_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Biweight_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Biweight_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Weighted_rank_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_NULL_res_BRCA[[i]] ))))
Scan.perturb_Weighted_rank_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Weighted_rank_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Cosine_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Cosine_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Cosine_NULL_res_BRCA[[i]] ))))
Scan.perturb_Cosine_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Cosine_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Cosine_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Cosine_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Cosine_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Cosine_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Euclidean_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Euclidean_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_NULL_res_BRCA[[i]] ))))
Scan.perturb_Euclidean_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Euclidean_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Euclidean_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Euclidean_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Euclidean_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Manhattan_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Manhattan_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_NULL_res_BRCA[[i]] ))))
Scan.perturb_Manhattan_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Manhattan_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Manhattan_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Manhattan_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Manhattan_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Canberra_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Canberra_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Canberra_NULL_res_BRCA[[i]] ))))
Scan.perturb_Canberra_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Canberra_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Canberra_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Canberra_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Canberra_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Canberra_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Chebyshev_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_NULL_res_BRCA[[i]] ))))
Scan.perturb_Chebyshev_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Chebyshev_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Chebyshev_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Chebyshev_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Dice_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dice_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dice_NULL_res_BRCA[[i]] ))))
Scan.perturb_Dice_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dice_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dice_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Dice_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Dice_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Dice_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Jaccard_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Jaccard_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_NULL_res_BRCA[[i]] ))))
Scan.perturb_Jaccard_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Jaccard_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Jaccard_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Jaccard_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Jaccard_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Mahalanobis_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_NULL_res_BRCA[[i]] ))))
Scan.perturb_Mahalanobis_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Mahalanobis_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_MI_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MI_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MI_NULL_res_BRCA[[i]] ))))
Scan.perturb_MI_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MI_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MI_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_MI_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MI_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MI_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_MIC_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MIC_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MIC_NULL_res_BRCA[[i]] ))))
Scan.perturb_MIC_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MIC_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MIC_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_MIC_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_MIC_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_MIC_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Lasso_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Lasso_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Lasso_NULL_res_BRCA[[i]] ))))
Scan.perturb_Lasso_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Lasso_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Lasso_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Lasso_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Lasso_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Lasso_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Elastic_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Elastic_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Elastic_NULL_res_BRCA[[i]] ))))
Scan.perturb_Elastic_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Elastic_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Elastic_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Elastic_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Elastic_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Elastic_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Ridge_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Ridge_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Ridge_NULL_res_BRCA[[i]] ))))
Scan.perturb_Ridge_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Ridge_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Ridge_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Ridge_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Ridge_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Ridge_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Phit_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phit_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phit_NULL_res_BRCA[[i]] ))))
Scan.perturb_Phit_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phit_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phit_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Phit_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phit_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phit_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Phis_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phis_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phis_NULL_res_BRCA[[i]] ))))
Scan.perturb_Phis_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phis_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phis_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Phis_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Phis_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Phis_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_Rhop_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Rhop_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Rhop_NULL_res_BRCA[[i]] ))))
Scan.perturb_Rhop_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Rhop_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Rhop_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_Rhop_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_Rhop_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_Rhop_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_IDA_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_IDA_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_IDA_NULL_res_BRCA[[i]] ))))
Scan.perturb_IDA_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_IDA_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_IDA_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_IDA_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_IDA_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_IDA_ENCORI_res_BRCA[[i]] ))))

# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.perturb_Pearson_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Pearson_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Pearson_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Pearson_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Pearson_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Pearson_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Pearson_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Pearson_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Pearson_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Spearman_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Spearman_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Spearman_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Spearman_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Spearman_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Spearman_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Spearman_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Spearman_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Spearman_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Kendall_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Kendall_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Kendall_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Kendall_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Kendall_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Kendall_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Kendall_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Kendall_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Kendall_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Dcor_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Dcor_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Dcor_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dcor_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Dcor_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Dcor_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dcor_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Dcor_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Dcor_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_RDC_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_RDC_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_RDC_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_RDC_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_RDC_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_RDC_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_RDC_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_RDC_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_RDC_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Hoeffd_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Hoeffd_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Hoeffd_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Hoeffd_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Hoeffd_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Hoeffd_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Hoeffd_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Hoeffd_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Hoeffd_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Zscore_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Zscore_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Zscore_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Zscore_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Zscore_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Zscore_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Zscore_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Zscore_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Zscore_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Biweight_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Biweight_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Biweight_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Biweight_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Biweight_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Biweight_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Biweight_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Biweight_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Biweight_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Weighted_rank_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Weighted_rank_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Weighted_rank_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Weighted_rank_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Weighted_rank_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Cosine_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Cosine_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Cosine_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Cosine_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Cosine_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Cosine_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Cosine_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Cosine_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Cosine_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Euclidean_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Euclidean_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Euclidean_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Euclidean_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Euclidean_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Euclidean_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Euclidean_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Euclidean_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Euclidean_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Manhattan_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Manhattan_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Manhattan_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Manhattan_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Manhattan_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Manhattan_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Manhattan_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Manhattan_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Manhattan_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Canberra_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Canberra_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Canberra_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Canberra_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Canberra_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Canberra_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Canberra_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Canberra_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Canberra_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Chebyshev_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Chebyshev_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Chebyshev_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Chebyshev_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Chebyshev_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Chebyshev_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Chebyshev_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Chebyshev_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Chebyshev_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Dice_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Dice_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Dice_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dice_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Dice_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Dice_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Dice_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Dice_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Dice_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Jaccard_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Jaccard_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Jaccard_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Jaccard_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Jaccard_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Jaccard_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Jaccard_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Jaccard_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Jaccard_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Mahalanobis_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Mahalanobis_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Mahalanobis_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Mahalanobis_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Mahalanobis_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_MI_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_MI_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_MI_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MI_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_MI_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_MI_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MI_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_MI_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_MI_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_MIC_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_MIC_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_MIC_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MIC_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_MIC_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_MIC_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_MIC_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_MIC_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_MIC_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Lasso_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Lasso_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Lasso_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Lasso_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Lasso_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Lasso_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Lasso_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Lasso_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Lasso_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Elastic_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Elastic_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Elastic_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Elastic_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Elastic_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Elastic_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Elastic_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Elastic_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Elastic_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Ridge_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Ridge_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Ridge_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Ridge_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Ridge_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Ridge_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Ridge_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Ridge_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Ridge_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Phit_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Phit_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Phit_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phit_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Phit_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Phit_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phit_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Phit_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Phit_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Phis_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Phis_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Phis_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phis_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Phis_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Phis_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Phis_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Phis_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Phis_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_Rhop_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_Rhop_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_Rhop_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Rhop_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_Rhop_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_Rhop_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_Rhop_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_Rhop_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_Rhop_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_IDA_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_IDA_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_IDA_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_IDA_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_IDA_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_IDA_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_IDAp_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_IDA_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_IDA_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.perturb_Pearson_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Pearson_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_NULL_res_BRCA[[i]]))))
Scan.perturb_Pearson_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Pearson_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Pearson_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Pearson_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Pearson_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Pearson_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Spearman_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Spearman_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_NULL_res_BRCA[[i]]))))
Scan.perturb_Spearman_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Spearman_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Spearman_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Spearman_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Spearman_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Spearman_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Kendall_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Kendall_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_NULL_res_BRCA[[i]]))))
Scan.perturb_Kendall_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Kendall_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Kendall_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Kendall_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Kendall_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Kendall_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Dcor_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Dcor_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_NULL_res_BRCA[[i]]))))
Scan.perturb_Dcor_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Dcor_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Dcor_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dcor_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Dcor_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dcor_ENCORI_res_BRCA[[i]]))))

Scan.perturb_RDC_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_RDC_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_NULL_res_BRCA[[i]]))))
Scan.perturb_RDC_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_RDC_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_TargetScan_res_BRCA[[i]]))))
Scan.perturb_RDC_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_RDC_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_RDC_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_RDC_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Hoeffd_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Hoeffd_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_NULL_res_BRCA[[i]]))))
Scan.perturb_Hoeffd_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Hoeffd_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Hoeffd_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Hoeffd_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Hoeffd_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Hoeffd_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Zscore_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Zscore_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_NULL_res_BRCA[[i]]))))
Scan.perturb_Zscore_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Zscore_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Zscore_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Zscore_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Zscore_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Zscore_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Biweight_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Biweight_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_NULL_res_BRCA[[i]]))))
Scan.perturb_Biweight_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Biweight_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Biweight_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Biweight_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Biweight_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Biweight_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Weighted_rank_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Weighted_rank_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_NULL_res_BRCA[[i]]))))
Scan.perturb_Weighted_rank_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Weighted_rank_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Weighted_rank_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Weighted_rank_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Weighted_rank_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Weighted_rank_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Cosine_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Cosine_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_NULL_res_BRCA[[i]]))))
Scan.perturb_Cosine_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Cosine_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Cosine_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Cosine_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Cosine_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Cosine_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Euclidean_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Euclidean_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_NULL_res_BRCA[[i]]))))
Scan.perturb_Euclidean_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Euclidean_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Euclidean_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Euclidean_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Euclidean_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Euclidean_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Manhattan_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Manhattan_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_NULL_res_BRCA[[i]]))))
Scan.perturb_Manhattan_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Manhattan_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Manhattan_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Manhattan_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Manhattan_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Manhattan_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Canberra_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Canberra_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_NULL_res_BRCA[[i]]))))
Scan.perturb_Canberra_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Canberra_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Canberra_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Canberra_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Canberra_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Canberra_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Chebyshev_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Chebyshev_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_NULL_res_BRCA[[i]]))))
Scan.perturb_Chebyshev_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Chebyshev_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Chebyshev_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Chebyshev_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Chebyshev_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Chebyshev_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Dice_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Dice_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_NULL_res_BRCA[[i]]))))
Scan.perturb_Dice_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Dice_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Dice_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Dice_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Dice_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Dice_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Jaccard_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Jaccard_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_NULL_res_BRCA[[i]]))))
Scan.perturb_Jaccard_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Jaccard_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Jaccard_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Jaccard_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Jaccard_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Jaccard_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Mahalanobis_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Mahalanobis_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_NULL_res_BRCA[[i]]))))
Scan.perturb_Mahalanobis_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Mahalanobis_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Mahalanobis_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Mahalanobis_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Mahalanobis_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Mahalanobis_ENCORI_res_BRCA[[i]]))))

Scan.perturb_MI_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MI_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_MI_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_NULL_res_BRCA[[i]]))))
Scan.perturb_MI_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MI_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_MI_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_TargetScan_res_BRCA[[i]]))))
Scan.perturb_MI_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MI_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_MI_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MI_ENCORI_res_BRCA[[i]]))))

Scan.perturb_MIC_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_MIC_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_NULL_res_BRCA[[i]]))))
Scan.perturb_MIC_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_MIC_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_TargetScan_res_BRCA[[i]]))))
Scan.perturb_MIC_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_MIC_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_MIC_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_MIC_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Lasso_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Lasso_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_NULL_res_BRCA[[i]]))))
Scan.perturb_Lasso_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Lasso_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Lasso_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Lasso_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Lasso_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Lasso_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Elastic_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Elastic_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_NULL_res_BRCA[[i]]))))
Scan.perturb_Elastic_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Elastic_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Elastic_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Elastic_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Elastic_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Elastic_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Ridge_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Ridge_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_NULL_res_BRCA[[i]]))))
Scan.perturb_Ridge_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Ridge_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Ridge_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Ridge_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Ridge_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Ridge_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Phit_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Phit_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_NULL_res_BRCA[[i]]))))
Scan.perturb_Phit_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Phit_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Phit_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phit_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Phit_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phit_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Phis_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Phis_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_NULL_res_BRCA[[i]]))))
Scan.perturb_Phis_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Phis_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Phis_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Phis_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Phis_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Phis_ENCORI_res_BRCA[[i]]))))

Scan.perturb_Rhop_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_Rhop_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_NULL_res_BRCA[[i]]))))
Scan.perturb_Rhop_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_Rhop_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_TargetScan_res_BRCA[[i]]))))
Scan.perturb_Rhop_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_Rhop_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_Rhop_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_Rhop_ENCORI_res_BRCA[[i]]))))

Scan.perturb_IDA_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_IDA_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_NULL_res_BRCA[[i]]))))
Scan.perturb_IDA_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_IDA_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_TargetScan_res_BRCA[[i]]))))
Scan.perturb_IDA_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_IDA_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_IDA_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_IDA_ENCORI_res_BRCA[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.perturb_Pearson_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Pearson_NULL_res_BRCA)
Scan.perturb_Pearson_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Pearson_TargetScan_res_BRCA)
Scan.perturb_Pearson_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Pearson_ENCORI_res_BRCA)

Scan.perturb_Spearman_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Spearman_NULL_res_BRCA)
Scan.perturb_Spearman_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Spearman_TargetScan_res_BRCA)
Scan.perturb_Spearman_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Spearman_ENCORI_res_BRCA)

Scan.perturb_Kendall_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Kendall_NULL_res_BRCA)
Scan.perturb_Kendall_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Kendall_TargetScan_res_BRCA)
Scan.perturb_Kendall_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Kendall_ENCORI_res_BRCA)

Scan.perturb_Dcor_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dcor_NULL_res_BRCA)
Scan.perturb_Dcor_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dcor_TargetScan_res_BRCA)
Scan.perturb_Dcor_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dcor_ENCORI_res_BRCA)

Scan.perturb_RDC_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_RDC_NULL_res_BRCA)
Scan.perturb_RDC_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_RDC_TargetScan_res_BRCA)
Scan.perturb_RDC_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_RDC_ENCORI_res_BRCA)

Scan.perturb_Hoeffd_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Hoeffd_NULL_res_BRCA)
Scan.perturb_Hoeffd_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Hoeffd_TargetScan_res_BRCA)
Scan.perturb_Hoeffd_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Hoeffd_ENCORI_res_BRCA)

Scan.perturb_Zscore_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Zscore_NULL_res_BRCA)
Scan.perturb_Zscore_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Zscore_TargetScan_res_BRCA)
Scan.perturb_Zscore_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Zscore_ENCORI_res_BRCA)

Scan.perturb_Biweight_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Biweight_NULL_res_BRCA)
Scan.perturb_Biweight_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Biweight_TargetScan_res_BRCA)
Scan.perturb_Biweight_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Biweight_ENCORI_res_BRCA)

Scan.perturb_Weighted_rank_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Weighted_rank_NULL_res_BRCA)
Scan.perturb_Weighted_rank_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Weighted_rank_TargetScan_res_BRCA)
Scan.perturb_Weighted_rank_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Weighted_rank_ENCORI_res_BRCA)

Scan.perturb_Cosine_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Cosine_NULL_res_BRCA)
Scan.perturb_Cosine_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Cosine_TargetScan_res_BRCA)
Scan.perturb_Cosine_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Cosine_ENCORI_res_BRCA)

Scan.perturb_Euclidean_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Euclidean_NULL_res_BRCA)
Scan.perturb_Euclidean_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Euclidean_TargetScan_res_BRCA)
Scan.perturb_Euclidean_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Euclidean_ENCORI_res_BRCA)

Scan.perturb_Manhattan_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Manhattan_NULL_res_BRCA)
Scan.perturb_Manhattan_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Manhattan_TargetScan_res_BRCA)
Scan.perturb_Manhattan_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Manhattan_ENCORI_res_BRCA)

Scan.perturb_Canberra_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Canberra_NULL_res_BRCA)
Scan.perturb_Canberra_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Canberra_TargetScan_res_BRCA)
Scan.perturb_Canberra_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Canberra_ENCORI_res_BRCA)

Scan.perturb_Chebyshev_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Chebyshev_NULL_res_BRCA)
Scan.perturb_Chebyshev_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Chebyshev_TargetScan_res_BRCA)
Scan.perturb_Chebyshev_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Chebyshev_ENCORI_res_BRCA)

Scan.perturb_Dice_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dice_NULL_res_BRCA)
Scan.perturb_Dice_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dice_TargetScan_res_BRCA)
Scan.perturb_Dice_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Dice_ENCORI_res_BRCA)

Scan.perturb_Jaccard_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Jaccard_NULL_res_BRCA)
Scan.perturb_Jaccard_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Jaccard_TargetScan_res_BRCA)
Scan.perturb_Jaccard_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Jaccard_ENCORI_res_BRCA)

Scan.perturb_Mahalanobis_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Mahalanobis_NULL_res_BRCA)
Scan.perturb_Mahalanobis_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Mahalanobis_TargetScan_res_BRCA)
Scan.perturb_Mahalanobis_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Mahalanobis_ENCORI_res_BRCA)

Scan.perturb_MI_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MI_NULL_res_BRCA)
Scan.perturb_MI_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MI_TargetScan_res_BRCA)
Scan.perturb_MI_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MI_ENCORI_res_BRCA)

Scan.perturb_MIC_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MIC_NULL_res_BRCA)
Scan.perturb_MIC_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MIC_TargetScan_res_BRCA)
Scan.perturb_MIC_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_MIC_ENCORI_res_BRCA)

Scan.perturb_Lasso_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Lasso_NULL_res_BRCA)
Scan.perturb_Lasso_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Lasso_TargetScan_res_BRCA)
Scan.perturb_Lasso_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Lasso_ENCORI_res_BRCA)

Scan.perturb_Elastic_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Elastic_NULL_res_BRCA)
Scan.perturb_Elastic_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Elastic_TargetScan_res_BRCA)
Scan.perturb_Elastic_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Elastic_ENCORI_res_BRCA)

Scan.perturb_Ridge_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Ridge_NULL_res_BRCA)
Scan.perturb_Ridge_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Ridge_TargetScan_res_BRCA)
Scan.perturb_Ridge_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Ridge_ENCORI_res_BRCA)

Scan.perturb_Phit_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phit_NULL_res_BRCA)
Scan.perturb_Phit_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phit_TargetScan_res_BRCA)
Scan.perturb_Phit_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phit_ENCORI_res_BRCA)

Scan.perturb_Phis_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phis_NULL_res_BRCA)
Scan.perturb_Phis_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phis_TargetScan_res_BRCA)
Scan.perturb_Phis_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Phis_ENCORI_res_BRCA)

Scan.perturb_Rhop_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Rhop_NULL_res_BRCA)
Scan.perturb_Rhop_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Rhop_TargetScan_res_BRCA)
Scan.perturb_Rhop_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_Rhop_ENCORI_res_BRCA)

Scan.perturb_IDA_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_IDA_NULL_res_BRCA)
Scan.perturb_IDA_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_IDA_TargetScan_res_BRCA)
Scan.perturb_IDA_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_IDA_ENCORI_res_BRCA)

save.image("Scan.perturb_BRCA.RData")
