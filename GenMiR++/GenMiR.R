######################################################################################################################################################################################
######################################################################################################################################################################################
############################## Network inference method is GenMiR++, and two strategies (linear interpolation and statistical perturbation) are used #################################
######################################################################################################################################################################################
######################################################################################################################################################################################

################################################################################# ################## #################################################################################
################################################################################# Application in K562 dataset ########################################################################
################################################################################# ################## #################################################################################
library(plyr)
library(igraph)
library(pracma)

ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

list.files <- paste("Scan.interp_K562_NULL", 1:19, ".csv", sep = "")
Scan.interp_GenMiR_NULL_res <- list()
for (i in seq(19)){
    interin <- read.csv(list.files[i], header = FALSE, sep = ",")
    colnames(interin) <- colnames(miRNA_scRNA_norm_filter)
    rownames(interin) <- colnames(mRNA_scRNA_norm_filter)
    Scan.interp_GenMiR_NULL_res[[i]] <- graph_from_incidence_matrix(t(interin))
}

Scan.interp_GenMiR_ENCORI_res <- lapply(seq(Scan.interp_GenMiR_NULL_res), function(i) Scan.interp_GenMiR_NULL_res[[i]] %s% ENCORI_graph)
Scan.interp_GenMiR_TargetScan_res <- lapply(seq(Scan.interp_GenMiR_NULL_res), function(i) Scan.interp_GenMiR_NULL_res[[i]] %s% TargetScan_graph)

list.files <- paste("Scan.perturb_K562_NULL", 1:19, ".csv", sep = "")
Scan.perturb_GenMiR_NULL_res <- list()
for (i in seq(19)){
    interin <- read.csv(list.files[i], header = FALSE, sep = ",")
    colnames(interin) <- colnames(miRNA_scRNA_norm_filter)
    rownames(interin) <- colnames(mRNA_scRNA_norm_filter)
    Scan.perturb_GenMiR_NULL_res[[i]] <- graph_from_incidence_matrix(t(interin))
}

Scan.perturb_GenMiR_ENCORI_res <- lapply(seq(Scan.perturb_GenMiR_NULL_res), function(i) Scan.perturb_GenMiR_NULL_res[[i]] %s% ENCORI_graph)
Scan.perturb_GenMiR_TargetScan_res <- lapply(seq(Scan.perturb_GenMiR_NULL_res), function(i) Scan.perturb_GenMiR_NULL_res[[i]] %s% TargetScan_graph)

# Number of predicted sample-specific miRNA-mRNA interactions
Scan_GenMiR_NULL_res_num <- unlist(lapply(seq(Scan_GenMiR_NULL_res), function(i) nrow(as_data_frame(Scan_GenMiR_NULL_res[[i]] ))))
Scan_GenMiR_TargetScan_res_num <- unlist(lapply(seq(Scan_GenMiR_TargetScan_res), function(i) nrow(as_data_frame(Scan_GenMiR_TargetScan_res[[i]] ))))
Scan_GenMiR_ENCORI_res_num <- unlist(lapply(seq(Scan_GenMiR_ENCORI_res), function(i) nrow(as_data_frame(Scan_GenMiR_ENCORI_res[[i]] ))))

Scan.perturb_GenMiR_NULL_res_num <- unlist(lapply(seq(Scan.perturb_GenMiR_NULL_res), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_NULL_res[[i]] ))))
Scan.perturb_GenMiR_TargetScan_res_num <- unlist(lapply(seq(Scan.perturb_GenMiR_TargetScan_res), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_TargetScan_res[[i]] ))))
Scan.perturb_GenMiR_ENCORI_res_num <- unlist(lapply(seq(Scan.perturb_GenMiR_ENCORI_res), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_ENCORI_res[[i]] ))))

# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.interp_GenMiR_NULL_res_validated <- lapply(seq(Scan.interp_GenMiR_NULL_res), function(i) as_data_frame(Scan.interp_GenMiR_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_GenMiR_TargetScan_res_validated <- lapply(seq(Scan.interp_GenMiR_TargetScan_res), function(i) as_data_frame(Scan.interp_GenMiR_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_GenMiR_ENCORI_res_validated <- lapply(seq(Scan.interp_GenMiR_ENCORI_res), function(i) as_data_frame(Scan.interp_GenMiR_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_GenMiR_NULL_res_validated <- lapply(seq(Scan.perturb_GenMiR_NULL_res), function(i) as_data_frame(Scan.perturb_GenMiR_NULL_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_GenMiR_TargetScan_res_validated <- lapply(seq(Scan.perturb_GenMiR_TargetScan_res), function(i) as_data_frame(Scan.perturb_GenMiR_TargetScan_res[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_GenMiR_ENCORI_res_validated <- lapply(seq(Scan.perturb_GenMiR_ENCORI_res), function(i) as_data_frame(Scan.perturb_GenMiR_ENCORI_res[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.interp_GenMiR_NULL_res_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_NULL_res), function(i) 100*nrow(Scan.interp_GenMiR_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_NULL_res[[i]]))))
Scan.interp_GenMiR_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_TargetScan_res), function(i) 100*nrow(Scan.interp_GenMiR_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_TargetScan_res[[i]]))))
Scan.interp_GenMiR_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_ENCORI_res), function(i) 100*nrow(Scan.interp_GenMiR_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_ENCORI_res[[i]]))))

Scan.perturb_GenMiR_NULL_res_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_NULL_res), function(i) 100*nrow(Scan.perturb_GenMiR_NULL_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_NULL_res[[i]]))))
Scan.perturb_GenMiR_TargetScan_res_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_TargetScan_res), function(i) 100*nrow(Scan.perturb_GenMiR_TargetScan_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_TargetScan_res[[i]]))))
Scan.perturb_GenMiR_ENCORI_res_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_ENCORI_res), function(i) 100*nrow(Scan.perturb_GenMiR_ENCORI_res_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_ENCORI_res[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.interp_GenMiR_NULL_res_Sim <- Sim.network(Scan.interp_GenMiR_NULL_res)
Scan.interp_GenMiR_TargetScan_res_Sim <- Sim.network(Scan.interp_GenMiR_TargetScan_res)
Scan.interp_GenMiR_ENCORI_res_Sim <- Sim.network(Scan.interp_GenMiR_ENCORI_res)

Scan.perturb_GenMiR_NULL_res_Sim <- Sim.network(Scan.perturb_GenMiR_NULL_res)
Scan.perturb_GenMiR_TargetScan_res_Sim <- Sim.network(Scan.perturb_GenMiR_TargetScan_res)
Scan.perturb_GenMiR_ENCORI_res_Sim <- Sim.network(Scan.perturb_GenMiR_ENCORI_res)

save.image("GenMiR.RData")

################################################################################# ######################## ############################################################################
################################################################################# Application in BRCA dataset #########################################################################
################################################################################# ######################## ############################################################################
library(plyr)
library(igraph)
library(pracma)

ENCORI <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
TargetScan <- read.csv("TargetScan_8.0.csv", header = TRUE, sep = ",")
ENCORI_graph <-make_graph(c(t(ENCORI)), directed = FALSE)
TargetScan_graph <-make_graph(c(t(TargetScan)), directed = FALSE)

list.files <- paste("Scan.interp_BRCA_NULL", 1:690, ".csv", sep = "")
Scan.interp_GenMiR_NULL_res_BRCA <- list()
for (i in seq(690)){
    interin <- read.csv(list.files[i], header = FALSE, sep = ",")
    colnames(interin) <- colnames(BRCA_miRNA_Exp_DEG[[2]])
    rownames(interin) <- colnames(BRCA_mRNA_Exp_DEG[[2]])
    Scan.interp_GenMiR_NULL_res_BRCA[[i]] <- graph_from_incidence_matrix(t(interin))
}

Scan.interp_GenMiR_ENCORI_res_BRCA <- lapply(seq(Scan.interp_GenMiR_NULL_res_BRCA), function(i) Scan.interp_GenMiR_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.interp_GenMiR_TargetScan_res_BRCA <- lapply(seq(Scan.interp_GenMiR_NULL_res_BRCA), function(i) Scan.interp_GenMiR_NULL_res_BRCA[[i]] %s% TargetScan_graph)

list.files <- paste("Scan.perturb_BRCA_NULL", 1:690, ".csv", sep = "")
Scan.perturb_GenMiR_NULL_res_BRCA <- list()
for (i in seq(690)){
    interin <- read.csv(list.files[i], header = FALSE, sep = ",")
    colnames(interin) <- colnames(BRCA_miRNA_Exp_DEG[[2]])
    rownames(interin) <- colnames(BRCA_mRNA_Exp_DEG[[2]])
    Scan.perturb_GenMiR_NULL_res_BRCA[[i]] <- graph_from_incidence_matrix(t(interin))
}

Scan.perturb_GenMiR_ENCORI_res_BRCA <- lapply(seq(Scan.perturb_GenMiR_NULL_res_BRCA), function(i) Scan.perturb_GenMiR_NULL_res_BRCA[[i]] %s% ENCORI_graph)
Scan.perturb_GenMiR_TargetScan_res_BRCA <- lapply(seq(Scan.perturb_GenMiR_NULL_res_BRCA), function(i) Scan.perturb_GenMiR_NULL_res_BRCA[[i]] %s% TargetScan_graph)

# Number of predicted sample-specific miRNA-mRNA interactions
Scan.interp_GenMiR_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.interp_GenMiR_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_GenMiR_NULL_res_BRCA[[i]] ))))
Scan.interp_GenMiR_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.interp_GenMiR_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_GenMiR_TargetScan_res_BRCA[[i]] ))))
Scan.interp_GenMiR_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.interp_GenMiR_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.interp_GenMiR_ENCORI_res_BRCA[[i]] ))))

Scan.perturb_GenMiR_NULL_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_GenMiR_NULL_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_NULL_res_BRCA[[i]] ))))
Scan.perturb_GenMiR_TargetScan_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_GenMiR_TargetScan_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_TargetScan_res_BRCA[[i]] ))))
Scan.perturb_GenMiR_ENCORI_res_BRCA_num <- unlist(lapply(seq(Scan.perturb_GenMiR_ENCORI_res_BRCA), function(i) nrow(as_data_frame(Scan.perturb_GenMiR_ENCORI_res_BRCA[[i]] ))))

# Experimentally validated sample-specific miRNA-mRNA interactions
miRTarget_groundtruth <- as.matrix(read.csv("miRTarBase_v9.0+TarBase_v8.0.csv", header = TRUE, sep=","))
miRTarget_groundtruth_graph <- make_graph(c(t(miRTarget_groundtruth[, 1:2])), directed = FALSE)

Scan.interp_GenMiR_NULL_res_BRCA_validated <- lapply(seq(Scan.interp_GenMiR_NULL_res_BRCA), function(i) as_data_frame(Scan.interp_GenMiR_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_GenMiR_TargetScan_res_BRCA_validated <- lapply(seq(Scan.interp_GenMiR_TargetScan_res_BRCA), function(i) as_data_frame(Scan.interp_GenMiR_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.interp_GenMiR_ENCORI_res_BRCA_validated <- lapply(seq(Scan.interp_GenMiR_ENCORI_res_BRCA), function(i) as_data_frame(Scan.interp_GenMiR_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

Scan.perturb_GenMiR_NULL_res_BRCA_validated <- lapply(seq(Scan.perturb_GenMiR_NULL_res_BRCA), function(i) as_data_frame(Scan.perturb_GenMiR_NULL_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_GenMiR_TargetScan_res_BRCA_validated <- lapply(seq(Scan.perturb_GenMiR_TargetScan_res_BRCA), function(i) as_data_frame(Scan.perturb_GenMiR_TargetScan_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))
Scan.perturb_GenMiR_ENCORI_res_BRCA_validated <- lapply(seq(Scan.perturb_GenMiR_ENCORI_res_BRCA), function(i) as_data_frame(Scan.perturb_GenMiR_ENCORI_res_BRCA[[i]] %s% miRTarget_groundtruth_graph))

## Percentage of experimentally validated sample-specific miRNA-mRNA interactions
Scan.interp_GenMiR_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_NULL_res_BRCA), function(i) 100*nrow(Scan.interp_GenMiR_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_NULL_res_BRCA[[i]]))))
Scan.interp_GenMiR_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_TargetScan_res_BRCA), function(i) 100*nrow(Scan.interp_GenMiR_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_TargetScan_res_BRCA[[i]]))))
Scan.interp_GenMiR_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.interp_GenMiR_ENCORI_res_BRCA), function(i) 100*nrow(Scan.interp_GenMiR_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.interp_GenMiR_ENCORI_res_BRCA[[i]]))))

Scan.perturb_GenMiR_NULL_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_NULL_res_BRCA), function(i) 100*nrow(Scan.perturb_GenMiR_NULL_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_NULL_res_BRCA[[i]]))))
Scan.perturb_GenMiR_TargetScan_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_TargetScan_res_BRCA), function(i) 100*nrow(Scan.perturb_GenMiR_TargetScan_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_TargetScan_res_BRCA[[i]]))))
Scan.perturb_GenMiR_ENCORI_res_BRCA_validated_per <- unlist(lapply(seq(Scan.perturb_GenMiR_ENCORI_res_BRCA), function(i) 100*nrow(Scan.perturb_GenMiR_ENCORI_res_BRCA_validated[[i]])/nrow(as_data_frame(Scan.perturb_GenMiR_ENCORI_res_BRCA[[i]]))))

## Network similarity of sample-specific miRNA reguatory network
Scan.interp_GenMiR_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_GenMiR_NULL_res)
Scan.interp_GenMiR_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_GenMiR_TargetScan_res)
Scan.interp_GenMiR_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.interp_GenMiR_ENCORI_res)

Scan.perturb_GenMiR_NULL_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_GenMiR_NULL_res)
Scan.perturb_GenMiR_TargetScan_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_GenMiR_TargetScan_res)
Scan.perturb_GenMiR_ENCORI_res_BRCA_Sim <- Sim.network_parallel(Scan.perturb_GenMiR_ENCORI_res)

save.image("GenMiR.RData")
