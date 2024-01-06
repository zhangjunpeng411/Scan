#######################################################################################################################################################################
#######################################################################################################################################################################
###################################################################################### Utility functions ##############################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
## Function for computing the average expression values of duplicate genes
# Exp_scRNA: Gene expression values of miRNAs or mRNAs, rows are samples and columns are miRNAs or mRNAs
# Output: temp is expression data without duplicate genes
Averg_Duplicate <- function(Exp_scRNA){
    
    uniqueNameList <- unique(colnames(Exp_scRNA))
    noOfgenes <- length(uniqueNameList)
    temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
    colnames(temp) <- uniqueNameList
    rownames(temp) <- rownames(Exp_scRNA)
    for(c in 1:noOfgenes){
        GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
    for(r in 1:nrow(temp)) {
        temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
  }
}
    return(temp)
}

## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
matrixzscore <- function(mat){

    mat.mean <- mean(mat[!is.na(mat)])
    mat.sd <- sd(mat[!is.na(mat)])
    mat.zscore <- (mat - mat.mean)/mat.sd

    return(mat.zscore)
}

## Function for querying gene expression data and miRNA-target interactions by combining gene expression data and putative miRNA-target interactions
# ExpData: Gene expression values of miRNAs and mRNAs, rows are samples and columns are miRNAs and mRNAs
# miRTarget: Putative miRNA-target interactions
# Output: ExpDataQuery is the queried ExpData, and miRTargetQuery is the queried miRTarget
querydata <- function(ExpData, miRTarget) {

    ExpDataNames <- colnames(ExpData)
    miRTarget <- as.matrix(miRTarget)
    
    miRTargetQuery <- miRTarget[intersect(which(miRTarget[, 1] %in% ExpDataNames),
        which(miRTarget[, 2] %in% ExpDataNames)), ] 
	
    ExpDataQuery <- ExpData[, union(which(ExpDataNames %in% unique(miRTargetQuery[, 1])),
        which(ExpDataNames %in% unique(miRTargetQuery[, 2])))]

    return(list(ExpDataQuery, miRTargetQuery))
}

## Function of network inference using Pearson 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Pearson <- function(miRExp, mRExp){    
        
    adj.Matrix <- corAndPvalue(mRExp, miRExp, method = "pearson")$cor
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Pearson.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Pearson_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Pearson(miRExp, mRExp)
    original_single <- Pearson(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Pearson.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Pearson.pvalue)
}

## Function of network inference using Spearman 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Spearman <- function(miRExp, mRExp){    
       
    adj.Matrix <- corAndPvalue(mRExp, miRExp, method = "spearman")$cor
     
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Spearman.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Spearman_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Spearman(miRExp, mRExp)
    original_single <- Spearman(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Spearman.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Spearman.pvalue)
}

## Function of network inference using Kendall 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Kendall <- function(miRExp, mRExp){    
        
    adj.Matrix <- corAndPvalue(mRExp, miRExp, method = "kendall")$cor
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Kendall.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Kendall_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Kendall(miRExp, mRExp)
    original_single <- Kendall(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Kendall.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Kendall.pvalue)
}

## Function of network inference using Dcor 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Dcor <- function(miRExp, mRExp){
    
      adj.Matrix <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            adj.Matrix[j, i] <- dcorT(mRExp[, j], miRExp[, i])
         }
      }
      colnames(adj.Matrix) <- colnames(miRExp)
      rownames(adj.Matrix) <- colnames(mRExp)
    
      return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Dcor.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Dcor_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Dcor(miRExp, mRExp)
    original_single <- Dcor(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Dcor.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Dcor.pvalue)
}


# Random Dependence Coefficient (RDC) calculate the dependence between the random samples as the highest canonical
# correlation between the k random non-linear projections of their copula transformations.
# 
# It also provides the p-value for the independence hypothesis, assuming the
# normality of the data through the Bartlett's approximation.
# x: A vector, matrix or numeric data frame
# y: A vector, matrix or numeric data frame
# k: Number of non-linear projections of the copula, by default k = 20
# s: Variance to draw i.i.d. projection coefficients in N ~(0, sI),
# by defect is 1/6
# f: Function that is used for the generation of random non-linear
# projections, if it is not indicated it uses the sinusoidal projections (sin)
# Output: rdc is random dependence coefficient between x and y
rdc_generic <- function(x, y, k = 20, s = 1/6, f = sin) {

  x <- cbind(apply(as.matrix(x), 2, function(u)rank(u)/length(u)), 1)
  y <- cbind(apply(as.matrix(y), 2, function(u)rank(u)/length(u)), 1)
  x <- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k), ncol(x))
  y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k), ncol(y))
  # can <- cancor(cbind(f(x), 1), cbind(f(y), 1))$cor
  try_can <- try(cancor(cbind(f(x), 1), cbind(f(y), 1))$cor, silent=TRUE)
        if ("try-error" %in% class(try_can)) {
        can <- 0
        } else {
        can <-  try_can
	}

  k <- length(can)
  chi <- (((2*k+3)/2)-nrow(x))*log(prod(1-can^2))
  rdc <- as.data.frame(cbind(can[1], pchisq(chi,k^2, lower.tail = FALSE)))
  colnames(rdc) <- c("rdc", "p.value")
  return(rdc)
}

## Function of network inference using RDC 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
RDC <- function(miRExp, mRExp){
    
      adj.Matrix <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            adj.Matrix[j, i] <- rdc_generic(mRExp[, j], miRExp[, i])$rdc
         }
      }
      colnames(adj.Matrix) <- colnames(miRExp)
      rownames(adj.Matrix) <- colnames(mRExp)
    
     return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: RDC.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
RDC_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- RDC(miRExp, mRExp)
    original_single <- RDC(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    RDC.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(RDC.pvalue)
}

## Function of network inference using Hoeffd 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Hoeffd <- function(miRExp, mRExp){
       
    adj.Matrix <- hoeffd(cbind(miRExp, mRExp))$D
    adj.Matrix <- adj.Matrix[which(rownames(adj.Matrix) %in% colnames(mRExp)), which(colnames(adj.Matrix) %in% colnames(miRExp))]
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Hoeffd.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Hoeffd_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Hoeffd(miRExp, mRExp)
    original_single <- Hoeffd(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Hoeffd.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Hoeffd.pvalue)
}

## Function of network inference using Zscore 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Zscore <- function(miRExp, mRExp){
    
      adj.Matrix <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            interin <- mRExp[which(miRExp[, i]==min(miRExp[, i])), j]
            interin <- median(interin)
	    sd.value <- sd(mRExp[, j])
            mean.value <- mean(mRExp[, j])
	    adj.Matrix[j, i] <- (interin - mean.value)/sd.value
         }
      }
      colnames(adj.Matrix) <- colnames(miRExp)
      rownames(adj.Matrix) <- colnames(mRExp)
   
      return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Zscore.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Zscore_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Zscore(miRExp, mRExp)
    original_single <- Zscore(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Zscore.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Zscore.pvalue)
}

## Function of network inference using MI 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
MI <- function(miRExp, mRExp){
    
    adj.Matrix <- knnmi.cross(t(mRExp), t(miRExp), k = 3, noise = 1e-12)
   
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: MI.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
MI_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- MI(miRExp, mRExp)
    original_single <- MI(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    MI.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(MI.pvalue)
}

## Function of network inference using MIC 
# Maximal information coefficient (MIC) calculates the dependence of pairs of variables. The MIC belongs to a larger
# classification of nonparametric exploration statistics based on maximum information (MINE).
#
# It also provides the p-value for the independence hypothesis, that are
# generated by permuting the data and seeing how likely it is that the observed
# MIC value arises from the perturbed data.
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
MIC <- function(miRExp, mRExp){
     
      adj.Matrix <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            adj.Matrix[j, i] <- mine(mRExp[, j], miRExp[, i])$MIC
         }
      }
      colnames(adj.Matrix) <- colnames(miRExp)
      rownames(adj.Matrix) <- colnames(mRExp)
    
      return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: MIC.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
MIC_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- MIC(miRExp, mRExp)
    original_single <- MIC(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    MIC.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(MIC.pvalue)
}

## Function for calculating regression between x and y
# x: variable 1
# y: variable 2
# Output: lasso.cor is the regression coefficient
lasso_generic <- function(x, y) {

    lambda.sel <- cv.glmnet(x, y, family = "mgaussian", alpha = 1)$lambda.min
    fit <- glmnet(x, y, alpha = 1, family = "mgaussian")
    lasso.cor <- as.matrix(do.call(cbind, coef(fit, s = lambda.sel))[-1, ])    
    colnames(lasso.cor) <- colnames(y)
    rownames(lasso.cor) <- colnames(x)
    return(lasso.cor)
}

## Function of network inference using Lasso 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Lasso <- function(miRExp, mRExp){
        
    adj.Matrix <- lasso_generic(mRExp, miRExp)       
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Lasso.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Lasso_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Lasso(miRExp, mRExp)
    original_single <- Lasso(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Lasso.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Lasso.pvalue)
}

## Function for calculating regression between x and y
# x: variable 1
# y: variable 2
# Output: elastic.cor is the regression coefficient
elastic_generic <- function(x, y) {

    lambda.sel <- cv.glmnet(x, y, family = "mgaussian", alpha = 0.5)$lambda.min
    fit <- glmnet(x, y, alpha = 0.5, family = "mgaussian")
    elastic.cor <- as.matrix(do.call(cbind, coef(fit, s = lambda.sel))[-1, ])
    colnames(elastic.cor) <- colnames(y)
    rownames(elastic.cor) <- colnames(x)
    return(elastic.cor)
}

## Function of network inference using Elastic 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Elastic <- function(miRExp, mRExp){
      
    adj.Matrix <- elastic_generic(mRExp, miRExp)       
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Elastic.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Elastic_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Elastic(miRExp, mRExp)
    original_single <- Elastic(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Elastic.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Elastic.pvalue)
}

## Function for calculating regression between x and y
# x: variable 1
# y: variable 2
# Output: ridge.cor is the regression coefficient
ridge_generic <- function(x, y) {

    lambda.sel <- cv.glmnet(x, y, family = "mgaussian", alpha = 0)$lambda.min
    fit <- glmnet(x, y, alpha = 0, family = "mgaussian")
    ridge.cor <- as.matrix(do.call(cbind, coef(fit, s = lambda.sel))[-1, ])
    colnames(ridge.cor) <- colnames(y)
    rownames(ridge.cor) <- colnames(x)
    return(ridge.cor)
}

## Function of network inference using Ridge 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Ridge <- function(miRExp, mRExp){
        
    adj.Matrix <- ridge_generic(mRExp, miRExp)       
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Ridge.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Ridge_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Ridge(miRExp, mRExp)
    original_single <- Ridge(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Ridge.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Ridge.pvalue)
}

## Function of network inference using IDA 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# pcmethod: Character string specifying the method of the PC algorithm, including "stable", "original", "stable.fast".
# alpha: Significance level (number in [0, 1]) for the individual conditional independence tests.
# Output: ida.adj is a correlation matrix between miRNAs and mRNAs
IDA <- function(miRExp, mRExp,                  
                pcmethod = c("stable", "original", "stable.fast"), 
		alpha = 0.01){    
         
    Exp <- scale(cbind(miRExp, mRExp))
    miRmR.comb <- expand.grid(colnames(miRExp), colnames(mRExp))
    colnames(miRmR.comb) <- c("mir", "gene") 
    miRNA.index <- unlist(lapply(seq(nrow(miRmR.comb)), function(i) which(colnames(Exp) %in% miRmR.comb[i, 1])))
    mRNA.index <- unlist(lapply(seq(nrow(miRmR.comb)), function(i) which(colnames(Exp) %in% miRmR.comb[i, 2])))
    index <- cbind(miRNA.index, mRNA.index, miRmR.comb)
    miRNA.index <- as.numeric(index[order(index[, 1]), ][, 1])
    mRNA.index <- as.numeric(index[order(index[, 1]), ][, 2])
    prior <- index[order(index[, 1]), ][, 3:4]
    num <- table(miRNA.index)
    num.names <- as.numeric(names(num)) 
   
    suffStat = list(C = stats::cor(Exp), n = nrow(Exp))
    cov.d <- cov(Exp)	
        
    pcFit <- pc(suffStat, indepTest = gaussCItest, p = ncol(Exp), alpha = alpha, skel.method = pcmethod)        
    
    interin <- list()
    for (i in seq(num.names)) {
        try_ida <- try(idaFast(num.names[i], mRNA.index[(sum(num[0:(i-1)])+1):sum(num[0:i])], cov.d, pcFit@graph), silent=TRUE)
        if ("try-error" %in% class(try_ida)) {
        interin[[i]] <- as.matrix(rep(0, length((sum(num[0:(i-1)])+1):sum(num[0:i]))))
        } else {
        interin[[i]] <-  try_ida
	}
    }     
    
    for (i in seq(num.names)) {
        interin[[i]][which(is.na(interin[[i]]) == TRUE)] <- 0
    }

    interinabs <- lapply(seq(num.names), function(i) 
	                 abs(interin[[i]]))

    ef <- unlist(lapply(seq(num.names), function(i) 
                       interin[[i]][cbind(seq(nrow(interinabs[[i]])), 
		       unlist(apply(interinabs[[i]], 1, function(x){which.min(x)})))]))     

    ida.res <- data.frame(prior, ef)

    G <- graph.data.frame(ida.res, directed = FALSE)
    adj.matrix <- as_adjacency_matrix(G, type="both", names=TRUE, sparse=FALSE, attr="ef")

    
    miRlist <- unique(prior[, 1])
    Tarlist <- unique(prior[, 2])
    ida.adj <- adj.matrix[which(rownames(adj.matrix) %in% Tarlist), which(colnames(adj.matrix) %in% miRlist)]
   
    return(ida.adj)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: IDA.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
IDA_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- IDA(miRExp, mRExp, pcmethod = "stable")
    original_single <- IDA(miRExp[-sample_index, ], mRExp[-sample_index, ], pcmethod = "stable")
    original <- nsamples * original_all - (nsamples - 1) * original_single
    IDA.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(IDA.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: euclidean.value is the distance
euclidean_generic <- function(x, y) {

    interin <- as.matrix(stats::dist(t(cbind(x, y)), method = 'euclidean'))
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    euclidean.value <- interin[y.index, x.index]
    
    return(euclidean.value)
}

## Function of network inference using Euclidean 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Euclidean <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(euclidean_generic(miRExp, mRExp) + eps(1))
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Euclidean.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Euclidean_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Euclidean(miRExp, mRExp)
    original_single <- Euclidean(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Euclidean.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Euclidean.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: manhattan.value is the distance
manhattan_generic <- function(x, y) {

    interin <- as.matrix(stats::dist(t(cbind(x, y)), method = 'manhattan'))
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    manhattan.value <- interin[y.index, x.index]
     
    return(manhattan.value)
}

## Function of network inference using Manhattan 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Manhattan <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(manhattan_generic(miRExp, mRExp) + eps(1))
   
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Manhattan.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Manhattan_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Manhattan(miRExp, mRExp)
    original_single <- Manhattan(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Manhattan.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Manhattan.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: canberra.value is the distance
canberra_generic <- function(x, y) {

    interin <- as.matrix(stats::dist(t(cbind(x, y)), method = 'canberra'))
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    canberra.value <- interin[y.index, x.index]
     
    return(canberra.value)
}

## Function of network inference using Canberra 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Canberra <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(canberra_generic(miRExp, mRExp) + eps(1))

    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Canberra.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Canberra_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Canberra(miRExp, mRExp)
    original_single <- Canberra(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Canberra.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Canberra.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: chebyshev.value is the distance
chebyshev_generic <- function(x, y) {

    interin <- as.matrix(stats::dist(t(cbind(x, y)), method = 'maximum'))
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    chebyshev.value <- interin[y.index, x.index]
     
    return(chebyshev.value)
}

## Function of network inference using Chebyshev 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Chebyshev <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(chebyshev_generic(miRExp, mRExp) + eps(1))
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Chebyshev.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Chebyshev_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Chebyshev(miRExp, mRExp)
    original_single <- Chebyshev(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Chebyshev.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Chebyshev.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: dice.value is the distance
dice_generic <- function(x, y) {
    
    interin <- philentropy::distance(t(cbind(x, y)), method = "dice")     
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    dice.value <- interin[y.index, x.index]
    colnames(dice.value) <- colnames(x)
    rownames(dice.value) <- colnames(y)

    return(dice.value)
}

## Function of network inference using Dice 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Dice <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(dice_generic(miRExp, mRExp) + eps(1))
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Dice.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Dice_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Dice(miRExp, mRExp)
    original_single <- Dice(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Dice.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Dice.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: jaccard.value is the distance
jaccard_generic <- function(x, y) {
    
    interin <- philentropy::distance(t(cbind(x, y)), method = "jaccard")     
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    jaccard.value <- interin[y.index, x.index]
    colnames(jaccard.value) <- colnames(x)
    rownames(jaccard.value) <- colnames(y)

    return(jaccard.value)
}

## Function of network inference using Jaccard 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Jaccard <- function(miRExp, mRExp){
   
    adj.Matrix <- 1/(jaccard_generic(miRExp, mRExp) + eps(1))
     
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Jaccard.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Jaccard_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Jaccard(miRExp, mRExp)
    original_single <- Jaccard(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Jaccard.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Jaccard.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: mahalanobis.value is the distance
mahalanobis_generic <- function(x, y) {
    
    mahalanobis.value <- StatMatch::mahalanobis.dist(t(y), t(x))
     
    return(mahalanobis.value)
}

## Function of network inference using Mahalanobis 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Mahalanobis <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(mahalanobis_generic(miRExp, mRExp) + eps(1))
     
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Mahalanobis.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Mahalanobis_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Mahalanobis(miRExp, mRExp)
    original_single <- Mahalanobis(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Mahalanobis.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Mahalanobis.pvalue)
}

## Function for calculating correlation between x and y
# x: variable 1
# y: variable 2
# Output: cosine.value is the correlation
cosine_generic <- function(x, y) {
    
    interin <- philentropy::distance(t(cbind(x, y)), method = "cosine")     
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    cosine.value <- interin[y.index, x.index]
    colnames(cosine.value) <- colnames(x)
    rownames(cosine.value) <- colnames(y)

    return(cosine.value)
}

## Function of network inference using Cosine
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Cosine <- function(miRExp, mRExp){
    
    adj.Matrix <- cosine_generic(miRExp, mRExp)
   
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Cosine.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Cosine_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Cosine(miRExp, mRExp)
    original_single <- Cosine(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Cosine.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Cosine.pvalue)
}

## Function for calculating correlation between x and y
# x: variable 1
# y: variable 2
# Output: biweight.value is the correlation
biweight_generic <- function(x, y){
    
    biweight.value <- WGCNA::bicor(y, x)
    
    return(biweight.value)
}

## Function of network inference using Biweight
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Biweight <- function(miRExp, mRExp){
    
    adj.Matrix <- biweight_generic(miRExp, mRExp)
   
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Biweight.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Biweight_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Biweight(miRExp, mRExp)
    original_single <- Biweight(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Biweight.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Biweight.pvalue)
}

## Function for calculating correlation between x and y
# x: variable 1
# y: variable 2
# Output: wtd_rank.value is the correlation
wtd_rank_generic <- function(x, y) {

    mat <- cbind(x, y)
    ranks <- apply(mat, 2, rank, ties = "average")
    # weight the ranks, calculate the savage scores  
    n <- nrow(mat)
    reciprocals <- 1 / seq_len(n)
    savage <- sapply(seq_len(n), function(i) sum(reciprocals[i:n]))
    # replace each rank with the savage score 
    savages <- ranks
    savages[] <- savage[ranks]
    # calculate pearson correlation
    wtd_rank.value <- WGCNA::corAndPvalue(savages, method = "pearson")$cor
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    wtd_rank.value <- wtd_rank.value[y.index, x.index]
    
    return(wtd_rank.value)

}

## Function of network inference using Weighted_rank
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Weighted_rank <- function(miRExp, mRExp){
    
    adj.Matrix <- wtd_rank_generic(miRExp, mRExp)
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Weighted_rank.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Weighted_rank_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Weighted_rank(miRExp, mRExp)
    original_single <- Weighted_rank(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Weighted_rank.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Weighted_rank.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: phit.value is the distance
phit_generic <- function(x, y){
    
    mat <- cbind(x, y)
    interin <- propr::phit(mat, select = colnames(mat))@matrix
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    phit.value <- interin[y.index, x.index]
    
    return(phit.value)
}

## Function of network inference using Phit
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Phit <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(phit_generic(miRExp, mRExp) + eps(1))
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Phit.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Phit_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Phit(miRExp, mRExp)
    original_single <- Phit(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Phit.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Phit.pvalue)
}

## Function for calculating distance between x and y
# x: variable 1
# y: variable 2
# Output: phis.value is the distance
phis_generic <- function(x, y){
    
    mat <- cbind(x, y)
    interin <- propr::phis(mat, select = colnames(mat))@matrix
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    phis.value <- interin[y.index, x.index]
    
    return(phis.value)
}

## Function of network inference using Phit
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Phis <- function(miRExp, mRExp){
    
    adj.Matrix <- 1/(phis_generic(miRExp, mRExp) + eps(1))
    
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Phis.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Phis_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Phis(miRExp, mRExp)
    original_single <- Phis(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Phis.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Phis.pvalue)
}

## Function for calculating correlation between x and y
# x: variable 1
# y: variable 2
# Output: rhop.value is the correlation
rhop_generic <- function(x, y){
    
    mat <- cbind(x, y)
    interin <- propr::perb(mat, select = colnames(mat))@matrix
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    rhop.value <- interin[y.index, x.index]
    
    return(rhop.value)
}

## Function of network inference using Phit
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# Output: adj.Matrix is a correlation matrix between miRNAs and mRNAs
Rhop <- function(miRExp, mRExp){
    
    adj.Matrix <- rhop_generic(miRExp, mRExp)
   
    return(adj.Matrix)
}

## Function for inferring the significant p-values between miRNAs and mRNAs in a specific sample of interest
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# sample_index: the index of specific sample
# Output: Rhop.pvalue is the significant p-values between miRNAs and mRNAs in a specific sample of interest
Rhop_Scan <- function(miRExp, mRExp, sample_index){
    
    nsamples <- nrow(miRExp)
    original_all <- Rhop(miRExp, mRExp)
    original_single <- Rhop(miRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Rhop.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Rhop.pvalue)
}

## Identifying sample-specific miRNA regulation using Scan (Sample-speCific miRNA regulAtioN) and linear interpolation strategy
# miRExp: miRNA expression data with rows are cells and columns are miRNAs.
# mRExp: mRNA expression data with rows are cells and columns are mRNAs.
# method: Methods for calculating correlations, distances or causal effects.
# pcmethod: Character string specifying the method of the PC algorithm, including "stable", "original", "stable.fast".
# alpha: Significance level (number in [0, 1]) for the individual conditional independence tests.
# p.value: Significance level for the identified miRNA-mRNA interactions.
# num.cores: Number of CPU cores when parallel computation.
# Output: An igraph object, graph of miRNA-mRNA regulatory network for each sample
Scan.interp <- function(miRExp, mRExp, 
                    method = c("IDA", "Pearson", "Spearman", "Kendall", 
		               "Dcor", "RDC", "Hoeffd", "Zscore", 
			       "MI", "MIC", "Lasso", "Elastic", "Ridge",
			       "Euclidean", "Manhattan", "Canberra",
			       "Chebyshev", "Dice", "Jaccard",
			       "Mahalanobis", "Cosine", "Biweight", 
			       "Weighted_rank", "Phit", 
			       "Phis", "Rhop"), 
                    pcmethod = c("stable", "original", "stable.fast"), 
		    alpha = 0.01,
		    p.value = 0.05,
		    num.cores = 32){
    
    nsamples <- nrow(miRExp)
    
    if (method == "IDA"){
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("pcalg", "igraph", "WGCNA"), .export = c("IDA_Scan", "IDA", "querydata", "matrixzscore")) %dopar% {
            IDA_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()    
    
    } else if (method == "Pearson"){
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Pearson_Scan", "Pearson", "querydata", "matrixzscore")) %dopar% {
            Pearson_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Spearman"){
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Spearman_Scan", "Spearman", "querydata", "matrixzscore")) %dopar% {
            Spearman_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Kendall"){
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Kendall_Scan", "Kendall", "querydata", "matrixzscore")) %dopar% {
            Kendall_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
     
    } else if (method == "Dcor"){
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "energy", "WGCNA"), .export = c("Dcor_Scan", "Dcor", "querydata", "matrixzscore")) %dopar% {
            Dcor_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "RDC"){
 	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("RDC_Scan", "RDC", "rdc_generic", "querydata", "matrixzscore")) %dopar% {
            RDC_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Hoeffd"){
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "Hmisc", "WGCNA"), .export = c("Hoeffd_Scan", "Hoeffd", "querydata", "matrixzscore")) %dopar% {
            Hoeffd_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()    

    } else if (method == "Zscore"){
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Zscore_Scan", "Zscore", "querydata", "matrixzscore")) %dopar% {
            Zscore_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "MI"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "parmigene", "WGCNA"), .export = c("MI_Scan", "MI", "querydata", "matrixzscore")) %dopar% {
            MI_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "MIC"){        
       
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "minerva", "WGCNA"), .export = c("MIC_Scan", "MIC", "querydata", "matrixzscore")) %dopar% {
            MIC_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Lasso"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet", "WGCNA"), .export = c("Lasso_Scan", "Lasso", "lasso_generic", "querydata", "matrixzscore")) %dopar% {
            Lasso_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Elastic"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet", "WGCNA"), .export = c("Elastic_Scan", "Elastic", "elastic_generic", "querydata", "matrixzscore")) %dopar% {
            Elastic_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Ridge"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet", "WGCNA"), .export = c("Ridge_Scan", "Ridge", "ridge_generic", "querydata", "matrixzscore")) %dopar% {
            Ridge_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Euclidean"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "WGCNA", "pracma"), .export = c("Euclidean_Scan", "Euclidean", "euclidean_generic", "querydata", "matrixzscore")) %dopar% {
            Euclidean_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Manhattan"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "WGCNA", "pracma"), .export = c("Manhattan_Scan", "Manhattan", "manhattan_generic", "querydata", "matrixzscore")) %dopar% {
            Manhattan_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Canberra"){        
       
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "WGCNA", "pracma"), .export = c("Canberra_Scan", "Canberra", "canberra_generic", "querydata", "matrixzscore")) %dopar% {
            Canberra_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Chebyshev"){        
       
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "WGCNA", "pracma"), .export = c("Chebyshev_Scan", "Chebyshev", "chebyshev_generic", "querydata", "matrixzscore")) %dopar% {
            Chebyshev_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
   
    } else if (method == "Dice"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy", "WGCNA", "pracma"), .export = c("Dice_Scan", "Dice", "dice_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Dice_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Jaccard"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy", "WGCNA", "pracma"), .export = c("Jaccard_Scan", "Jaccard", "jaccard_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Jaccard_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
   
    } else if (method == "Mahalanobis"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "StatMatch", "WGCNA", "pracma"), .export = c("Mahalanobis_Scan", "Mahalanobis", "mahalanobis_generic", "mahalanobis.dist", "querydata", "matrixzscore")) %dopar% {
            Mahalanobis_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Cosine"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy", "WGCNA"), .export = c("Cosine_Scan", "Cosine", "cosine_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Cosine_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Biweight"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Biweight_Scan", "Biweight", "biweight_generic", "querydata", "matrixzscore")) %dopar% {
            Biweight_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Weighted_rank"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Weighted_rank_Scan", "Weighted_rank", "wtd_rank_generic", "querydata", "matrixzscore")) %dopar% {
            Weighted_rank_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
    
    } else if (method == "Phit"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr", "WGCNA", "pracma"), .export = c("Phit_Scan", "Phit", "phit_generic", "phit", "querydata", "matrixzscore")) %dopar% {
             Phit_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Phis"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr", "WGCNA", "pracma"), .export = c("Phis_Scan", "Phis", "phis_generic", "phis", "querydata", "matrixzscore")) %dopar% {
             Phis_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Rhop"){        
        
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr", "WGCNA"), .export = c("Rhop_Scan", "Rhop", "rhop_generic", "perb", "querydata", "matrixzscore")) %dopar% {
            Rhop_Scan(miRExp, mRExp, i)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
    } 
    
    int.sample.specific.incident <- lapply(seq(res.single), function(i)
                                         ifelse(res.single[[i]] < p.value, 1, 0))   
   
    int.sample.specific.graph <- lapply(seq(res.single), function(i)
                                      graph_from_incidence_matrix(t(int.sample.specific.incident[[i]])))

    return(int.sample.specific.graph)
}


## Function for calculating similarity matrix between two list of networks (original version)
# net: List object, list of network
# Output: Sim is a similarity matrix between two list of networks
Sim.network <- function(net){
    
    net1 <- net
    net2 <- net
    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_interin <- nrow(as_data_frame(net1[[i]] %s% net2[[j]]))
	    Sim[i, j] <- overlap_interin/min(nrow(as_data_frame(net1[[i]])), nrow(as_data_frame(net2[[j]])))
	}
    }
        
    return(Sim)
}

## Function for calculating similarity matrix between two list of networks (parallel version)
# net: List object, list of network
# Output: Sim is a similarity matrix between two list of networks
Sim.network_parallel <- function(net, num.cores = 32){
    
    net1 <- net
    net2 <- net
    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    num <- m * n
    Sim <- matrix(NA, m, n)
    index <- matrix(NA, nrow = num, ncol = 2)
    
    for (i in seq(m)){
        for (j in seq(n)){
            index[(i-1)*n+j, 1] <- i
            index[(i-1)*n+j, 2] <- j
        }
    }

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl) 
    interin <- foreach(k = seq(num),
                       .packages = c("igraph")) %dopar% {
	               nrow(as_data_frame(net1[[index[k, 1]]] %s% net2[[index[k, 2]]]))/min(nrow(as_data_frame(net1[[index[k, 1]]])), nrow(as_data_frame(net2[[index[k, 2]]])))
		       }
    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()
    interin <- do.call(rbind, interin)
    
    for (i in seq(m)){
        for (j in seq(n)){
            Sim[i, j] <- interin[(i-1)*n+j]           
        }
    }
    
    return(Sim)
}


