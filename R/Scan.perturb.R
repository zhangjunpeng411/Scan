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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Pearson <- function(miRExp, mRExp, p.value = 0.05){    
     
    cor.pvalue <- corAndPvalue(mRExp, miRExp, method = "pearson")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
}

## Function of network inference using Spearman 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Spearman <- function(miRExp, mRExp, p.value = 0.05){
    
    cor.pvalue <- corAndPvalue(mRExp, miRExp, method = "spearman")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
}

## Function of network inference using Kendall 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Kendall <- function(miRExp, mRExp, p.value = 0.05){    
        
    cor.pvalue <- corAndPvalue(mRExp, miRExp, method = "kendall")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
}

## Function of network inference using Dcor 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: dcor.adj is a zero-one matrix between miRNAs and mRNAs
Dcor <- function(miRExp, mRExp, p.value = 0.05){
    
      cor.pvalue <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            cor.pvalue[j, i] <- dcorT.test(mRExp[, j], miRExp[, i])$p.value
         }
      }
      colnames(cor.pvalue) <- colnames(miRExp)
      rownames(cor.pvalue) <- colnames(mRExp)
      dcor.adj <- ifelse(cor.pvalue < p.value, 1, 0)    
    
      return(dcor.adj)
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
# p.value: A cutoff of p-value
# Output: rdc.adj is a zero-one matrix between miRNAs and mRNAs
RDC <- function(miRExp, mRExp, p.value = 0.05){
    
      cor.pvalue <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            cor.pvalue[j, i] <- rdc_generic(mRExp[, j], miRExp[, i])$p.value
         }
      }
      colnames(cor.pvalue) <- colnames(miRExp)
      rownames(cor.pvalue) <- colnames(mRExp)
      rdc.adj <- ifelse(cor.pvalue < p.value, 1, 0)    
    
      return(rdc.adj)
}

## Function of network inference using Hoeffd 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Hoeffd <- function(miRExp, mRExp, p.value = 0.05){
        
    cor.pvalue <- hoeffd(cbind(miRExp, mRExp))$P
    cor.pvalue <- cor.pvalue[which(rownames(cor.pvalue) %in% colnames(mRExp)), which(colnames(cor.pvalue) %in% colnames(miRExp))]
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
}

## Function of network inference using Zscore 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Zscore <- function(miRExp, mRExp, p.value = 0.05){
    
      cor.pvalue <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            interin <- mRExp[which(miRExp[, i]==min(miRExp[, i])), j]
            interin <- median(interin)
	    sd.value <- sd(mRExp[, j])
            mean.value <- mean(mRExp[, j])
	    cor <- (interin - mean.value)/sd.value
	    cor.pvalue[j, i] <- 1 - pnorm(abs(cor))
         }
      }
      colnames(cor.pvalue) <- colnames(miRExp)
      rownames(cor.pvalue) <- colnames(mRExp)
      adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
   
     return(adj.Matrix)
}

## Function for inferring mutual information between miRNAs and mRNAs
# mat1: the first gene expression 
# mat2: the second gene expression
# k and noise: Default parameters in MI method
# Output: original is the mutual information between miRNAs and mRNAs
mi_generic <- function(mat1, mat2, k = 3, noise = 1e-12){
    
    original <- knnmi.cross(t(mat1), t(mat2), k = k, noise = noise)
    return(original)
}

## Function of network inference using MI 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
MI <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- mi_generic(mRExp, miRExp)
    cor.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
}

## Function for calculating maximal information coefficient between x and y
# Maximal information coefficient (MIC) calculates the dependence of pairs of variables. The MIC belongs to a larger
# classification of nonparametric exploration statistics based on maximum information (MINE).
# x: variable 1
# y: variable 2
# Output: original is the maximal information coefficient
mic_generic <- function(x, y) {

    original <- mine(x = x, y = y)$MIC      
    return(original)
}

## Function of network inference using MIC 
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: mic.adj is a zero-one matrix between miRNAs and mRNAs
MIC <- function(miRExp, mRExp, p.value = 0.05){
        
      original <- matrix(nrow = ncol(mRExp), ncol = ncol(miRExp))
      for (i in seq(ncol(miRExp))){
        for (j in seq(ncol(mRExp))){
            original[j, i] <- mic_generic(mRExp[, j], miRExp[, i])
         }
      }
      mic.pvalue <- 1 - pnorm(matrixzscore(original))    
      colnames(mic.pvalue) <- colnames(miRExp)
      rownames(mic.pvalue) <- colnames(mRExp)
      mic.adj <- ifelse(mic.pvalue < p.value, 1, 0) 
    
      return(mic.adj)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Lasso <- function(miRExp, mRExp, p.value = 0.05){
        
    lasso.cor <- lasso_generic(mRExp, miRExp)
    lasso.pvalue <- 1 - pnorm(abs(matrixzscore(lasso.cor)))
    adj.Matrix <- ifelse(lasso.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Elastic <- function(miRExp, mRExp, p.value = 0.05){
       
      elastic.cor <- elastic_generic(mRExp, miRExp)
      elastic.pvalue <- 1 - pnorm(abs(matrixzscore(elastic.cor)))
      adj.Matrix <- ifelse(elastic.pvalue < p.value, 1, 0)    
    
      return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Ridge <- function(miRExp, mRExp, p.value = 0.05){
        
      ridge.cor <- ridge_generic(mRExp, miRExp)
      ridge.pvalue <- 1 - pnorm(abs(matrixzscore(ridge.cor)))
      adj.Matrix <- ifelse(ridge.pvalue < p.value, 1, 0)    
   
      return(adj.Matrix)
}

## Function of network inference using IDA
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# pcmethod: Character string specifying the method of the PC algorithm, including "stable", "original", "stable.fast".
# alpha: Significance level (number in [0, 1]) for the individual conditional independence tests.
# p.value: A cutoff of p-value
# Output: ida.adj is a zero-one matrix between miRNAs and mRNAs
IDA <- function(miRExp, mRExp,                  
                pcmethod = c("stable", "original", "stable.fast"), 
		alpha = 0.01, p.value = 0.05){    
       
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
    
    ef <- c()
    for (i in seq(num.names)) {
        try_ida <- try(idaFast(num.names[i], mRNA.index[(sum(num[0:(i-1)])+1):sum(num[0:i])], cov.d, pcFit@graph), silent=TRUE)
        if ("try-error" %in% class(try_ida)) {
        interin <- as.matrix(rep(0, length((sum(num[0:(i-1)])+1):sum(num[0:i]))))
        } else {
        interin <-  try_ida
	}

        interin[which(is.na(interin) == TRUE)] <- 0
        interinabs <- abs(interin)
        ef <- c(ef, interin[cbind(seq(nrow(interinabs)), 
		apply(interinabs, 1, function(x){which.min(x)}))])
    }     
     
    ef[which(ef>=1, arr.ind = TRUE)] <- 1
    ef[which(ef<=-1, arr.ind = TRUE)] <- -1
    ef.pvalue <- corPvalueFisher(ef, nrow(Exp))

    res <- data.frame(prior, ef, ef.pvalue)    
    int <- ifelse(res$ef.pvalue < p.value, 1, 0)
    ida.res <- data.frame(res, int)
    G <- graph.data.frame(ida.res[, c("mir", "gene", "int")], directed = FALSE)
    adj.matrix <- as_adjacency_matrix(G, type="both", names=TRUE, sparse=FALSE, attr="int")
    miRlist <- unique(prior[, 1])
    Tarlist <- unique(prior[, 2])
    ida.adj <- adj.matrix[which(rownames(adj.matrix) %in% Tarlist), which(colnames(adj.matrix) %in% miRlist)]

    return(ida.adj)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Euclidean <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(euclidean_generic(miRExp, mRExp) + eps(1))
    euclidean.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(euclidean.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Manhattan <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(manhattan_generic(miRExp, mRExp) + eps(1))
    manhattan.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(manhattan.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Canberra <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(canberra_generic(miRExp, mRExp) + eps(1))
    canberra.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(canberra.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Chebyshev <- function(miRExp, mRExp, p.value = 0.05){
   
    original <- 1/(chebyshev_generic(miRExp, mRExp) + eps(1))
    chebyshev.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(chebyshev.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Dice <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(dice_generic(miRExp, mRExp) + eps(1))
    dice.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(dice.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Jaccard <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(jaccard_generic(miRExp, mRExp) + eps(1))
    jaccard.pvalue <- 1 - pnorm(matrixzscore(original))   
    adj.Matrix <- ifelse(jaccard.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Mahalanobis <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- 1/(mahalanobis_generic(miRExp, mRExp) + eps(1))
    mahalanobis.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(mahalanobis.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Cosine <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- cosine_generic(miRExp, mRExp)
    cosine.pvalue <- 1 - pnorm(abs(matrixzscore(original)))   
    adj.Matrix <- ifelse(cosine.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Biweight <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- biweight_generic(miRExp, mRExp)
    biweight.pvalue <- 1 - pnorm(abs(matrixzscore(original)))
    adj.Matrix <- ifelse(biweight.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
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
    wtd_rank.value <- WGCNA::corAndPvalue(savages, method = "pearson")$p
    x.index <- seq(ncol(x))
    y.index <- seq(ncol(y)) + ncol(x)
    wtd_rank.value <- wtd_rank.value[y.index, x.index]
    
    return(wtd_rank.value)

}

## Function of network inference using Weighted_rank
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Weighted_rank <- function(miRExp, mRExp, p.value = 0.05){
    
    cor.pvalue <- wtd_rank_generic(miRExp, mRExp)
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
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
    phit.value <- phit.value[y.index, x.index]    
    
    return(phit.value)
}

## Function of network inference using Phit
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Phit <- function(miRExp, mRExp, p.value = 0.05){

    
    original <- 1/(phit_generic(miRExp, mRExp) + eps(1))
    phit.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(phit.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
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

## Function of network inference using Phis
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Phis <- function(miRExp, mRExp, p.value = 0.05){

    
    original <- 1/(phis_generic(miRExp, mRExp) + eps(1))
    phis.pvalue <- 1 - pnorm(matrixzscore(original))
    adj.Matrix <- ifelse(phis.pvalue < p.value, 1, 0)    
   
    return(adj.Matrix)
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

## Function of network inference using Rhop
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Rhop <- function(miRExp, mRExp, p.value = 0.05){
    
    original <- rhop_generic(miRExp, mRExp)
    phis.pvalue <- 1 - pnorm(abs(matrixzscore(original)))
    adj.Matrix <- ifelse(phis.pvalue < p.value, 1, 0)    
    
    return(adj.Matrix)
}

## Identifying sample-specific miRNA regulation using Scan (Sample-speCific miRNA regulAtioN) and statistical perturbation strategy
# miRExp: miRNA expression data with rows are cells and columns are miRNAs.
# mRExp: mRNA expression data with rows are cells and columns are mRNAs.
# method: Methods for calculating correlations, distances or causal effects.
# pcmethod: Character string specifying the method of the PC algorithm, including "stable", "original", "stable.fast".
# alpha: Significance level (number in [0, 1]) for the individual conditional independence tests.
# p.value: Significance level for the identified miRNA-mRNA interactions.
# num.cores: Number of CPU cores when parallel computation.
# Output: An igraph object, graph of miRNA-mRNA regulatory network for each sample
Scan.perturb <- function(miRExp, mRExp, 
                 method = c("IDA", "Pearson", "Spearman", "Kendall", 
		            "Dcor", "RDC", "Hoeffd", "Zscore", 
			    "MI", "MIC", "Lasso", "Elastic", "Ridge",
			    "Euclidean", "Manhattan", "Canberra",
			    "Chebyshev", "Dice", "Jaccard",
			    "Mahalanobis", "Cosine", 
			    "Biweight", "Weighted_rank", 
			    "Phit", "Phis", "Rhop"), 
                 pcmethod = c("stable", "original", "stable.fast"), 
		 alpha = 0.01,
		 p.value = 0.05,
		 num.cores = 32){
    
     if (method == "IDA"){
        res.all <- IDA(miRExp, mRExp, prior.information, pcmethod = pcmethod, alpha = alpha, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("pcalg", "igraph", "WGCNA"), .export = c("IDA", "querydata")) %dopar% {
            IDA(miRExp[-i, ], mRExp[-i, ], prior.information, pcmethod = pcmethod, alpha = alpha, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()    
    
    } else if (method == "Pearson"){
        res.all <- Pearson(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Pearson", "querydata")) %dopar% {
            Pearson(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Spearman"){
        res.all <- Spearman(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Spearman", "querydata")) %dopar% {
            Spearman(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Kendall"){
        res.all <- Kendall(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Kendall", "querydata")) %dopar% {
            Kendall(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
     

    } else if (method == "Dcor"){
        res.all <- Dcor(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "energy"), .export = c("Dcor", "querydata")) %dopar% {
            Dcor(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "RDC"){
        res.all <- RDC(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = "igraph", .export = c("RDC", "rdc_generic", "querydata")) %dopar% {
            RDC(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Hoeffd"){
        res.all <- Hoeffd(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "Hmisc"), .export = c("Hoeffd", "querydata")) %dopar% {
            Hoeffd(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()    

    } else if (method == "Zscore"){
        res.all <- Zscore(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = "igraph", .export = c("Zscore", "querydata")) %dopar% {
            Zscore(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "MI"){        
        res.all <- MI(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "parmigene"), .export = c("MI", "mi_generic", "querydata", "matrixzscore")) %dopar% {
            MI(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "MIC"){        
        res.all <- MIC(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "minerva"), .export = c("MIC", "mic_generic", "querydata", "matrixzscore")) %dopar% {
            MIC(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Lasso"){        
        res.all <- Lasso(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet"), .export = c("Lasso", "lasso_generic", "querydata", "matrixzscore")) %dopar% {
            Lasso(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Elastic"){        
        res.all <- Elastic(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet"), .export = c("Elastic", "elastic_generic", "querydata", "matrixzscore")) %dopar% {
            Elastic(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Ridge"){        
        res.all <- Ridge(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "glmnet"), .export = c("Ridge", "ridge_generic", "querydata", "matrixzscore")) %dopar% {
            Ridge(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Euclidean"){        
        res.all <- Euclidean(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "pracma"), .export = c("Euclidean", "euclidean_generic", "querydata", "matrixzscore")) %dopar% {
            Euclidean(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Manhattan"){        
        res.all <- Manhattan(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "pracma"), .export = c("Manhattan", "manhattan_generic", "querydata", "matrixzscore")) %dopar% {
            Manhattan(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Canberra"){        
        res.all <- Canberra(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "pracma"), .export = c("Canberra", "canberra_generic", "querydata", "matrixzscore")) %dopar% {
            Canberra(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Chebyshev"){        
        res.all <- Chebyshev(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "stats", "pracma"), .export = c("Chebyshev", "chebyshev_generic", "querydata", "matrixzscore")) %dopar% {
            Chebyshev(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()   

    } else if (method == "Dice"){        
        res.all <- Dice(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy", "pracma"), .export = c("Dice", "dice_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Dice(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Jaccard"){        
        res.all <- Jaccard(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy", "pracma"), .export = c("Jaccard", "jaccard_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Jaccard(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()   
    
    } else if (method == "Mahalanobis"){        
        res.all <- Mahalanobis(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "StatMatch", "pracma"), .export = c("Mahalanobis", "mahalanobis_generic", "mahalanobis.dist", "querydata", "matrixzscore")) %dopar% {
            Mahalanobis(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Cosine"){        
        res.all <- Cosine(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "philentropy"), .export = c("Cosine", "cosine_generic", "distance", "querydata", "matrixzscore")) %dopar% {
            Cosine(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Biweight"){        
        res.all <- Biweight(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "WGCNA"), .export = c("Biweight", "biweight_generic", "querydata", "matrixzscore")) %dopar% {
            Biweight(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Weighted_rank"){        
        res.all <- Weighted_rank(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph"), .export = c("Weighted_rank", "wtd_rank_generic", "querydata")) %dopar% {
            Weighted_rank(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()   
    
    } else if (method == "Phit"){        
        res.all <- Phit(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr", "pracma"), .export = c("Phit", "phit_generic", "phit", "querydata", "matrixzscore")) %dopar% {
            Phit(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Phis"){        
        res.all <- Phis(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr", "pracma"), .export = c("Phis", "phis_generic", "phis", "querydata", "matrixzscore")) %dopar% {
            Phis(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "Rhop"){        
        res.all <- Rhop(miRExp, mRExp, prior.information, p.value = p.value)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)  
        
	res.single <- foreach(i = seq(nrow(miRExp)), .packages = c("igraph", "propr"), .export = c("Rhop", "rhop_generic", "perb", "querydata", "matrixzscore")) %dopar% {
            Rhop(miRExp[-i, ], mRExp[-i, ], prior.information, p.value = p.value)
        }     
    
        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
    } 
    
    int.sample.specific.incident <- lapply(seq(res.single), function(i)
                                         abs(res.single[[i]] - res.all))    
    
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

