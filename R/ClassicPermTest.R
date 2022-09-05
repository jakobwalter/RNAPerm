#' Classic Permutation Test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means.
#' @param dge A dgeList object, created with edgeR, containing the normalization factors as computed by edgeR.
#' @param design A model matrix; The first column should be all 1s. The second column should have two unique values, 
#' corresponding to the groups
#' @param nPerm Number of random permutations used for the computation of the p-value
#' @import edgeR
#' @export
#' @examples
#' Y <- rnbinom(20*10, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(20, 10)))
#' X1 <- as.factor(rep(c("A", "B"), each = 20/2))
#' design <- model.matrix(~X1, contrasts.arg = list(X1 = "contr.sum"))
#' dge <- edgeR::DGEList(counts = t(Y), group = X1)
#' dge <- edgeR::calcNormFactors(dge)
#' pClassic <- classicPermTest(dge, design, 2000)



classicPermTest <- function(dge, design, nPerm){
  if (ncol(design) > 2  | length(unique(design[,2])) != 2){
    stop("Classical Permutation Test only works for two-group comparisons!")
  }
  if (class(dge)[1] != "DGEList"){
    stop("Class of DGE needs to be DGEList (edgeR)")
  }
  if (all(dge$samples$norm.factors == 1)){
    warning("Data might not be normalized!")
  }
  
  ## Normalize Data Using EdgeR
  Y <- dge$counts * dge$samples$norm.factors
  
  ## Create empty array for test statistics
  testStatistics <- array(dim = c(nPerm, nrow(Y)))
  
  ## first test-statistics are the ones observed
  testStatistics[1,] <- abs(rowMeans(Y[,design[,2] == design[1,2]]) - rowMeans(Y[,design[,2] != design[1,2]]))
  
  ## Compute permutation test-statistics using nPerm-1 permutations
  for (i in 2:nPerm){
    XPerm <- sample(design[,2])
    testStatistics[i,] <- abs(rowMeans(Y[,XPerm == design[1,2]]) - rowMeans(Y[,XPerm != design[1,2]]))
  }
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  
  return(pVals)
}