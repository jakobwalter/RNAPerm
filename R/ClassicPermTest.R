#' Classic Permutation Test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means. The data is normalized using edgeR.
#' @param X A vector encoded as a factor with two levels. 
#' The levels encode the two different classes that are to be tested for differential expression.
#' @param Y An array or dataframe 
#' @param nPerm Number of random permutations used for the computation of the p-value
#' @import edgeR
#' @export
#' @examples
#' Y <- rnbinom(40*100, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(40, 100)))
#' X <- as.factor(rep(c("A", "B"), each = 20))
#' ClassicPermTest(X, Y, 100)


ClassicPermTest <- function(X, Y, nPerm){
  ## Normalize Data Using EdgeR
  
  ### Compute Normalization Factors
  d <- edgeR::DGEList(counts = Y)
  d <- edgeR::calcNormFactors(d)
  
  ### Normalize Data
  Y <- d$counts * d$samples$norm.factors
  
  ## Create empty array for test statistics
  testStatistics <- array(dim = c(nPerm, ncol(Y)))
  
  ## first test-statistics are the ones observed
  testStatistics[1,] <- abs(colMeans(Y[X == levels(X)[1],]) - colMeans(Y[X == levels(X)[2],]))
  
  ## Compute permutation test-statistics using nPerm-1 permutations
  for (i in 2:nPerm){
    XPerm <- sample(X)
    testStatistics[i,] <- abs(colMeans(Y[XPerm == levels(X)[1],]) - colMeans(Y[XPerm == levels(X)[2],]))
  }
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  
  return(pVals)
}