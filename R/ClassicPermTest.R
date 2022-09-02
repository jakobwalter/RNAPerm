#' Classic Permutation Test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means. The data is normalized using edgeR.
#' @param X A vector encoded as a factor with two levels. 
#' The levels encode the two different classes that are to be tested for differential expression.
#' @param Y An array or dataframe 
#' @keyword absolute difference in means permutation test RNA-Seq
#' @importFrom base sum
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }



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
  testStatistics[1,] <- colMeans(Y[X == levels(X)[1],]) - colMeans(Y[X == levels(X)[2],])
  
  ## Compute permutation test-statistics using nPerm-1 permutations
  for (i in 2:nPerm){
    XPerm <- sample(X)
    testStatistics[i,] <- colMeans(Y[XPerm == levels(X)[1],]) - colMeans(Y[XPerm == levels(X)[2],])
  }
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  
  return(pVals)
}