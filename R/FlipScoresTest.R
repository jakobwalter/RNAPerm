#' Classic Permutation Test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means. The data is normalized using edgeR.
#' @param X A vector encoded as a factor with two levels. 
#' The levels encode the two different classes that are to be tested for differential expression.
#' @param Y An array or dataframe 
#' @param nPerm Number of random permutations used for the computation of the p-value
#' @keyword Internal
#' @import edgeR
#' @export
#' @examples
#' Y <- rnbinom(40*100, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(40, 100)))
#' X <- as.factor(rep(c("A", "B"), each = 20))
#' ClassicPermTest(X, Y, 100)


FlipScoresTestNoCovBasic <- function(X, Y, nPerm){
  ## Normalize Data Using EdgeR
  
  ### Compute Normalization Factors
  d <- edgeR::DGEList(counts = Y)
  d <- edgeR::calcNormFactors(d)
  os <- edgeR::getOffset(d)
  
  ### Recode X to be in (-1, 1)
  XRecoded <- (X == X[1])*2-1
  
  ### compute table of scores
  scores <- apply(Y, 2, function(YCol){
    modNull<- glm(YCol ~ 1, family = poisson)
    XRecoded*residuals(modNull, type = "response")
  }
  )
  
  ### create array to flip score contributions
  flipArray <- array(rbinom(nrow(Y)*nPerm, size = 1, prob = 0.5)*2-1, dim = c(20,40))
  
  ### first row is identity flip
  flipArray[1,] <- 1
  
  ### matrix multiplication to compute test statistics
  testStatistics <- flipArray %*% scores
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  return(pVals)
}