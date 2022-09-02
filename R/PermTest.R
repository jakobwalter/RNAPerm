ClassicPermTest <- function(X, Y, nPerm){
  ## Normalize Data Using EdgeR
  
  ### Compute Normalization Factors
  d <- edgeR::DGEList(counts = X)
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