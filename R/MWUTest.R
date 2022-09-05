#' MWU-Test for RNA-Seq Data
#'
#' Computes p-values using a Mann Whitney U-test. The data is normalized using edgeR.
#' @param X A vector encoded as a factor with two levels. 
#' The levels encode the two different classes that are to be tested for differential expression.
#' @param Y An array or dataframe 
#' @import edgeR
#' @export
#' @examples
#' Y <- rnbinom(40*100, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(40, 100)))
#' X <- as.factor(rep(c("A", "B"), each = 20))
#' MWUTest(X, Y, 100)


MWUTest <- function(dge, design){ 
  if (ncol(design) > 2){
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
  
  ### Compute Normalization Factors
  d <- edgeR::DGEList(counts = Y)
  d <- edgeR::calcNormFactors(d)
  
  ### Normalize Data
  Y <- d$counts * d$samples$norm.factors
  
  pVals <- apply(Y, 1, function(rowi){
    suppressWarnings(wilcox.test(rowi[design[,2] == design[1,2]], rowi[design[,2] == design[1,2]])$p.value)
    }
    )

  return(as.vector(pVals))
}