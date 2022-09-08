#' MWU test for RNA-Seq Data
#'
#' Computes p-values using a Mann-Whitney U Test for the null hypothesis \eqn{H_0: \beta_j = 0}.
#' @param dge A dgeList object, created with edgeR, containing the normalization factors as computed by edgeR.
#' @param design A model matrix; The first column should be all 1s. The second column should have two unique values, 
#' corresponding to the groups
#' @author Jakob Walter
#' @import edgeR
#' @importFrom stats wilcox.test
#' @export
#' @examples
#' Y <- rnbinom(20*10, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(20, 10)))
#' X1 <- as.factor(rep(c("A", "B"), each = 20/2))
#' design <- model.matrix(~X1, contrasts.arg = list(X1 = "contr.sum"))
#' dge <- edgeR::DGEList(counts = t(Y), group = X1)
#' dge <- edgeR::calcNormFactors(dge)
#' MWUTest(dge, design)


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
  
  pVals <- apply(Y, 1, function(rowi){
    suppressWarnings(wilcox.test(rowi[design[,2] == design[1,2]], rowi[design[,2] == design[1,2]])$p.value)
    }
    )

  return(as.vector(pVals))
}