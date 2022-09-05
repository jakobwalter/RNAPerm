#' FlipScores test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means for the null hypothesis \eqn{H_0: \beta_j = 0}.
#' @references Hemerik, Jesse, Jelle J. Goeman, and Livio Finos. 
#' "Robust testing in generalized linear models by sign flipping score contributions." 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology) 82.3 (2020): 841-864.
#' @param dge A dgeList object, created with edgeR, containing the normalization factors as computed by edgeR.
#' @param design A model matrix; The first column should be all 1s. The second column should have two unique values, 
#' corresponding to the groups
#' @param scoreType Type of Score contributions on which flipping is performed
#' @param toBeTested index of column of design matrix which is tested 
#' @param nPerm Number of random permutations used for the computation of the p-value
#' @author Jakob Walter
#' @import edgeR
#' @importFrom stats as.formula
#' @importFrom utils tail
#' @importFrom stats glm
#' @importFrom stats residuals
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
#' @export
#' @examples
#' Y <- rnbinom(20*10, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(20, 10)))
#' X1 <- as.factor(rep(c("A", "B"), each = 20/2))
#' design <- model.matrix(~X1, contrasts.arg = list(X1 = "contr.sum"))
#' dge <- edgeR::DGEList(counts = t(Y), group = X1)
#' dge <- edgeR::calcNormFactors(dge)
#' flipScoresTest(dge, design, scoreType = c("basic"), nPerm = 100)

flipScoresTest <- function(dge, design, scoreType = c("basic", "effective"),  toBeTested = 2, nPerm = 5000){
  scoreType <- match.arg(scoreType)
  ## Normalize Data Using EdgeR
  ### Compute Normalization Factors
  os <- edgeR::getOffset(dge)
  
  data <- data.frame(t(dge$counts), design)
  colnames(data) <- c(paste0("Y", 1:nrow(dge$counts)), paste0("X", 1:ncol(design)))
  
  
  ### compute table of scores
  scores <- sapply(1:nrow(dge$counts), function(i){
    ### get Formulas
    formH0 <- as.formula(paste(colnames(data)[i], "~",
                               paste(tail(colnames(data), 2)[-toBeTested], sep = "",  collapse =  " + "),
                               " + ",  "offset(os) - 1", sep = " "), 
    )
    
    formHA <- as.formula(paste(colnames(data)[i], "~",
                               paste(tail(colnames(data), 2), sep = "",  collapse =  " + "),
                               " + ",  "offset(os) - 1", sep = " "), 
    )
    
    
    ### Fit Models
    modH0 <- glm(formH0, family = "poisson", data = data, x = T)
    modHA <- glm(formHA, family = "poisson", data = data, x = T)
    
    if (scoreType == "basic"){
      scoresYi <- modHA$x[,toBeTested]*residuals(modH0, type = "response")
      
    } else if (scoreType == "effective"){
      ### Effective Score Test
      W <- exp(modH0$linear.predictors)          #weight matrix
      infMat <- t(modHA$x * W) %*% modHA$x[,]
      invInfMat <- solve(infMat)
      nu <- invInfMat %*% t(modHA$x*residuals(modH0, type = "response"))
      scoresYi <- nu[toBeTested,]
    } else {
      stop("Unexpected ScoreType")
    }
    return(scoresYi)
  })
  
  
  ### create array to flip score contributions
  flipArray <- array(rbinom(ncol(dge$counts)*nPerm, size = 1, prob = 0.5)*2-1, 
                     dim = c(nPerm,ncol(dge$counts))
  )
  
  ### first row is identity flip
  flipArray[1,] <- 1
  
  
  ### matrix multiplication to compute test statistics
  testStatistics <- abs(flipArray %*% scores)
  
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  
  
  return(pVals)
}
