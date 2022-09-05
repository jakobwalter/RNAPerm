#' Classic Permutation Test for RNA-Seq Data
#'
#' Computes p-values using a classic permutation test based on the absolute
#' difference in means. The data is normalized using edgeR.
#' @param X A vector encoded as a factor with two levels. 
#' The levels encode the two different classes that are to be tested for differential expression.
#' @param Y An array or dataframe 
#' @param nPerm Number of random permutations used for the computation of the p-value
#' @keywords Internal
#' @import edgeR
#' @export
#' @examples
#' Y <- rnbinom(40*100, mu = 10, size = 1/0.2)
#' Y <- data.frame(array(Y, dim = c(40, 100)))
#' X <- as.factor(rep(c("A", "B"), each = 20))
#' FlipScoresTestNoCovBasic (X, Y, 100)

FlipScoresTest <- function(dge, design, scoreType = "basic",  toBeTested = 2, nPerm = 5000){
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
  flipArray <- array(rbinom(nrow(Y)*nPerm, size = 1, prob = 0.5)*2-1, 
                     dim = c(nPerm,nrow(Y))
  )
  
  ### first row is identity flip
  flipArray[1,] <- 1
  
  
  ### matrix multiplication to compute test statistics
  testStatistics <- abs(flipArray %*% scores)
  
  
  ## compute proportion of permutation test-statistics larger than the observed ones
  pVals <- rowMeans(apply(testStatistics, 1, function(i) i >= testStatistics[1,]))
  
  
  return(pVals)
}