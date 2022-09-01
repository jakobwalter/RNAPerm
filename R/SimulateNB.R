#' Simulating RNA-Seq Data from NB distribution
#' 
#' This functions allows you to generate RNA-Seq like data from a negative 
#' binomial distribution.
#' @param nSamples number of generated samples. Must be divisible by two as our
#' group sizes are equal
#' @param nGenes number of generated genes.
#'  
#' @keywords cats
#' @export
#' @examples
#' 

SimulateNB <- function(nSamples,
                       nGenes,
                       offset = rep(1, nSamples),
                       mu1,
                       mu2,
                       phi1,
                       phi2,
                       beta = rep(0, nGenes)){


  ## We create the matrices of the two groups separately. Each Matrix
  ## must consist of n_genes rows and n_samples/2 columns.
  betaArray <- exp(t(apply(matrix(c(0, 1), nrow = nSamples/2, ncol = nGenes), 1, function(r) r * beta)))
  mu1array <- t(array(mu1, dim = c(nGenes, nSamples/2))) * betaArray
  phi1array <- t(array(phi1, dim = c(nGenes, nSamples/2)))
  
  ### sample from NB distribution and reshape it to matrix
  Y_1 <- rnbinom(nGenes*nSamples/2, mu = mu1array, size = 1/phi1array)
  Y_1 <- matrix(Y_1, nrow = nSamples/2)
  
  mu2array <- t(array(mu2, dim = c(nGenes, nSamples/2))) * betaArray
  phi2array <- t(array(phi2, dim = c(nGenes, nSamples/2)))

  ### sample from NB distribution and reshape it to matrix
  Y_2 <- rnbinom(nGenes*nSamples/2, mu = mu2array, size = 1/phi2array)
  Y_2 <- matrix(Y_2, nrow = nSamples/2)

  ### create different sample sizes and round back to integer
  Y_all <-rbind(Y_1, Y_2) * offset
  colnames(Y_all) <- names(mu_1)
  mode(Y_all) <- "integer"

  ### create covariate
  X_all <- rep(c(-1,1), each = nSamples/2)

  if(all(beta ==0)){
    dataList <- list("Y" = t(Y_all), "X" = X_all)
  } else{
    dataList <- list("Y" = t(Y_all), "X" = data.frame(list("beta" = X_all,
                                                      "gamma" = as.factor(rep(c(0,1), nSamples/2))
                                                      )
                                                   )
                     )
  }
  return(dataList)

}
