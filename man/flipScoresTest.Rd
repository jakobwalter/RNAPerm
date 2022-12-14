% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlipScoresTest.R
\name{flipScoresTest}
\alias{flipScoresTest}
\title{FlipScores test for RNA-Seq Data}
\usage{
flipScoresTest(
  dge,
  design,
  scoreType = c("basic", "effective"),
  toBeTested = 2,
  nPerm = 5000
)
}
\arguments{
\item{dge}{A dgeList object, created with edgeR, containing the normalization factors as computed by edgeR.}

\item{design}{A model matrix; The first column should be all 1s. The second column should have two unique values,
corresponding to the groups}

\item{scoreType}{Type of Score contributions on which flipping is performed}

\item{toBeTested}{index of column of design matrix which is tested}

\item{nPerm}{Number of random permutations used for the computation of the p-value}
}
\description{
Computes p-values using a classic permutation test based on the absolute
difference in means for the null hypothesis \eqn{H_0: \beta_j = 0}.
}
\examples{
Y <- rnbinom(20*10, mu = 10, size = 1/0.2)
Y <- data.frame(array(Y, dim = c(20, 10)))
X1 <- as.factor(rep(c("A", "B"), each = 20/2))
design <- model.matrix(~X1, contrasts.arg = list(X1 = "contr.sum"))
dge <- edgeR::DGEList(counts = t(Y), group = X1)
dge <- edgeR::calcNormFactors(dge)
flipScoresTest(dge, design, scoreType = c("basic"), nPerm = 100)
}
\references{
Hemerik, Jesse, Jelle J. Goeman, and Livio Finos.
"Robust testing in generalized linear models by sign flipping score contributions."
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 82.3 (2020): 841-864.
}
\author{
Jakob Walter
}
