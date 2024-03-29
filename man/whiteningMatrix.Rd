\name{whiteningMatrix}
\alias{whiteningMatrix}

\title{Compute Whitening Matrix}
\usage{
  whiteningMatrix(Sigma, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"))
}
\arguments{
  \item{Sigma}{Covariance matrix.}
  \item{method}{Determines the type of whitening transformation.}
}

\description{
  \code{whiteningMatrix} computes the whitening matrix \eqn{W}.
}

\value{
  \code{whiteningMatrix} returns a square whitening matrix \eqn{W}.
}

\details{
  Whitening is a linear transformation \eqn{z = W x} where the whitening matrix satisfies the constraint \eqn{W^T W = \Sigma^{-1}}  where \eqn{\Sigma = Cov(x)}.

  This function implements various natural whitening transformations discussed in  Kessy, Lewin, and Strimmer (2018).

  The following different whitening approaches can be selected:

  \code{method="ZCA"}: ZCA whitening, also known as Mahalanobis whitening, ensures that the average covariance between whitened and orginal variables is maximal.  

  \code{method="ZCA-cor"}: Likewise, ZCA-cor whitening leads to whitened variables that are maximally correlated (on average) with the original variables.  

  \code{method="PCA"}: In contrast, PCA whitening lead to maximally compressed whitened variables, as measured by squared covariance.

  \code{method="PCA-cor"}: PCA-cor whitening is similar to PCA whitening but uses squared correlations.  

  \code{method="Cholesky"}: computes a whitening matrix by applying Cholesky decomposition.  This yields both a lower triangular positive diagonal whitening matrix and lower triangular positive diagonal loadings (cross-covariance and cross-correlation).  

   Note that Cholesky whitening depends on the ordering of input variables. In the convention used here the first input variable is linked with the first latent variable only, the second input variable is linked to the first and second latent variable only, and so on, and the last variable is linked to all latent variables. 
   
  ZCA-cor whitening is implicitely employed in computing CAT and CAR scores used for variable selection in classification and regression, see the functions \code{catscore} in the \code{sda} package and \code{carscore} in the \code{care} package.

  In both PCA and PCA-cor whitening there is a sign-ambiguity in the eigenvector matrices. In order to resolve the sign-ambiguity we use eigenvector matrices with a positive diagonal so that PCA and PCA-cor cross-correlations and cross-covariances have a positive diagonal for the given ordering of the original variables.

  For details see Kessy, Lewin, and Strimmer (2018).
 
  Canonical correlation analysis (CCA) can also be understood as a special form of whitening (also implemented in this package).
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}) with Agnan Kessy and Alex Lewin.
}
\seealso{\code{\link{whiten}}, \code{\link{whiteningLoadings}}, \code{\link{scca}}.}


\references{
Kessy, A., A. Lewin, and K. Strimmer. 2018.
Optimal whitening and decorrelation. The American Statistician. 72: 309-314.
<DOI:10.1080/00031305.2016.1277159>

}


\examples{
# load whitening library
library("whitening")

# example data set
# E. Anderson. 1935.  The irises of the Gaspe Peninsula.
# Bull. Am. Iris Soc. 59: 2--5
data("iris")
X = as.matrix(iris[,1:4])
d = ncol(X) # 4
n = nrow(X) # 150
colnames(X) # "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"

# estimate covariance
S = cov(X)

# ZCA-cor whitening matrix
W.ZCAcor = whiteningMatrix(S, method="ZCA-cor")

# check constraint on the whitening matrix
crossprod(W.ZCAcor)
solve(S)

# whitened data
Z = tcrossprod(X, W.ZCAcor) 
Z
zapsmall( cov(Z) )
}
\keyword{multivariate}
