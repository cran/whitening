\name{whiteningLoadings}
\alias{whiteningLoadings}

\title{Compute Whitening Loadings}

\usage{
  whiteningLoadings(Sigma, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"))
}

\arguments{
  \item{Sigma}{Covariance matrix.}
  \item{method}{Determines the type of whitening transformation.}
}

\description{
  \code{whiteningLoading} computes the loadings (= cross-covariance matrix \eqn{\Phi=Cov(x, z)}) between the original and the whitened variables as well as the correlation loadings (= cross-correlation matrix  \eqn{\Psi=Cor(x, z)}).  The original variables are in the rows and the whitened variables in the columns. 
}

\value{
  \code{whiteningLoadings} returns a list with the following items:

  \code{Phi} - cross-covariance matrix \eqn{\Phi} - the loadings. 

  \code{Psi} - cross-correlation matrix \eqn{\Psi} - the correlation loadings. 
}

\details{
  \eqn{\Phi=Cov(x, z)} is the cross-covariance matrix between the original and the whitened variables. It satisfies \eqn{\Phi \Phi^{T} = \Sigma = Var(x)}. This cross-covariance matrix is the inverse of the whitening matrix so that \eqn{\Phi = W^{-1}}. The cross-covariance matrix is therefore relevant in inverse whitening transformations (=coloring transformations)  \eqn{x = \Phi z}.

\eqn{\Psi=Cor(x, z)} is the cross-correlation matrix between the original and the whitened variables.

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
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{explainedVariation}}, \code{\link{whiten}}, \code{\link{whiteningMatrix}}.}

\references{
Kessy, A., A. Lewin, and K. Strimmer. 2018.
Optimal whitening and decorrelation. The American Statistician. 72: 309-314.
<DOI:10.1080/00031305.2016.1277159>

}


\examples{
# load whitening library
library("whitening")

######

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

# ZCA-cor loadings
ldgs = whiteningLoadings(S, method="ZCA-cor")
ldgs

# cross-covariance matrix
Phi.ZCAcor = ldgs$Phi

# check constraint of cross-covariance matrix
tcrossprod(Phi.ZCAcor)
S

# cross-covariance matrix aka loadings is equal to the inverse whitening matrix
Phi.ZCAcor
solve(W.ZCAcor)
}
\keyword{multivariate}
