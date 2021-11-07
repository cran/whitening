\name{whiteningCrossCov}
\alias{whiteningCrossCov}

\title{Compute Whitening Cross-Covariance Matrix}
\usage{
whiteningCrossCov(Sigma, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", 
    "Chol-prec", "Chol-cov", "Cholesky"))

}
\arguments{
  \item{Sigma}{Covariance matrix.}
  
  \item{method}{Determines the type of whitening transformation.}
}

\description{
 \code{whiteningCrossCov} computes the cross-covariance \eqn{Phi=Cov(z, x)} between the whitened and the original variables.

}

\value{
  \code{whiteningCrossCov} returns a square cross-covariance matrix \eqn{Phi}.   

}

\details{
   \eqn{Phi=Cov(z, x)} is the cross-covariance matrix between the whitened and the original variables.
    It satisfies \eqn{Phi' Phi = Sigma}.   The transpose of the cross-covariance matrix is the inverse of the whitening matrix so that \eqn{Phi' = W^{-1}}.  The cross-covariance matrix is therefore relevant in inverse
whitening transformations (=coloring transformations)  \eqn{x = Phi' z}.

   The following six different whitening approaches can be selected:

   \code{method="ZCA"}: ZCA whitening, also known as Mahalanobis whitening, ensures that the average covariance between whitened and orginal variables is maximal.  

   \code{method="ZCA-cor"}: Likewise, ZCA-cor whitening leads to whitened variables that are maximally correlated (on average) with the original variables.  

   \code{method="PCA"}: In contrast, PCA whitening lead to maximally compressed whitened variables, as measured by squared covariance.

  \code{method="PCA-cor"}: PCA-cor whitening is similar to PCA whitening but uses squared correlations.  

  \code{method="Chol-prec"} and \code{method="Cholesky"}: computes a whitening matrix by applying Cholesky decomposition on the precision matrix.  This yields an upper triangular positive diagonal whitening matrix and lower triangular positive diagonal cross-covariance and cross-correlation matrices.
   
  \code{method="Chol-cov"}: computes a whitening matrix by applying Cholesky decomposition on the precision matrix.  This yields a lower triangular positive diagonal whitening matrix and upper triangular positive diagonal cross-covariance and cross-correlation matrices.
   

   ZCA-cor whitening is implicitely employed in computing CAT and CAR scores used for variable selection in classification and regression, see the functions \code{catscore} in the \code{sda} package and \code{carscore} in the \code{care} package.

   In both PCA and PCA-cor whitening there is a sign-ambiguity in the eigenvector matrices. In order to resolve the sign-ambiguity we use eigenvector matrices with a positive diagonal.  This has the effect to make cross-correlations and cross-correlations positive diagonal for PCA and PCA-cor.

   For details see Kessy, Lewin, and Strimmer (2018).
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{whiteningMatrix}}.}


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

# ZCA-cor cross-covariance matrix
Phi.ZCAcor = whiteningCrossCov(S, method="ZCA-cor")

# check constraint of cross-covariance matrix
crossprod(Phi.ZCAcor)
S

# transpose of cross-covariance matrix is equal to the inverse whitening matrix
t(Phi.ZCAcor)
solve(W.ZCAcor)
}
\keyword{multivariate}
