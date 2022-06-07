\name{whiten}
\alias{whiten}

\title{Whiten Data Matrix}
\usage{
  whiten(X, center=FALSE, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"))
}
\arguments{
  \item{X}{Data matrix, with samples in rows and variables in columns.}
 
  \item{center}{Center columns to mean zero.}

  \item{method}{Determines the type of whitening transformation (see Details).}
}

\description{
  \code{whiten} whitens a data matrix \eqn{X} using the empirical covariance matrix \eqn{cov(X)} as basis for computing the whitening transformation.
}

\value{
  \code{whiten} returns the whitened data matrix \eqn{Z = X W^T}.
}

\details{

  The following whitening approaches can be selected:

  \code{method="ZCA"} and \code{method="ZCA-cov"}: ZCA whitening, also known as Mahalanobis whitening, ensures that the average covariance between whitened and orginal variables is maximal.  

  \code{method="ZCA-cor"}: Likewise, ZCA-cor whitening leads to whitened variables that are maximally correlated (on average) with the original variables.  

  \code{method="PCA"} and \code{method="PCA-cov"}: In contrast, PCA whitening lead to maximally compressed whitened variables, as measured by squared covariance.

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
\seealso{\code{\link{whiteningMatrix}}, \code{\link{whiteningLoadings}}, \code{\link{scca}}.}


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

# whitened data
Z.ZCAcor = whiten(X, method="ZCA-cor")

# check covariance matrix
zapsmall( cov(Z.ZCAcor) )
}
\keyword{multivariate}
