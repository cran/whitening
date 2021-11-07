\name{whiten}
\alias{whiten}

\title{Whiten Data Matrix}
\usage{
whiten(X, center=FALSE, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor",
    "Chol-prec", "Chol-cov", "Cholesky"))
}
\arguments{
  \item{X}{Data matrix, with samples in rows and variables in columns.}
 
 \item{center}{Center columns to mean zero.}

  \item{method}{Determines the type of whitening transformation.}

}

\description{

  \code{whiten} whitens a data matrix \eqn{X} using the empirical covariance matrix \eqn{cov(X)} as basis for computing the whitening transformation.

}

\value{

  \code{whiten} returns the whitened data matrix \eqn{Z = X W'}.

}

\details{

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
  Korbinian Strimmer (\url{https://strimmerlab.github.io}) with Agnan Kessy and Alex Lewin.
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

# whitened data
Z.ZCAcor = whiten(X, method="ZCA-cor")

# check covariance matrix
zapsmall( cov(Z.ZCAcor) )
}
\keyword{multivariate}
