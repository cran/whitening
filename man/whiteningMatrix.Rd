\name{whiteningMatrix}
\alias{whiteningMatrix}
\alias{whiten}

\title{Compute Whitening Matrix and Whiten Data}
\usage{
whiteningMatrix(Sigma, method=c("ZCA", "PCA", "Cholesky",
                                "ZCA-cor", "PCA-cor"))
whiten(X, method=c("ZCA", "PCA", "Cholesky", "ZCA-cor", "PCA-cor"))
}
\arguments{
  \item{Sigma}{Covariance matrix.}
  
  \item{method}{Determines the type of whitening.}

  \item{X}{Data matrix, with samples in rows and variables in columns.}
}

\description{
  \code{whiteningMatrix} computes the whitening matrix \code{W}  corresponding to
   the five natural whitening procedures discussed in  Kessy, Lewin, and Strimmer (2018).

  \code{whiten} whitens data \code{X} using the empirical covariance matrix \eqn{cov(X)} as basis for
   computing the whitening transformation.

}

\value{
  \code{whiteningMatrix} returns a square whitening matrix \code{W}.

  \code{whiten} returns the whitened data matrix \eqn{Z = X W'}.
}

\details{
   ZCA whitening, or Mahalanobis whitening ensures that the average covariance
   between whitened and orginal variables is maximal.  Likewise, ZCA-cor whitening 
   leads to whitened variables that are maximally correlated (on average) with the 
   original variables.  

   In contrast, PCA and PCA-cor whitening lead to maximally compressed 
   whitened variables, as measured by squared covariance and correlation, respectively.

   Cholesky whitening is the unique whitening procedure that results from 
   lower-triangular positive diagonal cross-covariance and cross-correlations   matrices.

   In PCA and PCA-cor eigenvector matrices with positive diagonal are used, 
   in order to resolve the sign-ambiguity and also to make cross-correlations
   and cross-correlations positive diagonal.

   For details see Kessy, Lewin, and Strimmer (2018).

   ZCA-cor whitening is implicitely employed in computing CAT and CAR scores
   (cf. \code{\link[sda]{catscore}} and \code{\link[care]{carscore}}).

   Canonical correlation analysis (CCA) can also be understood as a
   special form form of whitening.
}

\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}) with Agnan Kessy and Alex Lewin.
}
\seealso{\code{\link[sda]{catscore}} and \code{\link[care]{carscore}}.}


\references{
Kessy, A., A. Lewin, and K. Strimmer. 2018.
Optimal whitening and decorrelation. The American Statistician.
\url{https://doi.org/10.1080/00031305.2016.1277159}

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

# whitened data
Z.ZCAcor.1 = tcrossprod(X, W.ZCAcor) 
zapsmall( cov(Z.ZCAcor.1) )

# directly compute whitened data from X
Z.ZCAcor.2 = whiten(X, method="ZCA-cor")
zapsmall( cov(Z.ZCAcor.2) )
}
\keyword{multivariate}
