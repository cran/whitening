\name{explainedVariation}
\alias{explainedVariation}

\title{Compute Explained Variation from Loadings}

\usage{
  explainedVariation(Phi)
}

\arguments{
  \item{Phi}{Loading matrix (with columns referring to whitened variables).}
}


\description{
  \code{explainedVariation} computes the explained variation for each whitened variables from the loadings (both covariance loadings and correlation loadings).
}

\value{
  \code{explainedVariation} returns a vector with the explained variation contributed by each whitened variable.

}

\details{
  \code{explainedVariation} computes for each column of
the loading matrix the sum of squares of the elements in that column.
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{whiteningLoadings}}}

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

# PCA-cor loadings
ldgs = whiteningLoadings(S, method="PCA-cor")

# Explained variation from correlation loadings
explainedVariation( ldgs$Psi )


}
\keyword{multivariate}
