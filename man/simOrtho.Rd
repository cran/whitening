\name{simOrtho}
\alias{simOrtho}

\title{Simulate Random Orthogonal Matrix}
\usage{
simOrtho(d, nonNegDiag = FALSE)
}

\arguments{
  \item{d}{The dimension of the orthogonal matrix.} 
  \item{nonNegDiag}{force the elements on the diagonal to be nonnegative (default: FALSE).} 
}

\description{
  \code{simOrtho} generates a random orthogonal matrix.
}

\value{
  \code{simOrtho} returns a real matrix of size \eqn{d} times \eqn{d}.
}

\details{

The algorithm is based on QR decomposition of a random matrix, see Section 3 page 404 in Stewart (1980).
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{whiteningMatrix}}}


\references{
G.W. Stewart. 1980. The efficient generation of random orthogonal matrices with an application to condition estimators. SIAM J. Numer. Anal.  17:403-409. <DOI:10.1137/0717034>

}


\examples{
# load whitening library
library("whitening")

# simulate random orthogonal matrix
Q = simOrtho(4) # matrix of dimension 4x4
Q

zapsmall( crossprod(Q) )
zapsmall( tcrossprod(Q) )


}
\keyword{multivariate}
