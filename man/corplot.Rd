\name{corplot}
\alias{corplot}
\alias{loadplot}


\title{Plots of Correlations and Loadings}

\usage{
corplot(cca.out, X, Y)
loadplot(cca.out, numScores)
}

\description{
  \code{corplot} computes the correlation within and between \code{X} and \code{Y}
                 and displays the three corresponding matrices visusally.

  \code{loadplot} computes the squared loadings for \code{X} and \code{Y} and plots the
                  resulting matrices.  
}

\arguments{
  \item{cca.out}{output from the \code{\link{scca}} or \code{\link{cca}} function.}
    
  \item{X, Y}{input data matrices.}

  \item{numScores}{number of CCA scores shown in plot.}
}

\value{
  A plot.
}

\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).

  Part of the plot code was adapted from the \code{\link[CCA]{img.matcor}} function 
  in the CCA package and from the \code{\link[fields]{image.plot}} function in the fields package.
}
\seealso{\code{\link{scca}}.}

\keyword{plot}

