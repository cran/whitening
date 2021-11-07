\encoding{latin1}
\name{whitening-package}
\alias{whitening-package}
\docType{package}
\title{The whitening Package}

\description{
The "whitening" package implements the whitening methods (ZCA, PCA, Cholesky,
 ZCA-cor, and PCA-cor) discussed in Kessy, Lewin, and Strimmer (2018)
 as well as the whitening approach to canonical correlation analysis allowing negative canonical correlations described in Jendoubi and Strimmer (2019).
}

\author{Korbinian Strimmer (\url{https://strimmerlab.github.io/}) with Takoua Jendoubi, Agnan Kessy, and Alex Lewin.
}
\references{

Kessy, A., A. Lewin, and K. Strimmer. 2018.
Optimal whitening and decorrelation. The American Statistician. 72: 309-314.
<DOI:10.1080/00031305.2016.1277159>

Jendoubi, T., and K. Strimmer 2019. A whitening approach to probabilistic canonical correlation analysis for omics data integration. 
BMC Bioinformatics 20: 15.
<DOI:10.1186/s12859-018-2572-9>

Website: \url{https://strimmerlab.github.io/software/whitening/}
}
\keyword{multivariate}
\seealso{
\code{\link{whiteningMatrix}},
\code{\link{whiten}},
\code{\link{whiteningCrossCov}},
\code{\link{cca}}, and
\code{\link{scca}}.
}
