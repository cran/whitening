\name{pitprops14}
\alias{pitprops14}

\title{Pitprops Correlation Data for 14 Variables}


\usage{
data(pitprops14)
}

\description{
Pit prop timber is used in construction to build mining tunnels. 
The \code{pitprops14} data is described in Jeffers (1967) and is
 a correlation matrix that was calculated from measuring 14 physical properties
of 180 pit props made from wood from Corsican pines grown in East Anglia, UK.
}
\details{
The 14 variables are described in Jeffers (1967) as follows:

\code{topdiam}: the top diameter of the prop in inches;
\code{length}: the length of the prop in inches;
\code{moist}: the moisture content of the prop, expressed as a percentage of the dry weight;
\code{testsg}: the specific gravity of the timber at the time of the test;
\code{ovensg}: the oven-dry specific gravity of the timber;
\code{rinotop}: the number of annual rings at the top of the prop;
\code{ringbut}: the number of annual rings at the base of the prop;
\code{bowmax}: the maximum bow in inches;
\code{bowdist}: the distance of the point of maximum bow from the top of the prop in inches;
\code{whorls}: the number of knot whorls;
\code{clear}: the length of clear prop from the top of the prop in inches;
\code{knots}: the average number of knots per whorl;
\code{diaknot}: the average diameter of the knots in inches;
\code{maxcs}: the maximum compressive strength in lb per square inch.
}

\format{
A correlation matrix of dimension 14 times 14.
}

\source{
The data set is printed in Jeffers (1967) in Table 2 and Table 5.
}

\references{
Jeffers, J. N. R. 1967. Two case studies in the application of principal component analysis.    
JRSS C (Applied Statistics) 16: 225-236. <DOI:10.2307/2985919>
}


\examples{
# load whitening library
library("whitening")

# load pitprops14 data set
data(pitprops14)
colnames(pitprops14)

# correlation matrix for the first 13 variables
pitprops13 = pitprops14[1:13, 1:13]

# correlation loadings for PCA whitening
Psi = whiteningLoadings(pitprops13, "PCA")$Psi

# corresponding explained variation
Psi.explained = explainedVariation(Psi)

# the first six whitened variables account for 87% of the variation
cumsum(Psi.explained)/13*100

}


\keyword{datasets}
