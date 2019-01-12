\name{nutrimouse}
\alias{nutrimouse}

\title{Nutrimouse Data}


\usage{
data(nutrimouse)
}

\description{
The nutrimouse dataset is a collection of gene expression and lipid measurements 
collected in a nutrigenomic study in the mouse studying 40 animals by Martin et al. (2007).
}
\details{This data set is used to illustrate CCA-based data integration
 in Jendoubi and Strimmer (2019) and is also described in Gonzalez et al. (2008).
}

\format{
A list containing the following components:

  \code{gene}  collects gene expression of 120 genes in liver tissue for 40 mice.

  \code{lipid} collects concentrations of 21 lipids for 40 mice.

  \code{diet} describes the diet of each mouse ("coc", "fish", "lin", "ref", or  "sun").

  \code{genotype} describes the genotype of each mouse: wild type ("wt") or PPARalpha deficient ("ppar").
}

\source{
The original data are available in the CCA R package by Gonzalez et al. (2008), 
see \code{\link[CCA]{nutrimouse}}.
}

\references{
Gonzalez, I., Dejean, S., Martin, P.G.P,  Baccini, A.
CCA: an R package to extend canonical correlation analysis. J. Statist. Software 23:1--13 (2008)

Jendoubi, T., Strimmer, K.: A whitening approach to probabilistic canonical correlation analysis for omics data integration. BMC Bioinformatics 20:15 
<DOI:10.1186/s12859-018-2572-9>

Martin, P.G.P., Guillou, H., Lasserre, F., Dejean, S., Lan, A., Pascussi, J.-M., Cristobal, M.S., Legrand, P., Besse, P., Pineau, T.: Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed
through a multigenomic study. Hepatology 54, 767--777 (2007)
 <DOI:10.1002/hep.21510>
}


\examples{
# load whitening library
library("whitening")

# load nutrimouse data set
data(nutrimouse)

dim(nutrimouse$gene) # 40 120
dim(nutrimouse$lipid) # 40 21
levels( nutrimouse$diet ) #  "coc"  "fish" "lin"  "ref"  "sun"
levels( nutrimouse$genotype ) # "wt"   "ppar"
}


\keyword{datasets}
