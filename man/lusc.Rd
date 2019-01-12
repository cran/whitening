\name{lusc}
\alias{lusc}

\title{TCGA LUSC Data}


\usage{
data(lusc)
}

\description{
A preprocessed sample of gene expression and methylation data as well as selected clinical covariates for 130 patients with lung squamous cell carcinoma (LUSC) as available from The Cancer Genome Atlas (TCGA) database (Kandoth et al. 2013).
 }

\details{This data set is used to illustrate CCA-based data integration
 in Jendoubi and Strimmer (2019) and also described in Wan et al. (2016).
}

\format{
  \code{lusc$rnaseq2} is a 130 x 206 matrix containing the calibrated gene expression 
    levels of 206 genes for 130 patients.

  \code{lusc$methyl} is a 130 x 234 matrix containing the methylation levels
    of 234 probes for 130 patients.

  \code{sex} is a vector recording the sex (male vs. female) of the 130 patients.
 
  \code{packs} is the number of cigarette packs per year smoked by each patient.

  \code{survivalTime} is number of days to last follow-up or the days to death.

  \code{censoringStatus} is the vital status (0=alive, 1=dead).
}

\source{
The data were retrieved from TCGA (Kandoth et al. 2014) using the TCGA2STAT tool following the guidelines and the preprocessing steps detailed in Wan et al. (2016).
}

\references{
Jendoubi, T., Strimmer, K.: A whitening approach to probabilistic canonical correlation analysis for omics data integration. BMC Bioinformatics 20:15 
<DOI:10.1186/s12859-018-2572-9>

Kandoth, C., McLellan, M.D., Vandin, F., Ye, K., Niu, B., Lu, C., Xie, M., andJ. F. McMichael, Q.Z., Wyczalkowski, M.A., Leiserson, M.D.M., Miller, C.A., Welch, J.S., Walter, M.J., Wendl, M.C., Ley, T.J., Wilson, R.K., Raphael, B.J., Ding, L.: Mutational landscape and significance across 12 major cancer types.
Nature 502, 333--339 (2013). <DOI:10.1038/nature12634>

Wan, Y.-W., Allen, G.I., Liu, Z.: TCGA2STAT: simple TCGA data access for integrated statistical analysis in R. Bioinformatics 32, 952--954 (2016). <DOI:10.1093/bioinformatics/btv677>
}


\examples{
# load whitening library
library("whitening")

# load TGCA LUSC data set
data(lusc)

names(lusc)
#"rnaseq2"         "methyl"          "sex"             "packs"          
#"survivalTime"    "censoringStatus" 

dim(lusc$rnaseq2) # 130 206 gene expression
dim(lusc$methyl)  # 130 234 methylation level

\dontrun{
library("survival")
s = Surv(lusc$survivalTime, lusc$censoringStatus)
plot(survfit(s ~ lusc$sex), xlab = "Years", ylab = "Probability of survival", lty=c(2,1), lwd=2)
legend("topright", legend = c("male", "female"), lty =c(1,2), lwd=2)
}

}

\keyword{datasets}
