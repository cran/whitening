\name{forina1986}
\alias{forina1986}

\title{Forina 1986 Wine Data - Extended UCI Wine Data}


\usage{
data(forina1986)
}

\description{
The forina1986 dataset describes 27 properties of 178 samples of wine from three 
grape varieties (59 Barolo, 71 Grignolino, 48 Barbera) as reported in Forina et al. (1986).   
}
\details{This data set contains the full set of covariates described in Forina et al. (1986) 
except for Sulphate (variable 15 in Forina et al. 1986). These are:
1) Alcohol,
2) Sugar-free extract,
3) Fixed acidity,
4) Tartaric acid,
5) Malic acid,
6) Uronic acids,
7) pH,
8) Ash,
9) Alcalinity of ash,
10) Potassium,
11) Calcium,
12) Magnesium,
13) Phosphate,
14) Chloride,
15) Total phenols,
16) Flavanoids,
17) Nonflavanoid phenols,
18) Proanthocyanins,
19) Color intensity,
20) Hue,
21) OD280/OD315 of diluted wines,
22) OD280/OD315 of flavonoids,
23) Glycerol,
24) 2-3-butanediol,
25) Total nitrogen,
26) Proline,
and 27) Methanol.

The UCI wine data set (\url{https://archive.ics.uci.edu/ml/datasets/wine}) 
is a subset of the Forina et al. (1986) data set comprising only 13 variables.
}

\format{
A list containing the following components:

  \code{attribs}  collects measurements for 27 attributes of 178 wine samples.

  \code{type} describes the variety ("Barolo", "Grignolino", or "Barbera").
}

\source{
The original data matrix is available from 
\url{https://www.researchgate.net/publication/271908647_Wines_MForina_CArmanino_MCastino_MUbigli_Multivariate_data_analysis_as_discriminating_method_of_the_origin_of_wines_Vitis_25_189-201_1986}.
}

\references{
Forina, M., Armanino, C., Castino, M., and Ubigli, M.
Multivariate data analysis as a discriminating method of the
origin of wines. Vitis 25:189-201 (1986). \url{https://ojs.openagrar.de/index.php/VITIS/article/view/5950}.
}


\examples{
# load whitening library
library("whitening")

# load Forina 1986 wine data set
data(forina1986)

table(forina1986$type)
#    Barolo Grignolino    Barbera 
#        59         71         48

dim(forina1986$attrib)
# 178  27

colnames(forina1986$attrib)
# [1] "Alcohol"                      "Sugar-free extract"          
# [3] "Fixed acidity"                "Tartaric acid"               
# [5] "Malic acid"                   "Uronic acids"                
# [7] "pH"                           "Ash"                         
# [9] "Alkalinity of ash"            "Potassium"                   
#[11] "Calcium"                      "Magnesium"                   
#[13] "Phosphate"                    "Chloride"                    
#[15] "Total phenols"                "Flavanoids"                  
#[17] "Nonflavanoid phenols"         "Proanthocyanins"             
#[19] "Color intensity"              "Hue"                         
#[21] "OD280/OD315 of diluted wines" "OD280/OD315 of flavonoids"   
#[23] "Glycerol"                     "2-3-butanediol"              
#[25] "Total nitrogen"               "Proline"                     
#[27] "Methanol"  

# PCA-cor whitened data
Z = whiten(forina1986$attrib, method="PCA-cor") 

wt = as.integer(forina1986$type)
plot(Z[,1], Z[,2], xlab=expression(paste(Z[1])), ylab=expression(paste(Z[2])), 
  main="Forina 1986 Wine Data", sub="PCA-cor Whitening", col=wt, pch=wt+14)
legend("topright", levels(forina1986$type)[1:3], col=1:3, pch=(1:3)+14 )



## relationship to UCI wine data

# UCI wine data is a subset
uciwine.attrib = forina1986$attrib[, c("Alcohol", "Malic acid", "Ash", 
  "Alcalinity of ash", "Magnesium", "Total phenols", "Flavanoids", 
  "Nonflavanoid phenols", "Proanthocyanins", "Color intensity", "Hue", 
  "OD280/OD315 of diluted wines", "Proline")]

# two small differences compared to UCI wine data matrix
uciwine.attrib[172,"Color intensity"]  # 9.9 but 9.899999 in UCI matrix
uciwine.attrib[71,"Hue"] # 0.91 but 0.906 in UCI matrix
}


\keyword{datasets}
