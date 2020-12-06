### whiten.R  (2020-12-05)
###
###    Whitening data matrix
###
### Copyright 2018-20 Korbinian Strimmer
###
###
### This file is part of the `whitening' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# whiten data using empirical covariance matrix
whiten = function(X, center=FALSE, method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Chol-prec", "Chol-cov", "Cholesky"))
{
    method = match.arg(method)

    S = cov(X)
    W = whiteningMatrix(S, method=method)
    Z = tcrossprod(X, W) # whitened data

    if(center) Z = sweep(Z, 2, colMeans(Z))

    return(Z)
}
