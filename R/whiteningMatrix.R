### whiteningMatrix.R  (2020-12-05)
###
###    Compute whitening matrix
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


makePositivDiagonal = function(U) return( sweep(U, 2, sign(diag(U)), "*") ) # U %*% diag( sign(diag(U)) ) 

# create whitening matrix W from given covariance matrix Sigma
whiteningMatrix = function(Sigma, 
  method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Chol-prec", "Chol-cov", "Cholesky"))
{
  method=match.arg(method)

  if(method=="Cholesky") method="Chol-prec"

  if(method=="ZCA")
  {
    eSigma = eigen(Sigma, symmetric=TRUE) 
    U = eSigma$vectors
    lambda = eSigma$values
    
    W = U %*% diag(1/sqrt(lambda)) %*% t(U)
  }

  if(method=="PCA")
  {
    eSigma = eigen(Sigma, symmetric=TRUE) 
    U = eSigma$vectors
    lambda = eSigma$values
    
    # fix sign ambiguity in eigenvectors by making U positive diagonal
    U = makePositivDiagonal(U)

    W = diag(1/sqrt(lambda)) %*% t(U)
  }

  if(method=="Chol-prec")
  {
     W = chol(solve(Sigma))
  }

  if(method=="Chol-cov")
  {
     W = solve(t(chol(Sigma)))
  }

  if(method=="ZCA-cor")
  {
    v = diag(Sigma)
    R = cov2cor(Sigma)
    eR = eigen(R, symmetric=TRUE) 
    G = eR$vectors
    theta = eR$values

    W = G %*% diag(1/sqrt(theta)) %*% t(G) %*% diag(1/sqrt(v))
  }

  if(method=="PCA-cor")
  {
    v = diag(Sigma)
    R = cov2cor(Sigma)
    eR = eigen(R, symmetric=TRUE) 
    G = eR$vectors
    theta = eR$values

    # fix sign ambiguity in eigenvectors by making G positive diagonal
    G = makePositivDiagonal(G)

    W = diag(1/sqrt(theta)) %*% t(G) %*% diag(1/sqrt(v))
  }

  return (W)
}

