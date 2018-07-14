### scca.R  (2018-07-09)
###
###    Estimate canonical correlations, canonical directions, and loading
###
### Copyright 2018 Korbinian Strimmer
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



# empirical version - same as scca but with lambda.cor=0 and verbose=FALSE
cca = function(X, Y, scale=TRUE)
{
  out = scca(X, Y, lambda.cor=0, scale=scale, verbose=FALSE)

  return(out)
}

scca = function(X, Y, lambda.cor, scale=TRUE, verbose=TRUE)
{
  if (!is.matrix(X)) stop("X needs to be a matrix!")
  if (!is.matrix(Y)) stop("Y needs to be a matrix!")

  p = dim(X)[2]
  q = dim(Y)[2]
  n = dim(X)[1]
  if (dim(Y)[1] !=n ) stop("X and Y needs to have the same number of rows!")

  # find regularization parameter for the correlation matrix
  if( missing(lambda.cor) )
  {
    lambda.cor.estimated = TRUE
    # regularize the joint correlation matrix (X and Y combined)
    lambda.cor = estimate.lambda( cbind(X,Y), verbose=verbose )
  }
  else
  {
    lambda.cor.estimated = FALSE
  }

 
  # memory and time-saving way of computing cross products of shrinkage correlation matrix
  # RX^{alpha} %*% m
  cpLeft = function(x, m, alpha, lambda.cor)
  {
    a = crossprod.powcor.shrink(x, m, alpha=alpha, lambda=lambda.cor, verbose = FALSE)
    attr(a, "class") = NULL
    attr(a, "lambda.estimated") = NULL
    attr(a, "lambda") = NULL
    return(a)
  }

  # m %*% RX^{alpha}
  cpRight = function(x, m, alpha, lambda.cor)
    t( cpLeft(x, t(m), alpha, lambda.cor ) )


  # shrinkage estimate of cross-correlation
  RXY = (1-lambda.cor)*cor(X, Y)

  # adjusted cross-correlations pxq    K = RXX^(-1/2) %*% RXY %*% isqrtRYY^(-1/2) 
  # this avoids to compute the full matrices RXX^(-1/2) and RYY^(-1/2)
  K = cpRight(Y, cpLeft(X, RXY, -1/2, lambda.cor), -1/2, lambda.cor)  


  # decompose adjusted cross-correlations K into rotation matrices and canonical correlation
  #svd.out = fast.svd(K) 
  svd.out = svd(K) 
  QX = t(svd.out$u) # t(U) mxp
  QY = t(svd.out$v) # t(V) mxq
  lambda = svd.out$d  # canonical correlations (singular values of K)  m
  # check decomposition: sum( (K - t(QX) %*% diag(lambda) %*% QY)^2 )

  # fix sign ambiguity by making QX and QY positive diagonal
  sQX = sign(diag(QX)) # m
  sQY = sign(diag(QY)) # m
  QX = sweep(QX, 1, sQX, "*") # diag(sQX) %*% QX 
  QY = sweep(QY, 1, sQY, "*") # diag(sQY) %*% QY 
  lambda = lambda*sQX*sQY # some canonical correlations may now be negative
  # SVD decomposition is still valid: sum( (K - t(QX) %*% diag(lambda) %*% QY)^2 )
 

  # compute canonical directions

  # canonical directions for X  (rows contain directions) = t(Xcoef)
  WX = cpRight(X, QX, -1/2, lambda.cor) # QX %*% RXX^(-1/2) (mxp)

  # canonical directions for Y  (rows contain directions) = t(Ycoef)
  WY = cpRight(Y, QY, -1/2, lambda.cor) # QY %*% RYY^(-1/2) (mxq)


  # bring back scale into the canonical directions if data should not be scaled
  if( !scale)
  {
    SDX = apply(X, 2, sd) # standard deviations for X
    SDY = apply(Y, 2, sd) # standard deviations for Y
    WX = sweep(WX, 2, SDX, "/")
    WY = sweep(WY, 2, SDY, "/")
  }

  colnames(WX) = colnames(X)
  colnames(WY) = colnames(Y)

  # compute loadings

  # loadings for X
  PhiX = cpRight(X, QX, 1/2, lambda.cor) # QX %*% RXX^(-1/2) (mxp)

  # loadings for Y
  PhiY = cpRight(Y, QY, 1/2, lambda.cor) # QY %*% RYY^(-1/2) (mxq)

  # bring back scale into the loadings if data should not be scaled
  if( !scale)
  {
    WX = sweep(PhiX, 2, SDX, "*")
    WY = sweep(PhiY, 2, SDY, "*")
  }

  colnames(PhiX) = colnames(X)
  colnames(PhiY) = colnames(Y)

  return( list(K=K, lambda=lambda, WX=WX, WY=WY, PhiX=PhiX, PhiY=PhiY, scale=scale,
          lambda.cor.estimated=lambda.cor.estimated,
          lambda.cor=lambda.cor) )
}

