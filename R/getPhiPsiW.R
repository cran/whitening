### getPhiPsiW.R  (2022-06-02)
###
###    Compute loadings (Phi, Psi) and whitening matrix (W)
###
### Copyright 2018-22 Korbinian Strimmer
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


makePosDiag = function(U) return( sweep(U, 2, sign(diag(U)), "*") ) # U %*% diag( sign(diag(U)) ) 

# compute whitening matrix W from given covariance matrix Sigma
# or the loadings matrix Phi  (= inverse of W)

getPhiPsiW = function(Sigma, 
  method=c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"),
  returnW=FALSE, returnPhiPsi=FALSE)
{
  method=match.arg(method)

  v = diag(Sigma)
  if(method=="ZCA" | method=="PCA")
  {
    eSigma = eigen(Sigma, symmetric=TRUE) 
    U = eSigma$vectors
    lambda = eSigma$values
  }
  if(method=="ZCA-cor" | method=="PCA-cor")
  {
    R = cov2cor(Sigma)
    eR = eigen(R, symmetric=TRUE) 
    G = eR$vectors
    theta = eR$values
  }

  if(method=="ZCA")
  {
    if (returnW)
      W = U %*% diag(1/sqrt(lambda)) %*% t(U)
    if (returnPhiPsi)
    {
      Phi = U %*% diag(sqrt(lambda)) %*% t(U)
      Psi = diag(1/sqrt(v)) %*% Phi
    }
  }

  if(method=="PCA")
  {    
    # fix sign ambiguity in eigenvectors by making U positive diagonal
    U = makePosDiag(U)

    if (returnW)
      W = diag(1/sqrt(lambda)) %*% t(U)
    if (returnPhiPsi)
    {
      Phi = U %*% diag(sqrt(lambda))
      Psi = diag(1/sqrt(v)) %*% Phi
    }
  }

  if(method=="Cholesky") # Chol-cov  (lower triangular loadings and whitening matrix)
  {
    if (returnW)
      W = solve(t(chol(Sigma))) # t() because R chol() returns upper triangular matrix
    if (returnPhiPsi)
    {
      Phi = t(chol(Sigma))
      Psi = diag(1/sqrt(v)) %*% Phi
    }
  }

  if(method=="ZCA-cor")
  {
    if (returnW)
      W = G %*% diag(1/sqrt(theta)) %*% t(G) %*% diag(1/sqrt(v))
    if (returnPhiPsi)
    {
      Psi = G %*% diag(sqrt(theta)) %*% t(G)
      Phi = diag(sqrt(v)) %*% Psi
    }
  }

  if(method=="PCA-cor")
  {
    # fix sign ambiguity in eigenvectors by making G positive diagonal
    G = makePosDiag(G)

    if (returnW)
      W = diag(1/sqrt(theta)) %*% t(G) %*% diag(1/sqrt(v))
    if (returnPhiPsi)
    {
      Psi = G %*% diag(sqrt(theta))
      Phi = diag(sqrt(v)) %*% Psi
    }
  }

  result=list()

  if (returnW)
  {
    colnames(W) = colnames(Sigma)
    rownames(W) = paste0("L", 1:ncol(Sigma))
    attr(W, "method") = method
    result$W = W
  }

  if (returnPhiPsi)
  {
    rownames(Phi) = colnames(Sigma)
    colnames(Phi) = paste0("L", 1:ncol(Sigma))
    attr(Phi, "method") = method
    result$Phi = Phi

    rownames(Psi) = colnames(Sigma)
    colnames(Psi) = paste0("L", 1:ncol(Sigma))
    attr(Psi, "method") = method
    result$Psi = Psi
  }

  return(result)
}

