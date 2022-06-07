### whiteningCrossCov.R  (2022-06-01)
###
###    Compute cross-covariance matrix
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


# create whitening cross-covariance matrix Phi.old from given covariance matrix Sigma

# note that t(Phi.old) = W^{-1}

whiteningCrossCov = function(Sigma, 
  method=c("ZCA", "ZCA-cor",
           "PCA", "PCA-cor", 
           "Cholesky"))
{
  .Deprecated("whiteningLoadings")  # notify user to use whiteningLoadings instead.
  
  Phi.old = t( getPhiPsiW(Sigma=Sigma, method=method, returnPhiPsi=TRUE)$Phi )

  return (Phi.old)
}

