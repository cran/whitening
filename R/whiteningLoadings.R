### whiteningLoadings.R  (2022-06-02)
###
###    Compute (correlation) loadings Phi and Psi
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


# compute whitening loadings (cross-covariance matrix Phi)
# and correlation loadings (cross-correlation matrix Psi)
# from a given covariance matrix Sigma

whiteningLoadings = function(Sigma, 
  method=c("ZCA", "ZCA-cor",
           "PCA", "PCA-cor", 
           "Cholesky"))
{
  method=match.arg(method)

  PhiPsi = getPhiPsiW(Sigma, method, returnPhiPsi=TRUE) 

  return(PhiPsi)
}

