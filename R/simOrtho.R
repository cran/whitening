### simOrtho.R  (2022-05-24)
###
###    Simulate random orthogonal matrix
###
### Copyright 2021-2022 Korbinian Strimmer
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



# simulate random orthogonal matrix of dimension d

# G.W. Stewart. 1980. The efficient generation of random orthogonal
# matrices with an application to condition estimators. 
# SIAM J. Numer. Anal.  17:403-409.
# https://doi.org/10.1137/0717034
#
# See Section 3 page 404

simOrtho = function(d, nonNegDiag = FALSE)
{
  A = matrix(rnorm(d*d), d, d)

  # QR decomposition of A
  qr.out = qr(A)  
  Q = qr.Q(qr.out)
  R = qr.R(qr.out) # upper triangular, with arbitrary signs on diagonal

  if (nonNegDiag)
  {
    sgn = sign(diag(Q))  
  }
  else # default
  {
    sgn = sign(diag(R))
    #R.new = diag(sgn) %*% R  # now signs on diagonal are all positive 
  }

  Q.new = Q %*% diag(sgn)  # adjust Q correspondingly

  return( Q.new ) 
}


