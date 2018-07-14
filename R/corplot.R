### corplot.R  (2018-07-03)
###
###    Plot correlations and loadings
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


# helper function for plotting correlations

# adapted from img.matcor() in CCA
corplot = function(cca.out, X, Y)
{
  p = dim(X)[2]
  q = dim(Y)[2]

  lambda.cor = cca.out$lambda.cor
  RXY = (1-lambda.cor)*cor(X, Y)
  Xcor = (1-lambda.cor)*cor(X)+lambda.cor*diag(p)
  Ycor = (1-lambda.cor)*cor(Y)+lambda.cor*diag(q)
 
  lX = ncol(Xcor)
  lY = ncol(Ycor)
  def.par = par(no.readonly = TRUE)
    
  layout(matrix(c(1, 2, 3, 3, 0, 0), ncol = 2, nrow = 3, 
      byrow = TRUE), widths = 1, heights = c(0.8, 1, 0.06))
  par(pty = "s", mar = c(2, 2, 2, 1))
  image(1:lX, 1:lX, t(Xcor[lX:1, ]), zlim = c(-1, 1), 
            main = "X correlation", axes = FALSE, xlab = "", 
            ylab = "", col = tim.colors(64))
  box()
        
  image(1:lY, 1:lY, t(Ycor[lY:1, ]), zlim = c(-1, 1), 
  main = "Y correlation", col = tim.colors(64), axes = FALSE, 
            xlab = "", ylab = "")
  box()

  partXY = RXY
  if (lX > lY)
  {
     partXY = t(RXY)
     lX = ncol(Ycor)
     lY = ncol(Xcor)
  }

  par(pty = "m", mar = c(5, 4, 2, 3), mai = c(0.8, 0.5, 0.3, 0.4))
  image(1:lY, 1:lX, t(partXY), zlim = c(-1, 1), main = "Cross-correlation", 
            axes = FALSE, xlab = "", ylab = "", col = tim.colors(64))
  box()

  image.plot(legend.only = TRUE, zlim = c(-1, 1), horizontal = TRUE, legend.width = 2.5)
    
  par(def.par)
}

loadplot = function(cca.out, numScores)
{
  if (cca.out$scale ==TRUE) c = "Correlation" else c=""

  maxcor=max(max(cca.out$PhiX^2), max(cca.out$PhiY^2))
  def.par = par(no.readonly = TRUE)

  par(mfrow=c(1,2))
  par(mar=c(5.1, 4.5 ,4.1, 2.1))
  image( (t(cca.out$PhiX^2)[,numScores:1]), axes=FALSE, col=tim.colors(64), zlim=c(0,maxcor), 
    main=paste("Squared", c, "Loadings X"), 
    xlab=expression(paste(X[i])), ylab=expression(paste(tilde(X[i]))))

  par(mar=c(5.1, 4.5 ,4.1, 10))
  image.plot((t(cca.out$PhiY^2)[,numScores:1]), axes=FALSE, col=tim.colors(64), zlim=c(0,maxcor),
    main=paste("Squared", c, "Loadings Y"),
    xlab=expression(paste(Y[i])), ylab=expression(paste(tilde(Y)[i])))
 
  par(def.par)
}



