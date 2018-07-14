# 2018 modified from fields package version 9.6
# see https://CRAN.R-project.org/package=fields 


tim.colors = function (n) 
{
    orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
        "#AF0000", "#9F0000", "#8F0000", "#800000")
    if (n == 64 ) 
        return(orig)
    else
        stop("only n=64 supported")
}

image.plot = function (..., add = FALSE, breaks = NULL, nlevel = 64, col = NULL, 
    horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
    legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL, 
    legend.line = 2, graphics.reset = FALSE, bigplot = NULL, 
    smallplot = NULL, legend.only = FALSE, lab.breaks = NULL, 
    axis.args = NULL, legend.args = NULL, legend.cex = 1, midpoint = FALSE, 
    border = NA, lwd = 1) 
{
    old.par <- par(no.readonly = TRUE)
    if (is.null(col)) {
        col <- tim.colors(nlevel)
    }
    else {
        nlevel <- length(col)
    }
    info <- imagePlotInfo(..., breaks = breaks, nlevel = nlevel)
    breaks <- info$breaks
    
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., breaks = breaks, add = add, col = col)
        }
        else {

            stop("poly.image() not supported")

        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1:2
    iy <- breaks
    nBreaks <- length(breaks)
    midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
    iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
    
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!horizontal) {
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
    else {
        image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
    if (!is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
            axis.args)
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
            1, 4), line = legend.line, cex = legend.cex)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
imagePlotInfo = function (..., breaks = NULL, nlevel) 
{
    temp <- list(...)
    xlim <- NA
    ylim <- NA
    zlim <- NA
    poly.grid <- FALSE
    if (is.list(temp[[1]])) {
        xlim <- range(temp[[1]]$x, na.rm = TRUE)
        ylim <- range(temp[[1]]$y, na.rm = TRUE)
        zlim <- range(temp[[1]]$z, na.rm = TRUE)
        if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) & 
            is.matrix(temp[[1]]$z)) {
            poly.grid <- TRUE
        }
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
            poly.grid <- TRUE
        }
    }
    if (is.matrix(temp[[1]]) & !poly.grid) {
        xlim <- c(0, 1)
        ylim <- c(0, 1)
        zlim <- range(temp[[1]], na.rm = TRUE)
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[3]])) {
            xlim <- range(temp[[1]], na.rm = TRUE)
            ylim <- range(temp[[2]], na.rm = TRUE)
            zlim <- range(temp[[3]], na.rm = TRUE)
        }
    }
    if (!is.na(zlim[1])) {
        if (zlim[1] == zlim[2]) {
            if (zlim[1] == 0) {
                zlim[1] <- -1e-08
                zlim[2] <- 1e-08
            }
            else {
                delta <- 0.01 * abs(zlim[1])
                zlim[1] <- zlim[1] - delta
                zlim[2] <- zlim[2] + delta
            }
        }
    }
    if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
        poly.grid <- TRUE
    }
    xthere <- match("x", names(temp))
    ythere <- match("y", names(temp))
    zthere <- match("z", names(temp))
    if (!is.na(zthere)) 
        zlim <- range(temp$z, na.rm = TRUE)
    if (!is.na(xthere)) 
        xlim <- range(temp$x, na.rm = TRUE)
    if (!is.na(ythere)) 
        ylim <- range(temp$y, na.rm = TRUE)
    if (!is.null(temp$zlim)) 
        zlim <- temp$zlim
    if (!is.null(temp$xlim)) 
        xlim <- temp$xlim
    if (!is.null(temp$ylim)) 
        ylim <- temp$ylim
    if (is.null(breaks)) {
        midpoints <- seq(zlim[1], zlim[2], , nlevel)
        delta <- (midpoints[2] - midpoints[1])/2
        breaks <- c(midpoints[1] - delta, midpoints + delta)
    }
    list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid, 
        breaks = breaks)
}

imageplot.setup =function (x, add = FALSE, legend.shrink = 0.9, legend.width = 1, 
    horizontal = FALSE, legend.mar = NULL, bigplot = NULL, smallplot = NULL, 
    ...) 
{
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot)) 
        stick <- TRUE
    else stick <- FALSE
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2], 
        par()$cin[1]/par()$din[1])
    offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
    legend.width <- char.size * legend.width
    legend.mar <- legend.mar * char.size
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- legend.mar
            smallplot[4] <- legend.width + smallplot[3]
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - legend.mar
            smallplot[1] <- smallplot[2] - legend.width
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.space <- old.par$mar[1] * char.size
            bigplot[3] <- smallplot[4] + offset
        }
    }
    if (stick & (!horizontal)) {
        dp <- smallplot[2] - smallplot[1]
        smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
        smallplot[2] <- smallplot[1] + dp
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
}

