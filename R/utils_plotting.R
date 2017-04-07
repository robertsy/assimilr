## additional plotting functions



# covariance matrices -----------------------------------------------------


#' Plot covariance
#' 
#' @description 
#' plot a squared matrix as a heatmap
#' take care of putting things in the right direction and add a legend
#' @param A the matrix to plot
#' @param xlab,ylab, plotlegend plotting options
#' @param ... additional parameters passed to image.plot=>image
image_mat <- function(A,xlab='columns', ylab='rows',plotlegend=TRUE,...){
  Nrow <- nrow(A)
  Ncol <- ncol(A)
  if (plotlegend){
    image.plot(1:Nrow, 1:Ncol, t(A[Nrow:1,]),
               xaxt='n', yaxt='n',
               xlab=xlab, ylab=ylab,...)
  } else {
    image(     1:Nrow, 1:Ncol, t(A[Nrow:1,]),
               xaxt='n', yaxt='n',
               xlab=xlab, ylab=ylab,
               col=tim.colors(), ...)
  }
}


#' Combine multiple ggplots
#' @description 
#' convenient function to combine multiple ggplot plots on one page
#' @param ... can be any number of plots as returned by ggplot
#' @param plotlist is alternatively a list of plots
#' @param cols is the number of columns in which to organize the plots
#' @param layout can alternatively be set manually
multiplot <- function (..., plotlist = NULL, cols = 1, layout = NULL)
{
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  }
  else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

