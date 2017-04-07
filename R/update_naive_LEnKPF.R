

#' naive LEnKPF analysis
#' 
#' @description 
#' Apply an independent EnKPF analysis at each location using only neighboring information.
#' 
#' Special cases:
#' setting gam.fix=1: LEnKF
#' setting gam.fix=0: LPF
#' @inheritParams EnKPF
#' @param l the localization radius
#' @param get_neighbours the function to get neighboring points, depending on geometry (typically ring_neighbours or sweq_neighbours)
naive_LEnKPF <- function(xb, y, H, R, l,
                         get_neighbours=ring_neighbours,     ## geometry information
                         ndim=nrow(xb),                      ## domain dimension
                         taper=1,
                         gam.fix=NA,                         ## if passed then gamma not chosen adaptively
                         d=length(y),K=ncol(xb),q=nrow(xb),  ## dimensions
                         unif=runif(1),                      ## uniform used for balanced sampling (if NULL different everywhere)
                         R2=msqrt(R),
                         eps1=R2%*%matrix(rnorm(d*K),d,K),
                         eps2=R2%*%matrix(rnorm(d*K),d,K),
                         ...                                 ## extra arguments passed to EnKPF
){

  ## to save outputs:
  gam <- ess <-  rep(NA, ndim)

  ## initialize solution:
  xa <- xb

  ## background covariance:
  P <- cov(t(xb)) *taper


  for (i in 1:ndim){
    neighbours <- get_neighbours(i, l, ndim)
    x.ind <- neighbours$x.ind
    y.ind <- which(apply(H[, x.ind, drop=FALSE], 1, function(x) {!all(x==0)}))
    if (length(y.ind)==0)
      next
    x.local <- xb[x.ind, ,drop=FALSE]
    y.local <- y[y.ind, drop=FALSE]
    H.local <- H[y.ind, x.ind, drop=FALSE]
    R.local <- R[y.ind, y.ind, drop=FALSE]

    mod.local <- EnKPF(x.local, y=y.local, H=H.local, R=R.local, unif=unif,
                       taper=1, ## P is already tapered, so not needed!
                       eps1=eps1[y.ind,], eps2=eps2[y.ind,], P=P[x.ind,x.ind], gam.fix=gam.fix,...)


    xa[neighbours$center,] <- mod.local$xa[neighbours$local.center,]

    gam[i] <- mod.local$gam
    ess[i] <- mod.local$ess


  }

  return(list(xa=xa, gam=gam, ess=ess, index=mod.local$index) )

}




