
#' block-LEnKPF analysis
#' 
#' @description 
#' Assimilate observations sequentially in blocks. 
#' Smooth discontinuities by conditional resampling using the empirical covariance structure.
#' The method does not work without a taper for P. 
#' The specific update of one block is in one_block_update
#' REM: for historical reasons the indices  referred to as u and v in the paper are named in the code as x and v.
#'
#' @inheritParams EnKPF
#' @param l the localization radius (not used directly)
#' @param block_size size of domain considered to include observations in a block
#' @param get_partition function to partition the observations in blocks, depending on geometry (typically ring_partition or sweq_partition)
block_LEnKPF  <- function(xb, y, H, R,l,
                          block_size=l/2,                     ## how many grid points to take observations from in one block
                          get_partition=ring_partition,       ## geometry information
                          ndim=nrow(xb),                      ## domain dimension
                          taper=1,                            ## does not work if set to 1
                          gam.fix=NA,                         ## if passed then gamma not chosen adaptively
                          d=length(y),K=ncol(xb),q=nrow(xb),  ## dimensions
                          unif=runif(1),                      ## uniform used for balanced sampling (if NULL different everywhere)
                          ...                                 ## extra arguments passed to one_block_update
){

  ## initialize solution
  xa <- xb

  ## partition observations and get geometry of blocks
  partition <- get_partition( H, taper, block_size )

  ## how many colours
  ncols <- max(partition$colors)

  ## Record gam local:
  gam_loc <- numeric(ndim)
  ess_loc <- rep(NA,ndim)#numeric(ndim)

  for (colour in 1:ncols){
    # update the covariance matrix:
    P <- cov(t(xa))
    P <- P * taper

    # which are the blocks of this colour?
    blocks <- which(partition$colors==colour)

    ## update the blocks of this colour:
    for (block.i in blocks){
      ind <- partition$ind[[block.i]]

      block_fit <- one_block_update(   xa, y, P, H, R,
                                            ind,
                                            gam.fix=gam.fix,
                                            unif=unif,
                                            ...)
      xa[ind$x.ind, ] <- block_fit$x.update
      xa[ind$v.ind, ] <- block_fit$v.update

      ## find unique locations corresponding to x.ind:
      loc_ind <- unique( (ind$x.ind -1) %% ndim +1)#unique(ind$x.ind %% (ndim+1))
      gam_loc[loc_ind] <- block_fit$gam
      ess_loc[loc_ind] <- block_fit$ess



    }
  }

  return(list(xa=xa, gam=gam_loc, ess=ess_loc))
}







#' assimilation of one block with block-LEnKPF
#' 
#' @description 
#' Assimilation of one block of observations. 
#' Gamma is chosen adaptively if not specified through gam.fix with adaptive gamma.
#'
#' @inheritParams EnKPF
#' @param ind are the indices of used for the current block, containing x.ind, v.ind and y.ind
one_block_update        <- function(xa, y, P, H, R,
                                    ind,
                                    K=ncol(xa),
                                    gam.fix=NA,
                                    unif=runif(1),
                                    ...){ ## additional parameters to pass to enkf_step

  ## get local indices:
  x.ind <- ind$x.ind
  v.ind <- ind$v.ind

  ## get local vectors:
  x.local <- xa[x.ind, ,drop=F]
  v.local <- xa[v.ind, ,drop=F]

  ## local observations:
  y.ind <- ind$y.ind
  y.local <- y[y.ind]
  d.local <- length(y.local)

  if (d.local==0) {
    message( paste('no observation in block'))
    return( list(x.update=x.local, v.update=v.local) )
  }

  ## get local matrices:
  H.local <- matrix(H[y.ind, x.ind, drop=F], nrow=length(y.ind), ncol=length(x.ind))
  R.local <- matrix(R[y.ind, y.ind, drop=F], nrow=length(y.ind), ncol=length(y.ind))
  Pxx <- P[x.ind, x.ind, drop=F]
  Pvx <- P[v.ind, x.ind, drop=F]

  ## get local residuals (y.res (y - Hx))
  y.resx <- y.local - H.local%*%x.local

  ## optimize gam and compute K, Q and w:
  ## adaptive gamma (or not):
  if (is.na(gam.fix)){
    gam.an  <- adaptive_gamma(xb=x.local,y=y.local,
                              P=Pxx, H=H.local, R=R.local,
                              y.resx=y.resx,
                              ...)
  } else{

    gam.an <- enkf_step(y.local, x.local, gam.fix, Pxx, H.local, R.local, y.resx)
  }
  ## unpack optimal results:
  gam <- gam.an$gam
  w   <- gam.an$w
  ess <-  gam.an$ess
  Kx  <- gam.an$K
  Qxx <- gam.an$Q
  nux <- gam.an$nux
  y.resnu <- gam.an$y.resnu

  ## perturbations:
  ## if gam=0 then eps1/2 are matrices of 0 (to avoid 1/0 and for better clarity)
  R2 <- msqrt(R.local)
  eps1 <- if(gam!=0) 1/sqrt(gam)   * R2%*%matrix(rnorm(d.local*K),d.local,K) else matrix(0,d.local,K)
  eps2 <- if(gam!=0) 1/sqrt(1-gam) * R2%*%matrix(rnorm(d.local*K),d.local,K) else matrix(0,d.local,K)

  ## special case if gam=1, pure EnKF
  if (gam == 1) {
    x.update <- nux + Kx %*% eps1
    ess <- 1
    gam <- 1
  } else{
    ## mux:
    matq <- chol(H.local %*% Qxx %*% t(H.local) + R.local/(1-gam))
    matq.inv <- chol2inv(matq)
    KQx <- Qxx %*% t(H.local) %*% matq.inv
    mux <- nux + KQx %*% y.resnu

    ## noisy increments:
    xinc <- Kx %*% eps1 -
      KQx %*% H.local %*% Kx %*% eps1 +
      KQx %*% eps2

    ## resampling
    # index <- bal.sample(w, R=K, unif=unif)$index
    # index <- leftmatch(1:K, index) # permute indices to reduce discontinuities
    index <- bal_sample_ordered(w, unif=unif)$index
    

    ## update x:
    x.update <- mux[,index] + xinc
  }

  ## update v conditionally on x:
  # v.update <- v.local + Pvx %*% solve(Pxx,(x.update - x.local))
  v.update <-
    tryCatch(
    {v.local + Pvx %*% solve(Pxx,(x.update - x.local))},
    error=function(e){
      warning('problem to compute x_v: fall back on background value.')
      v.local
    })

  return( list(x.update=x.update, v.update=v.update, gam=gam, ess=ess) )


}




