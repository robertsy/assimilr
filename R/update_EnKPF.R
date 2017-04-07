

#' Global EnKPF
#' 
#' @description 
#' Special cases: 
#' setting gam.fix=1: EnKF
#' setting gam.fix=0: PF
#' setting gam.fix=1 and kloc=TRUE recovers the LEnKF
#' 
#' Otherwise gamma is chosen adaptively (then pass e.0 and e.1 as arguments)
#' see also utils_enkpf for the exact procedure.
#'
#' @param xb the background ensemble
#' @param y the observations
#' @param H the observation linear operator
#' @param R the observations error covariance
#' @param taper the tapering correlation matrix applied to P
#' @param P the background covariance
#' @param gam.fix fixed gamma value
#' @param d,K,q dimension of observation, ensemble and state space
#' @param R2 matrix square-root of R
#' @param eps1,eps2 random perturbations N(0,R)
#' @param unif used for balanced sampling
#' @param ... additional parameters passed to adaptive_gamma, typically e.0 and e.1
#' @return a list with xa, gam and ess
EnKPF        <- function(xb, y, H, R,
                         taper=1,
                         P=cov(t(xb))*taper,                 ## background covariance
                         gam.fix=NA,                         ## if passed then gamma not chosen adaptively
                         d=length(y),K=ncol(xb),q=nrow(xb),  ## dimensions
                         R2=msqrt(R),
                         ## random perturbations (typically passed as argument for local algorithms)
                         eps1=R2%*%matrix(rnorm(d*K),d,K),
                         eps2=R2%*%matrix(rnorm(d*K),d,K),
                         unif=runif(1),                      ## uniform used for balanced sampling
                         ...                                 ## extra arguments passed to adaptive_gamma
){

  ## innovations:
  y.resx <- y - H%*%xb

  ## adaptive gamma (or not):
  if (is.na(gam.fix)){
    gam.an  <- adaptive_gamma(xb=xb,y=y,
                              P=P, H=H, R=R,
                              y.resx=y.resx,
                              ...)
  } else{
    gam.an <- enkf_step(y, xb, gam.fix,P, H, R, y.resx,...)
  }
  ## unpack optimal results:
  gam <- gam.an$gam
  w   <- gam.an$w
  ess <-  gam.an$ess
  Kal  <- gam.an$K
  Q <- gam.an$Q
  nu <- gam.an$nux
  y.resnu <- gam.an$y.resnu

  ## special case if gam=1, pure EnKF
  if (gam == 1) {
    eps <- Kal %*% eps1
    xa <- nu + eps
    ess <- 1
    gam <- 1
    index <- 1:ncol(xb)
  } else{
    ## resampling
    # index <- bal.sample(w, R=ncol(xb), unif=unif)$index
    # index <- leftmatch(1:ncol(xb), index)
    index <- bal_sample_ordered(w, unif=unif)$index

    ## update:
    matq <- chol(H %*% Q %*% t(H) + R/(1-gam))
    matq.inv <- chol2inv(matq)
    KQ <- Q %*% t(H) %*% matq.inv
    mu <- nu + KQ %*% y.resnu

    ## noisy increments:
    if (gam!=0){
      xinc <- Kal %*% eps1/sqrt(gam) -
        KQ %*% H %*% Kal %*% eps1/sqrt(gam) +
        KQ %*% eps2/sqrt(1-gam)
    } else{
      xinc <- 0
    }

    xa <- mu[,index] + xinc
  }
  return(list(xa= xa, gam=gam, ess=ess, index=index))
}

