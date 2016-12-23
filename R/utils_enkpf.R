## various utility functions related to EnKPFs



#' Compute the non-random elements of the EnKPF update
#' typically used for choosing gamma adaptively
#'
#' @param y the observations
#' @param xb the background ensemble
#' @param gam
#' @param P background covariance matrix
#' @param H observation operator
#' @param R observation covariance
#' @param y.resx = y - Hxb the innovations
#' @param kloc=T used the localized Kalman gain
#' @param ... additional arguments used in case kloc=T
#' @return a list with Kalman gains, weights and other quantities of interest
enkf_step <- function (y, xb, gam, P, H, R, y.resx, kloc=FALSE, ...) {
  K <- ncol(xb)

  Kal <- if (kloc) kloc_gain(gam*P, H, R, ...) else kalman.gain(gam*P, H, R)

  ## special case if gam == 1, pure EnKF
  if (gam == 1) {
    Q <- Kal %*% R %*% t(Kal)
    ## ensure Positive-definiteness (for pathological cases)
    Q <- ensure_PD(Q)
    nux <- xb + Kal %*% y.resx
    y.resnu <- y - H%*%nux
    w <- rep(1/K, K)
    ess <- 1
  } else{
    ## special case if gam == 0, pure PF
    ## ensure_PD is for pathological cases
    Q <- if(gam==0)  matrix(0, nrow=nrow(xb), ncol=nrow(xb)) else ensure_PD(1/gam * Kal %*% R %*% t(Kal))

    nux <- xb + Kal %*% y.resx
    y.resnu <- y - H%*%nux

    ## weights and ESS
    w <- apply(y.resnu * solve(H %*% Q %*% t(H) + R/(1-gam), y.resnu),2,sum)
    w <- exp(-0.5*(w-min(w)))
    w <- w/sum(w)
    ess <- 1/sum(w^2)/K
  }

  output <- list(ess=ess, w=w,
                 K=Kal, Q=Q,
                 nux=nux, gam=gam,
                 y.resnu=y.resnu)

  return(output)

}



#' Compute the Kalman Gain
#' Uses the cholesky decomposition of (HPH' + R)
#'
#' @param P the covariance matrix
#' @param H the observation linear operator
#' @param R the observations error covariance
#' @return K
kalman.gain <- function(P, H, R){
  mat <- try(chol(H %*% P %*% t(H) + R), silent=TRUE)
  mat <- chol2inv(mat)
  P %*% t(H) %*% mat
}


#' project covariance P on closed positive-definite matrix
#' used in enkf_step to ensure that Q is pd, which happens to not be the case
#' in a few pathological cases
ensure_PD <- function(P){
  ## ensure positive-definiteness:
  P <- Matrix::nearPD(P)
  P <- as.matrix(P$mat)
  return(P)
}


#' computes a symmetric matrix square root of a positive definite matrix
#' use eigenvalue decompostion
#'
#' @param A a symmetric matrix
#' @return S such that S'S = A
msqrt <- function(A){
  if(ncol(A)==1) return(sqrt(A))
  e <- eigen(A)
  V <- e$vectors
  return(V %*% diag(sqrt(e$values)) %*% t(V))
}

#' Purpose: Balanced sampling
#' Multiplicities N_1, ..., N_n, i.e. element j is sampled N_j times
#' index = numbers of sampled elements
#' Author: Hans-Rudolf Kuensch, Date: 21 Apr 2010.
#' Modified by Sylvain Robert (added unif as argument and some names)
#'
#' @param w = probabilities (must sum to one !)
#' @param R = sample size
#' @param unif = uniform (pass if deterministic behaviour wished)
#' @return list(N, index)
bal.sample <- function(w, R=length(w), unif=runif(1)){
  n <- length(w)
  M <- floor(R*cumsum(w) + unif)
  N <- M - c(0,M[-n])
  index <- rep(1:n,N)
  list(N=N,index=index)
}




#' Select gamma adaptively
#' Binary search such that e.0 <= ess <= e.1
#'
#' @param e.0 ESS lower bound
#' @param e.1 ESS upper bound
#' @param imax, imin: gammas considered go from imin/(imax-imin+1) to imax/(imax-imin+1),
#' by default from 1/15 to 1
#' @param delta.gam gammas increments, can be changed
#' @param y, can be used potentially, but not currently.
#' @param ... additional parameters passed to enkf_step
#' @return best_model list containing computed quantities for the optimal gamma chosen
adaptive_gamma <- function(e.0=0.5, e.1=0.8, imax=15, imin=1, delta.gam=1/(imax-imin+1), y=y, ...){
  best_model <- NA
  best_ess <- 2 # best ess (max=1)
  while (imax >= imin){
    imid <- ceiling( imin + (imax-imin)/2 )
    gam <- imid * delta.gam
    model <- enkf_step(gam=gam,y=y,...)
    ess.mid <- model$ess

    ## if ess.mid is smaller than necessary, increase gam:
    if(ess.mid < e.0){ #increase gam, no change in best_ess
      imin <- imid + 1
    } else{
      ## if ess.mid is bigger than best ess so far, decrease gam:
      if(ess.mid > best_ess){
        imax <- imid - 1
        ## otherwise, update best_ess and best_model
      } else{ #update ess.0
        best_ess <- ess.mid
        best_model <- model
        ## if current ess > e.1, decrease gam
        if (ess.mid > e.1){
          imax <- imid - 1
        } else{ ## set imax such that it breaks out of the loop:
          imax <- -1
        }
      }
    }
  }
  return(best_model)
}






#' Rearrange indices on the right to match as many as
#' possible on the left... Default for left should be 1:length(rightind)
#' @param leftind: indices on the left (reference)
#' @param rightind: indices on the right
#' @return the right indices rearranged
#' @examples
#' leftind <- c(1,2,3,3,3)
#' rightind<- c(1,3,3,4,5)
#'
#' leftmatch(leftind, rightind)
## should be 1,4/5,3,3,4/5
leftmatch <- function(leftind, rightind){
  ## permute indices on the right such that
  ## there are a minimum of mismatches with
  ## the indices on the left.

  if (is.null(leftind)) return(rightind)

  k <- length(leftind)

  ## Cost matrix:
  ## 1 if different, 0 if same particle:
  C <- outer(leftind,rightind,  FUN='!=')*1

  ## solve assignment problem:
  hungarian <- solve_LSAP(C)

  rightind[hungarian]
}




#' equivalent local Kalman gain
#' compute a Kalman gain locally at each location and construct a global
#' Kalman gain with the appropriate values and zeroes
#' @param P the covariance matrix
#' @param H the observation linear operator
#' @param R the observations error covariance
#' @param l the localization radius
#' @param get_neighbours the function to get neighboring points
#' (typically ring_neighbours or sweq_neighbours)
#' @param ndim different from q in case of sweq
kloc_gain <- function(P, H, R, l,
                      get_neighbours, ndim=nrow(P) ){

  Kloc <- matrix(0, nrow=nrow(P), ncol=nrow(H))

  for (i in 1:ndim){
    neighbours <- get_neighbours(i, l, ndim)
    ## all x in neighborhood:
    x.ind <- neighbours$x.ind
    ## all related y:
    y.ind <- which(apply(H[, x.ind,drop=F], 1, function(x) {!all(x==0)}))
    if (length(y.ind)==0)
      next
    ## restrict again the x's to the necessary ones:
    x.ind <- which(apply(H[y.ind,, drop=F], 2, function(x){ any(x!=0)}))

    mat <- chol(H[y.ind, x.ind, drop=F] %*%
                  P[x.ind, x.ind, drop=F] %*% t(H[y.ind, x.ind, drop=F]) +
                  R[y.ind, y.ind, drop=F])

    Kloc[neighbours$center, y.ind] <- P[neighbours$center,x.ind] %*%
      t(H[y.ind, x.ind, drop=F]) %*%
      chol2inv(mat)
  }

  return(Kloc)
}

