## some useful functions for geometric manipulations, for example finding
## neighbours on a ring, or manipulating sweq multiple fields.
## None of these should be used directly except for debugging and exploration purpose.
## the following are passed as argument to other functions: ring/sweq_neigbours/partition/GC


# ring geometry tools ---------------------------------------------------

#' Special modulo for indices on ring ( a mod b)
#' with special treatment for a==0
#' @param a,b a mod b
#' @examples
#' ring_mod(-2, 10)
ring_mod <- function(a,b){
  c <- a %% b
  ifelse(c==0, b, c)
}


#' Compute minimal distances (discrete) on a ring...
#'
#' @param i from
#' @param j to
#' @param N size of the ring
#' @return distance on a ring between i and j
#' @examples
#' ring_dist(1, 9, 10)
ring_dist <- function(i, j, N){
  abs(N/2 - abs(abs(i-j) - N/2))
}


#' Compute all pairwise distances on a rings
#'
#' @param N the size of the ring
#' @return the matrix of pairwise distances on a ring of size N.
#' @examples
#' A <- ring_alldist(9)
ring_alldist <- function(N){
  ind <- 1:N
  A <- matrix(NA, nrow=N, ncol=N)
  for (i in ind){
    A[i,] <- ring_dist(i, ind, N)
  }
  return(A)
}



#' Neighborhood on a ring
#' get the l neighbours on each side of
#' i on a ring of size ndim
#'
#' @param i center
#' @param l number of neighbours in each direction
#' @param ndim size of the ring
#' @return list(x.ind, center, local.center)
#' local.center are the indices of the centers in the local space (x.local) to apply the update
ring_neighbours <- function(i, l, ndim){
  x.ind <- ring_mod((i-l):(i+l), ndim)
  center  <- i
  local.center <- ring_mod(l + 1, ndim)
  return(list(x.ind=unique(x.ind), center=center, local.center=local.center))
}



#' Partition the observations for the block-LEnKPF
#' taking into account the ring geometry
#' We assume that R is diagonal and don't take it into account
#'
#' @param H the observation operator
#' @param taper the tapering matrix (qxq)
#' @param block.size chunk size in domain space to take into account for one block
#' @examples
#' ndim <- 100
#' taper <- ring_GC(ndim, 5)
#' H <- diag(ndim)
#' partition <- ring_partition(H, taper, 5)
#' partition <- ring_partition( matrix( c(rep(0,ndim-1),1),1,ndim), taper, 5)
#' H <- matrix(0, 1, ndim)
#' H[1,round(ndim/2)] <- 1
#' partition <- ring_partition(H, taper, 5)
ring_partition <- function(H, taper, block.size){
  ## split y in blocks of size ~ block.size
  ## and assign color such that they can
  ## be updated in parallel

  ## from ring geometry
  ndim <- ncol(H)

  ## cut the domain in blocks of approximate size:
  domain <- 1:ndim
  n.blocks <- ceiling(ndim/block.size)
  ind.blocks <- split(domain, cut(domain, n.blocks, labels = FALSE)) ## in domain space

  ## for each block, find the corresponding y's
  yind.blocks <- vector('list', n.blocks)
  for (i in 1:n.blocks){
    x.ind <- ind.blocks[[i]]
    yind.blocks[[i]] <- which(apply(H[, x.ind, drop=FALSE], 1, function(x) {!all(x==0)}))
  }

  ## find x.ind, v.ind and w.ind:
  all.ind <- vector('list', n.blocks)
  for (i in 1:n.blocks){
    y.ind <- yind.blocks[[i]]
    x.ind <- which(apply(H[y.ind,, drop=F], 2, function(x){ any(x!=0)}))
    v.ind <- setdiff( which(apply( taper[x.ind,,drop=FALSE], 2, function(x){ any(x!=0)})), x.ind )

    all.ind[[i]]$x.ind <- x.ind
    all.ind[[i]]$v.ind <- v.ind
    all.ind[[i]]$y.ind <- y.ind
  }

  ## assign the colors:
  col.vec <- rep(1, n.blocks) ## original guess : all 1
  ## list of indices associated with each color:
  col.pool <- list( c(all.ind[[1]]$x.ind,all.ind[[1]]$v.ind) )
  ind_path <-  2:n.blocks
  for (i in ind_path){
    new.ind <- c(all.ind[[i]]$x.ind,all.ind[[i]]$v.ind,all.ind[[i]]$w.ind)
    ## is there a color with which there is no intersection of indices?
    which.block.indep <- sapply( col.pool, function(x) { ! any( new.ind %in% x )} )
    if (all(!which.block.indep)){ # if no such color:
      col.vec[i] <- max(col.vec)+1 # assign new color
      col.pool[[col.vec[i]]] <- new.ind # initate pool of indices for this color
    } else {
      col.vec[i] <- which(which.block.indep)[1] #take the first color matching the criterion
      col.pool[[col.vec[i]]] <- c( col.pool[[col.vec[i]]], new.ind) # add indices to its color pool
    }
  }

  return( list( ind=all.ind, colors=col.vec) )
}




# sweq geometry tools -----------------------------------------------------

#' Find the neighbours around i, l in each direction
#' It takes into account the ring structure of the space and the fact that the states has three components (h,u and r)
#'
#' @param i the reference position
#' @param l the number of neighbours in each direction
#' @return a list with x.ind the indices of the neighbours and center the position i for h, u and r (in the vector state)
#' local.center are the indices of the centers in the local space (x.local) to apply the update
#' @export
#' @examples
#' sweq_neighbours(5, 2, 10)
sweq_neighbours <- function(i, l, ndim){
  x.ind <- ring_mod((i-l):(i+l), ndim)
  all.ind <- c(x.ind, x.ind + ndim, x.ind + 2*ndim)
  center  <- c(i, i + ndim, i + 2*ndim)
  w.size <- 2*l + 1
  local.center <- c(l+1, l+1 + w.size, l+1 + 2*w.size)
  return(list(x.ind=all.ind, center=center, local.center=local.center))
}



#' Split a state vector into the h,u and r components
#' Basically the oposite of sweq_split.
#'
#' @param h,u, r
#' @return state vector or matrix
#' @examples
#' state0 <- .sweq_init(10, umean=0, unoise=FALSE, topo=1)
#' state <- sweq_stack(state0$h, state0$u, state0$r)
sweq_stack <- function(h, u, r){
  if(is.null(dim(h))) {  #vectors:
    state=c(h, u, r)
  } else {  #matrix of ensembles:
    state=rbind(h, u, r)
  }
  return(state)
}



#' Split a state vector into the h,u and r components
#' Basically the opposite of sweq_stack
#'
#' @param state Vector (3ndim) or matrix (3ndim x N) if ensemble
#' @param ndim
#' @param names.only T/F create factors (useful in  sweq_plot)
#' @return a list with components h, u and r, each of size ndim
#' @examples
#' state0 <- .sweq_init(10, umean=0, unoise=FALSE, topo=1)
#' hur <- sweq_split(state0$state, 10)
#' h <- hur$h
#' u <- hur$u
#' r <- hur$r
sweq_split <- function(state, ndim, names.only=FALSE){
  if (names.only){
    return( rep(c('h', 'u', 'r'), each=ndim) )
  }
  if (is.null(dim(state))){#vectors:
    h <- state[1:ndim]
    u <- state[(ndim+1):(2*ndim)]
    r <- state[(2*ndim+1):(3*ndim)]

  } else{# matrix of ensembles:
    h <- state[1:ndim,]
    u <- state[(ndim+1):(2*ndim),]
    r <- state[(2*ndim+1):(3*ndim),]
  }

  return(list(h=h, u=u, r=r))
}



#' project back onto the manifold  where rain>0:
#'
#' @param xa matrix of analysis ensemble to project
#' @param xb matrix of background ensemble (currently not used)
#' @param model.run object (currently not used)
#' @param ndim domain dimension
#' @return xa with rain set to zero
sweq_proj <- function(xa, xb, model.run, ndim=model.run$ndim){
  ## project negative rain to zero:
  rind <- sweq_split(xa, ndim, names.only = TRUE) == 'r'
  xa[(xa  < 0 & rind)] <- 0
  xa
}



#' Partition the observations for the block-LEnKPF
#' taking into account the sweq geometry
#' We assume that R is diagonal and don't take it into account
#'
#' @param H the observation operator
#' @param taper the tapering matrix (qxq)
#' @param block.size chunk size in domain space to take into account for one block
#' @examples
#' ndim <- 168
#' ## observe all rain + some wind:
#' y.ind <- c( (ndim+1):(ndim+10), (2*ndim+1):(3*ndim) )
#' H <- diag(ndim*3)[y.ind, ]
#' taper <- sweq_GC(ndim, 5)
#' partition <- sweq_partition(H, taper, 21)
sweq_partition <- function(H, taper, block.size){
  ## split y in blocks of size ~ block.size
  ## and assign color such that they can
  ## be updated in parallel

  ## from SWEQ
  ndim <- ncol(H)/3


  ## cut the domain in blocks of approximate size:
  domain <- 1:ndim
  n.blocks <- ceiling(ndim/block.size)
  ind.blocks <- split(domain, cut(domain, n.blocks, labels = FALSE)) ## in domain space

  ## for each block, find the corresponding y's
  yind.blocks <- vector('list', n.blocks)
  for (i in 1:n.blocks){
    ## get the indices for h,r and w:
    x.ind <- c( ind.blocks[[i]], ndim+ind.blocks[[i]], 2*ndim+ind.blocks[[i]])
    yind.blocks[[i]] <- which(apply(H[, x.ind], 1, function(x) {!all(x==0)}))
  }

  ## find x.ind, v.ind and w.ind:
  all.ind <- vector('list', n.blocks)
  for (i in 1:n.blocks){
    y.ind <- yind.blocks[[i]]
    x.ind <- which(apply(H[y.ind,, drop=F], 2, function(x){ any(x!=0)}))
    v.ind <- setdiff( which(apply( taper[x.ind,,drop=FALSE], 2, function(x){ any(x!=0)})), x.ind )

    all.ind[[i]]$x.ind <- x.ind
    all.ind[[i]]$v.ind <- v.ind
    all.ind[[i]]$y.ind <- y.ind
  }

  ## assign the colors:
  col.vec <- rep(1, n.blocks) ## original guess : all 1
  ## list of indices associated with each color:
  col.pool <- list( c(all.ind[[1]]$x.ind,all.ind[[1]]$v.ind) )
  for (i in 2:n.blocks){
    new.ind <- c(all.ind[[i]]$x.ind,all.ind[[i]]$v.ind)
    ## is there a color with which there is no intersection of indices?
    which.block.indep <- sapply( col.pool, function(x) { ! any( new.ind %in% x )} )
    if (all(!which.block.indep)){ # if no such color:
      col.vec[i] <- max(col.vec)+1 # assign new color
      col.pool[[col.vec[i]]] <- new.ind # initate pool of indices for this color
    } else {
      col.vec[i] <- which(which.block.indep)[1] #take the first color matching the criterion
      col.pool[[col.vec[i]]] <- c( col.pool[[col.vec[i]]], new.ind) # add indices to its color pool
    }
  }

  return( list( ind=all.ind, colors=col.vec) )
}




additive_error <- function(rho, xa, xb, model.run){



  ## values chosen for freq=360:
  # pm_r <- 0.0005; pm_h <- 0.05; pm_u <- 0.00025
  pm_r <- 0.00075; pm_h <- 0.0075; pm_u <- 0.00075
  pm_r <- 0.0025; pm_h <- 0.025; pm_u <- 0.0025
  pm_r <- 0.0025; pm_h <- 0.05; pm_u <- 0.0025

  ## adapt to other frequencies:
  freqfac <- model.run$freq/360
  pm_r <- freqfac * pm_r; pm_h <- freqfac * pm_h; pm_u <- freqfac * pm_u


  ## iid model error:
  k <- ncol(xa)
  q <- model.run$ndim
  eps_r <- pm_r * matrix( rnorm(q*k), q, k)
  eps_h <- pm_h * matrix( rnorm(q*k), q, k)
  eps_u <- pm_u * matrix( rnorm(q*k), q, k)

  ## add to state:
  state <- sweq_split(xa, model.run$ndim)
  newr <- state$r + eps_r
  newh <- state$h + eps_h
  newu <- state$u + eps_u

  xa <- sweq_stack(newh, newu,newr)

  ## based on cov of x: (introduces weird artefacts)
  # Pb <- cov(t(xa))*taper
  # Pb2 <- msqrt(Pb)
  # p <- nrow(xa)
  # xa <- ens0 + Pb2 %*% matrix(rnorm(p*k),p,k)

  # browser()
  return(xa)
}




# Gaspari Cohn covariance tools -------------------------------------------

#' correlation function from Gaspari & Coh 1999, eq. 4.10
#' used in GC_taper
#'
#' @param z is a distance
#' @param c the support half-length
#' @family sweq.manip
#' @examples
#' GC_function(2, 5)
#' sapply(1:10, GC_function, c=5)
GC_function <- function(z, c){
  ## correlation function
  ##
  ##
  ## c is the support half-length
  zc <- z/c
  if (z >= 0 & z <= c){
    -1/4 * zc^5 + 1/2 * zc^4 + 5/8 * zc^3 - 5/3 * zc^2 + 1
  } else if (z > c & z <= 2*c) {
    1/12 * zc^5 - 1/2 * zc^4 + 5/8 * zc^3 + 5/3 * zc^2 - 5 * zc + 4 - 2/3 * c/z
  } else {
    0
  }
}



#' Compute the GC taper on a ring, with support half-length c
#'
#' @param ndim
#' @param c the support half-length
#' @return taper matrix
#' @examples
#' taper <- ring_GC(100, 10)
#' image_mat(taper)
ring_GC <- function(ndim, c){
  if (c*4 > ndim) warning("no zeroes in taper")
  dd <- ring_alldist(ndim)
  apply(dd, c(1,2), GC_function, c=c)
}




#' Compute the GC taper on the sweq geometry for h,u and r, with support half-length c
#' Cross correlations between fields can be defined (and should be non-zero)
#'
#' @return taper matrix of size 3ndim x 3ndim
#' @examples
#' taper <- sweq_GC(100, 10)
#' image_mat(taper)
sweq_GC <- function(ndim, c, cross.vec=rep(0.9,3)){
  if(length(cross.vec)==1) cross.vec <- rep(cross.vec, 3)
  uni.taper <- ring_GC(ndim, c)
  Tmat <- diag(3)
  Tmat[1,2] <- Tmat[2,1] <- cross.vec[1]
  Tmat[1,3] <- Tmat[3,1] <- cross.vec[3]
  Tmat[2,3] <- Tmat[3,2] <- cross.vec[2]
  global.taper <- kronecker(Tmat, uni.taper)
  return(global.taper)
}


