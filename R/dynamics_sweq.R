## functions necessary to simulate the modified SWEQ of Wursch&Craig 2014.

## not all functions are well documented
## one should mainly use sweq_simulate and sweq_ens0 for generating simulated data and
## an initial ensemble. Other useful functions are toward the end, sweq_as_df, sweq_run_as_df to
## create nice data frames and sweq_ggplot to plot the ensemble with the true field
## the other functions are subaltern and should not be used directly except for debugging and exploration.



# main functions: ---------------------------------------------------------

#' Simulate a time series of states and observations from the modified sweq model
#' @param duration the total integration time
#' @param freq the frequency at which the process is recorded and observed
#' @param ndim dimension of the 1-dim domain
#' @param umean Initial wind field
#' @param unoise: TRUE/FALSE  wind perturbation
#' @param topo 1=nothing, 2=one mountain, 3=5mountains
#' @param noiseswitch=1 for adding random plumes
#' @param kh,ku,kr diffusion parameters
#' @param noisefreq=1 probability of random plume at each integration cycle
#' @param thres.rain threshold over which rain is observed
#' @param sigr, sigu standard deviation of r and u observations
#' @param alpha, beta model parameters
#' @param norain: should no-rain observations be created
#' @param R.sigr, R.sigu are used to create the R matrix (not necessarily equal to sigr and sigu)
#' @param seed random seed (pass as NA if set externally)
#' @param hc, hr critical heights for clouds and rain formation
#' @param state0 possible initial state at time 0
#' @return a list with everything needed for assimilation and plotting
#' @export
#' @examples
#' myseed <- 1
#' ndim <- 300
#' freq <- 60
#' duration <- 12*freq
#' kh <- 25000; ku<- 25000; kr<- 200
#' alpha <- 1/4000; beta <- 3; hr <- 90.4; hc <- 90.02
#' train <- 0.005;
#' sigr <- 0.1
#' sigu <- 0.0025
#' R.sigu <- sigu
#' R.sigr <- 0.025
#' norain <- TRUE ## TRUE=assimilate norain observations
#' set.seed(myseed)
#' sweq_run <- sweq_simulate(duration, freq, ndim,
#'                           alpha=alpha, beta=beta,
#'                           noisefreq = 1, umean = 0,
#'                           hr=hr, hc=hc,
#'                           sigr=sigr, sigu=sigu,
#'                           R.sigr=R.sigr, R.sigu=R.sigu,
#'                           norain=norain,
#'                           thres.rain=train)
sweq_simulate <- function(duration, freq, ndim, umean=0, unoise=F, topo=1,
                          noiseswitch=1, kh=25000, ku=kh, kr=200, noisefreq=1,
                          thres.rain=0.005, sigr=0.0025, sigu=0.01,
                          alpha=1/4000, beta=3,
                          norain=TRUE, aircraft=0, rgauss=FALSE,
                          R.sigr=sigr, R.sigu=sigu, seed=sample.int(1000,1),
                          hc=90.02, hr=90.4,
                          state0=NA){

  params <- list(kh=kh, ku=kh, kr=kr, noisefreq=noisefreq,
                 thres.rain=thres.rain, sigr=sigr, sigu=sigu,
                 alpha=alpha, beta=beta,
                 R.sigr=R.sigr, R.sigu=R.sigu, seed=seed,
                 norain=norain, aircraft=aircraft, rgauss=rgauss,
                 hc=hc, hr=hr)

  ## pass seed=NA to not set it internally (useful if set externally)
  if(!is.na(seed)) set.seed(seed)

  bigR <- .sweq_bigR(ndim, sigu=R.sigu, sigr=R.sigr)

  if (missing(state0)) {
    state0 <- .sweq_init(ndim, umean=umean, unoise=unoise, topo=topo)
    ## burnin period:
    state1 <- .sweq_integrate(state0$state,
                              nsteps=10000,
                              elev=state0$elev,
                              noiseswitch=noiseswitch, noisefreq=noisefreq,
                              kh=kh, ku=ku, kr=kr, ndim=ndim,
                              alpha=alpha, beta=beta,
                              hc=hc, hr=hr)
  } else {
    state1 <- state0$state
  }

  ## integrate for tot.time and record every freq:
  state.ts <- .sweq_ts(state1,
                       duration=duration, freq=freq,
                       elev=state0$elev,
                       noiseswitch=noiseswitch, noisefreq=noisefreq,
                       kh=kh, ku=ku, kr=kr, ndim=ndim,
                       alpha=alpha, beta=beta,
                       hc=hc, hr=hr)

  ## observations:
  ## y.ts is a list with y, H, R, u.ind, etc.
  y.ts <- .sweq_allobs(state.ts, ndim, sigr, sigu, thres.rain, bigR,
                       norain=norain, aircraft = aircraft, rgauss=rgauss)


  ## propagate function closure:
  ## (include all the parameters, ku, kr, alpha, etc!)
  f.propagate <- function(state, nsteps, nonoise=FALSE){
    .sweq_integrate(state=state, nsteps=nsteps, elev=state0$elev,
                    noiseswitch=ifelse(nonoise, 0, 1),
                    noisefreq=noisefreq,
                    kh=kh, ku=ku, kr=kr, alpha=alpha, beta=beta,
                    hc=hc, hr=hr)
  }

  return(list(state.ts=state.ts, y.ts=y.ts, elev=state0$elev, f.propagate=f.propagate,
              ndim=ndim, duration=duration, freq=freq, params=params))
}




#' Generate an initial ensemble for assimilation:
#'
#' @param K Ensemble size
#' @param sweq.run model run list returned by sweq_simulate
#' @param klag number of integration steps between two members. Increase for more independence (but slower).
#' @return ens0
#' @export
#' @examples
#' ndim <- 168
#' freq <- 60
#' duration <- 1*freq
#' sweq.run <- sweq_simulate(duration, freq, ndim)
#' ens0 <- sweq_ens0(K=40, sweq.run)
#' sweq_plot(sweq.run$state.ts[1,], ens0)
sweq_ens0 <- function(K, sweq.run, klag=1000, umean=0, ...){
  # klag time steps between ensembles (increase to have more independence)
  ens1 <- .sweq_init(sweq.run$ndim, umean=umean, unoise=T, state=T, kens=1, elev=sweq.run$elev)

  ens0 <- matrix(NA, nrow=sweq.run$ndim*3, ncol=K)
  for (i in 1:K){
    ## integrate current solution klag steps ahead:
    ens1 <- sweq.run$f.propagate(ens1, klag, ...)
    ens0[,i] <- ens1
  }

  return(ens0)
}







# time integration --------------------------------------------------------

#' Integrate the state in time
#'
#' @param state A vector (\eqn{R^q}) or a matrix of ensemble vectors (in \eqn{R^{qxN}}).
#' @param nsteps number of integration steps
#' @param noiseswitch=1 to add random plumes
#' @param elev possible topography
#' @param noisefreq probability of random plume in each integration step
#' @param kh, ku, kr diffusion parameters
#' @param alpha, beta model parameters
#' @param hc, hr critical height for cloud and rain formation
#' @return The state after integration.
#' @examples
#' set.seed(1)
#' ndim <- 168
#' state0 <- .sweq_init(ndim, umean=0, unoise=F, topo=1)
#' state1 <- .sweq_integrate(state0$state, elev=state0$elev, noiseswitch=1, nsteps = 1000)
#' sweq_plot(state1)
.sweq_integrate <- function(state,  nsteps=1, noiseswitch=0, elev=rep(0, ndim),
                            ndim=length(state)/3, noisefreq=1,
                            kh=25000, ku=25000, kr=200, alpha=1/4000, beta=3,
                            hc=90.02, hr=90.4){

  if (!is.null(dim(state))) ndim <- nrow(state)/3
  ## split state into components:
  vari <- sweq_split(state, ndim)
  h <- vari$h
  u <- vari$u
  r <- vari$r

  ## if it is an ensemble:
  kens <- 1
  if (!is.null(dim(h))) {
    kens <- ncol(h)

    h <- c(h)
    u <- c(u)
    r <- c(r)
  }

  ## Uniforms for the fortran routine:
  uniforms1 <- runif(kens*nsteps)
  uniforms2 <- runif(kens*nsteps)

  output <- .Fortran('rsweqintegrate',
                     h = as.double(h),
                     u = as.double(u),
                     r = as.double(r),
                     elev = as.double(elev),
                     nsteps = as.integer(nsteps),
                     kens=as.integer(kens),
                     ndim=as.integer(ndim),
                     noiseswitch=as.integer(noiseswitch),
                     kh=as.double(kh),
                     ku=as.double(ku),
                     kr=as.double(kr),
                     noise_freq=as.double(noisefreq),
                     alpha=as.double(alpha),
                     beta=as.double(beta),
                     h_c=as.double(hc),
                     h_r=as.double(hr),
                     uniforms1=as.double(uniforms1),
                     uniforms2=as.double(uniforms2))

  ## if it is an ensemble:
  if (kens > 1){
    output$h <- t(matrix(output$h, nrow=kens, byrow=TRUE))
    output$u <- t(matrix(output$u, nrow=kens, byrow=TRUE))
    output$r <- t(matrix(output$r, nrow=kens, byrow=TRUE))
  }


  state <- sweq_stack(output$h, output$u, output$r)
  return(state)
}





#' Integrate the state in time and save the time series
#' DOES NOT work with ensemble!
#'
#' @inheritParams .sweq_integrate
#' @return The a matrix with the state at each time as rows (in \eqn{R^{tot x q}}).
#' @examples
#' set.seed(1)
#' ndim <- 168
#' freq <- 60
#' duration <- 20*freq
#' state0 <- .sweq_init(ndim, umean=0, unoise=F, topo=1)
#' state1 <- .sweq_ts(state0$state, duration=duration, freq=freq, noiseswitch=1)
#' sweq_plot(state1[duration/freq+1,])
.sweq_ts <- function(state, duration, freq, elev=rep(0, length(state)/3), ...){
  state.ts <- matrix(NA, nrow=ceiling(duration/freq)+1, ncol=length(state))
  state.ts[1,] <- state
  for (i in 2:((duration/freq)+1)){
    ## integrate:
    ## (don't pass the seed, otherwise random noise always the same)
    state.ts[i,] <- .sweq_integrate(state.ts[i-1,], elev=elev,
                                    nsteps=freq,...)
  }
  return(state.ts)
}






#' Initialize a random state or ensemble of state (if kens > 1)
#'
#' @inheritParams sweq_simulate
#' @return if(state): vector of state; else: list with h,u,r,elev and state (stacked)
#' @examples
#' ndim <- 168
#' state0 <- .sweq_init(ndim, umean=0, unoise=F, topo=1)
#' #with an ensemble
#' set.seed(111)
#' K <- 10
#' ens0 <- .sweq_init(ndim, umean=0, unoise=T, kens=K, elev=sweq.run$elev)
#' ens0  <- .sweq_integrate(ens0, elev=rep(0, ndim), nsteps=3000, noiseswitch=1)
#' sweq_plot(ens0[,1])
.sweq_init <- function(ndim, umean=0, unoise=TRUE, topo=1, state=FALSE, kens=1,
                       elev=rep(0,ndim), h0=90){

  ## if ensemble:
  if (kens > 1){
    ens <- replicate(kens, .sweq_init(ndim, umean=umean, unoise=unoise,
                                      topo=topo, state=T, elev=elev, kens=1))
    return(ens)
  }

  ## create topography:
  msize <- 10
  if (topo==1) {
    elev <- rep(0, ndim)
  } else if (topo==2){ #one mountain in the middle
    elev <- dnorm(1:ndim, (ndim)/2, sd=msize)
  } else if (topo==3){ #5 mountains at random position
    nmount <- 5
    rand_pos <- sample(seq(round(3*msize), ndim-round(3*msize)), 5)
    elev <- rep(0, ndim)
    for (pos in rand_pos){
      elev <- elev + dnorm(1:ndim, pos, sd=msize)
    }
  }


  ## h field:
  h <- h0 - elev


  ## u field:
  u <- rep(umean,ndim)


  noisefun <- function(ndim, amplitude=0.005){
    sig <- 3
    ngauss <- rpois(1, 3)
    noise <- rep(0,ndim)

    for (i in 1:ngauss){
      ## random center of noise (not too close to boundaries)
      pos <-  3*sig + runif(1) * (ndim - 6*sig)
      noise <- noise + amplitude * dnorm(1:ndim, pos, sig)
    }
    return(noise)
  }

  if (unoise) u <- u + noisefun(ndim)

  ## r field:
  r <- rep(0, ndim)

  ## state representation:
  state.vec <- sweq_stack(h,u,r)

  if (state) {
    return(state.vec)
  } else {
    return(list(u=u, h=h, r=r, elev=elev, state=state.vec))
  }

}





# Artificial observations -------------------------------------------------

#' Generate a time series of observations for state.ts
#'
#' @inheritParams sweq_simulate
#' @return a list of list, each as returned by .sweq_obs
.sweq_allobs <- function(state.ts, ndim, sigr, sigu, thres.rain, bigR,
                         norain=TRUE, aircraft=0, rgauss=FALSE){
  bigObs <- list()
  for (i in 1:nrow(state.ts)){
    obs.i <- .sweq_obs(state.ts[i,], sigr, sigu, thres.rain, ndim, bigR, norain = norain, aircraft = aircraft, rgauss=rgauss)
    bigObs[[i]] <- obs.i
  }
  return(bigObs)
}



#' Generate an observation of state with the defined parameters
#'
#' @inheritParams sweq_simulate
#' @param bigR the full size R matrix (as if everything was observed)
#' @return a list with elements y (yu, yr), H, R, u.ind, rlind, yu and yr.
.sweq_obs <- function(state, sigr, sigu, thres.rain, ndim, bigR,
                      norain=TRUE, aircraft=0, rgauss=FALSE){
  vari <- sweq_split(state, ndim)
  h <- vari$h
  u <- vari$u
  r <- vari$r

  ## rain observations:
  yr <- .sweq_rainobs(r, sigr= sigr, tol=thres.rain, rgauss=rgauss)
  r.ind <- 1:ndim

  ## wind observations:
  u.ind <- yr != 0                         # only where it rains
  ## additional wind obs:
  if (aircraft!=0){
    random_loc <- sample((1:ndim)[!u.ind], round(ndim*aircraft))
    u.ind[random_loc] <- TRUE
  }
  yu <- u[u.ind] + rnorm(sum(u.ind), 0, sigu)


  ## adjust r for norain observations:
  if (!norain) {
    r.ind <- r.ind[yr >= thres.rain]
    yr <- yr[r.ind]
  }

  ## compute H
  H <- .sweq_H(ndim, u.ind=u.ind, r.ind=r.ind)

  R <- .sweq_R(H, bigR)
  obs <- list(y=c(yu, yr), H=H, R=R, u.ind=which(u.ind), r.ind=r.ind, yu=yu, yr=yr)
  return(obs)
}



#' Generate an observation of state with the defined parameters
#'
#' @param rain a vector of rain
#' @param tol the threshold under which rain is not observed
#' @return a vector of observations, with 0 if rain < tol and an error process otherwise
.sweq_rainobs <- function(rain, sigr, tol=10^-4, rgauss=FALSE){
  ## contains:
  boxcox <- function(x, lambda=.5) (x^lambda - 1)/lambda
  invboxcox <- function(y, lambda=.5) (lambda*y + 1)^(1/lambda)
  gdist <- function(x, sigr=sigr, lambda=.5){
    inparenthesis <- boxcox(x,lambda) + rnorm(length(x), 0, sigr)
    ifelse(inparenthesis < -1/lambda, 0, invboxcox(inparenthesis, lambda))
  }

  y <- rep(0, length(rain))
  rain[rain < tol] <- 0    # no detection under tol
  ind0 <- rain==0

  y[!ind0] <- gdist(rain[!ind0], sigr)


  ## Truncated Gaussian observations
  if (rgauss){
    rainobs <- rain[!ind0] + rnorm(sum(!ind0), 0, sigr)
    rainobs[rainobs<0] <- 0
    y[!ind0] <- rainobs
  }
  return(y)
}



#' take y, an observation object, and change its R matrix
#' (useful if needs to be changed a posteriori)
#' @inheritParams sweq_simulate
#' @param y an observation object as retur by .sweq_allobs
sweq_change_R <- function(y, R.sigr, R.sigu, ndim){
  R <- .sweq_bigR(ndim, R.sigu, R.sigr)
  R <- .sweq_R(y$H, R)
  y$R <- R
  return(y)
}



#' Create H matrix: dx 3ndim [ Hu' Hr' ]'
#'
#' @inheritParams sweq_simulate
#' @param u.ind locations where the wind is observed
#' @param r.ind locations where the rain is observed (usually everywhere)
#' @return the H matrix
#' @examples
#' ndim <- 10
#' .sweq_H(ndim, u.ind = 1:5, r.ind =1:ndim)
.sweq_H <- function(ndim, u.ind, r.ind=1:nx){
  if (is.logical(u.ind)) u.ind <- which(u.ind)
  r.ind <- r.ind + 2 * ndim
  u.ind.abs <- u.ind + ndim                 #indice of observed wind u in state
  d <- length(r.ind) + length(u.ind.abs)  #dimension of observations (all rain + some u)
  H.nonzero <- c(u.ind.abs, r.ind)        #what is observed
  H <- matrix(0, nrow=d, ncol=3*ndim)
  for (i in 1:d){
    H[i, H.nonzero[i]] <- 1
  }
  return(H)
}





#' Extract the part of bigR only where wind is actually observed
#' rain is observed everywhere
#'
#' @param H
#' @param u.ind locations where the wind is observed
#' @param r.ind locations where the rain is observed (usually everywhere)
#' @return the R matrix
#' @examples
#' ndim <- 10
#' H <- .sweq_H(ndim, u.ind = 1:5, r.ind =1:ndim)
#' bigR <- .sweq_bigR(ndim, 0.1, 0.01)
#' R <- .sweq_R(H, bigR)
.sweq_R <- function(H, bigR){
  nx <- ncol(H)/3
  obs.ind <- H %*% 1:(3*nx) - nx #shifted back to remove h...
  return(bigR[obs.ind, obs.ind])
}




#' create the observation error covariance R
#' assuming that wind is observed at every location
#'
#' @param ndim the size of the 1-dim domain
#' @param sigu the variance of the wind observation process
#' @param sigr the variance of the rain observation process
#' @examples
#' ndim <- 10
#' .sweq_bigR(ndim, 0.1, 0.01)
.sweq_bigR <- function(ndim, sigu, sigr){
  bigR <- diag(2*ndim)
  uind <- 1:ndim
  rind <- (ndim+1):(2*ndim)
  bigR[uind,uind] <- bigR[uind, uind] * sigu^2
  bigR[rind, rind] <- bigR[rind, rind]  * sigr^2
  return(bigR)
}



# utility functions -------------------------------------------------------


#' take a state object from a sweq_simulate list and return a nice data frame
#' @param state as in sweq_simulate$state.ts
#' @param field names to assign
#' @examples
#' sweq.run <- sweq_simulate(1, 60, 168)
#' ens0 <- sweq_ens0(40, sweq.run)
#' state_df <- sweq_as_df(sweq.run$state.ts[1,])
#' ens_df <- sweq_as_df(ens0)
#' test:
#' all.equal(ens0[,1],  filter(ens_df, ensemble=='ens_1')$value)
sweq_as_df <- function(state, field_names=c('fluid height', 'rain content', 'wind')){
  ## ensemble or state?
  is_ens <- is.matrix(state)

  if(is_ens) {
    state0 <- state[,1]
  } else {
    state0 <- state
  }

  ## domain information:
  ndim <- length(state0)/3
  field <- sweq_split(NA, ndim, names.only=TRUE)
  field_names <- factor(field_names, levels=field_names)
  field <- rep(field_names[c(1,3,2)], each=ndim)

  if (is_ens){
    output <- data.frame(state)
    colnames(output) <- paste('ens_',1:ncol(state), sep='')
    output$field <- field; output$x <- 1:ndim
    output <- gather(output, ensemble, value, -c(x,field))
  } else{
    output <- data.frame( value=state, field, x=1:ndim)
  }

  return( tbl_df(output ))
}



#' same as sweq_as_df but applied to a cycled assimilation model run
#' as returned by da_cycle
#' @param model_run returned by da_cycle
#' @param duration, freq information from the run
#' @param sweq_run true run, if passed then added with method name=state
sweq_run_as_df <- function(model_run, duration, freq, sweq_run=NULL){
  time_vec <- seq(0,duration, by=freq)

  ens_as_df <- function(ens, ens_name){
    ens_all <- lapply( ens, sweq_as_df )
    n_one_ens <- nrow(ens_all[[1]])
    ens_all <- bind_rows(ens_all)
    ens_all %>%
      mutate( time=rep(time_vec, each=n_one_ens ) ) %>%
      mutate( type=ens_name)
  }


  ## ENSEMBLE
  ## make each ensemble a df (with time as variable):
  ensB_df <- ens_as_df(model_run$ensB, 'ensB')
  ensA_df <- ens_as_df(model_run$ensA, 'ensA')

  ## bind together
  ens_all <- bind_rows(ensB_df, ensA_df)

  ## TRUE STATE
  ## make state a list:
  if(!missing(sweq_run)){
    state_ls <- split(t(sweq_run$state.ts), rep(1:nrow(sweq_run$state.ts), each = ncol(sweq_run$state.ts)))
    state_ts <- ens_as_df( state_ls , 'state') %>%
      rename(true_value=value)

    ## JOIN TOGETHER:
    ens_all <- ens_all %>% left_join( select(state_ts, -type) , by=c('field', 'x', 'time'))
  }


  return(ens_all)
}



# plotting ----------------------------------------------------------------

#' plot the state, the ensemble and the observations all together
#' @param the true state
#' @param ens the ensemble matrix
#' @param obs the observation object
#' @param hlim the vector of hr,hc to plot as dotted lines
#' @param selected_field to plot only some
#' @param params as from sweq_simulate$params
#' @param tit optional title
#' @param norain=T means that we plot the no-rain observations
#' @param ._ylim for enforcing plotting limits
#' @param field_names for plotting
sweq_ggplot <- function(state, ens=NULL, obs, h_lim=c(90.02, 90.3),
                        selected_field=c(1,2,3),
                        params=NULL,
                        tit=NULL, norain=FALSE,
                        h_ylim=c(89.9, 90.6), r_ylim=c(0, 0.075), u_ylim=c(-0.045,0.045),
                        field_names=c('fluid height', 'rain content', 'wind'),
                        psize=1)
{
  # browser()
  ## domain information:
  ndim <- ifelse(is.list(state), length(state[[1]])/3, length(state)/3)

  field <- sweq_split(NA, ndim, names.only=TRUE)
  field_names <- factor(field_names, levels=field_names)
  field <- rep(field_names[c(1,3,2)], each=ndim)


  if (!is.list(ens)){## plot state in 3 panels (usual)
    ## state data.frame:
    if (!is.list(state)){ ## only one
      state_df <- data.frame( value=state, field, x=1:ndim, model='state')
    } else{ ## various models
      state_df <- NULL
      for (i in 1:length(state)){
        state_add <- data.frame( value=state[[i]], field, x=1:ndim, model=names(state)[i])
        state_df <- rbind(state_df, state_add)
      }
    }


    ## plot the state:
    if (!is.list(state)){ ## only one
      g <- ggplot( state_df, aes(x=x/2, y=value))+ geom_line()
    } else{ ## various models:
      g <- ggplot( state_df, aes(x=x/2, y=value, colour=model))+ geom_line()
    }
    g <- g + facet_wrap(~field, scales='free_y', ncol=1)


    ## add lines:
    ## Hlevels:
    lines.data <- data.frame( value=h_lim, field=rep(field_names[1],2))
    g <- g+ geom_hline(data=lines.data, aes(yintercept=value), linetype="dotted")
    ## rain threshold:
    if (!is.null(params)){
      lines.data.rain <- data.frame( value=params$thres.rain, field=rep(field_names[2],1))
      g <- g+ geom_hline(data=lines.data.rain, aes(yintercept=value), linetype="dotted")
    }

    ## add limits:
    ## h:
    fake_data <- data_frame(x=c(1,1), value=h_ylim, field=factor(field_names[1], levels = field_names), ensemble='water')
    g <- g + geom_point(data= fake_data, aes(x=x/2, y=value), color='white', alpha=0)
    ## rain
    fake_data <- data_frame(x=c(1,1), value=r_ylim, field=factor(field_names[2], levels = field_names), ensemble='rain')
    g <- g + geom_point(data= fake_data, aes(x=x/2, y=value), color='white', alpha=0)
    ## wind:
    fake_data <- data_frame(x=c(1,1), value=u_ylim, field=factor(field_names[3], levels = field_names), ensemble='wind')
    g <- g + geom_point(data= fake_data, aes(x=x/2, y=value), color='white', alpha=0)

    ## add observations:
    if (!missing(obs)){
      obs.data <- data.frame( value=obs$y, x=c(obs$u.ind, obs$r.ind),
                              field=rep(field_names[c(3,2)], c(length(obs$u.ind),length(obs$r.ind))))
      if (norain) {
        g <- g + geom_point(data= obs.data, aes(x=x/2, y=value), color='#e41a1c', size=psize)
      } else {
        g <- g + geom_point(data=subset( obs.data, value!=0), aes(x=x/2, y=value), color='#e41a1c', size=psize)
      }

    }


    ## ensemble:
    if (!missing(ens)){ ## add ensemble members:
      ens_df <- data.frame(ens, field=field, x=1:ndim, model='ensemble')
      ens_df <- gather(ens_df, ensemble, value, -c(x,field, model))
      g <- g + geom_line( data=ens_df, aes(x=x/2, y=value, group=ensemble),
                          colour='#984ea3', alpha=.5)
    }
  }
  else{## plot ensemble for different models:
    ens_df <- NULL
    for (i in 1:length(ens)){
      ens_add <- data.frame(ens, field=field, x=1:ndim, model=names(ens)[i])
      ens_add <- gather(ens_add, ensemble, value, -c(x,field, model))
      ens_df <- rbind(ens_df, ens_add)
    }

    message('discard state and plot ensemble for different models')
    if (length(selected_field) !=1 ) message('for this type of plots use only one field (e.g. selected_field=2')
    g <- ggplot( data=subset(ens_df, field==field_names[selected_field]),
                 aes(x=x/2, y=value, group=ensemble)) + geom_line(colour='#984ea3', alpha=.5)
    g <- g + facet_wrap(~model, ncol=1)

  }


  ## Cosmetic:
  g  <- g + xlab('Distance [km]') + ylab(NULL)
  g <- g + theme_bw()
  g <- g + guides(colour = guide_legend(override.aes = list(size = 6)))
  g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 plot.background=element_blank(),
                 legend.position=c(.85, .57),
                 legend.title=element_blank(),
                 legend.background=element_blank(),
                 legend.key=element_blank())
  ## remove padding of plot on both sides:
  g <- g + scale_x_continuous(expand=c(0,0), limits = c(0, ndim/2))
  #g <- g + scale_colour_manual(values=c('black', col1,col3))


  ## add global title:
  if (!missing(tit)) g <- g + ggtitle(tit)
  g
}


