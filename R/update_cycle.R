#' Function for cycled data assimilation
#' @param xb the initial ensemble
#' @param model.run the true model time series
#' @param f.update the analysis function
#' @param rho the covariance inflation parameter
#' @param nonoise=T doesn't inject noise in forecast step (bad idea)
#' @param verbose=T for minimal print info
#' @param manifold_proj is to project the ensemble on a manifold after analysis, typicall sweq_proj
#' @param cov_inflation for creating customized inflation function. The default is just multiplicative inflation
#' @param ... extra argument for f.update
#' @return a list with the xa time series, and other informations returned by f.update
da_cycle <- function(xb, model.run, f.update, rho=1.00, nonoise=FALSE, verbose=TRUE,
                     manifold_proj=NULL, cov_inflation=NULL,
                     ...## extra arguments passed to f.update
){

  if (missing(cov_inflation) & rho!=1) cov_inflation <- function(rho,xa,xb,model.run) rowMeans(xa) + rho*(xa-rowMeans(xa))


  ## extract some parameters:
  nsteps <- nrow(model.run$state.ts)
  K <- ncol(xb)

  ## initialize output
  ens.run <- list()

  for (ti in 1:nsteps){
    if (verbose) cat(".")

    ## record xb:
    ens.run$ensB[[ti]] <- xb

    ## extract current observation:
    obs.i <- model.run$y.ts[[ti]]


    ## Analysis:
    fit <- f.update(xb=xb, y=obs.i$y, H=obs.i$H, R=obs.i$R,...)

    xa <- fit$xa

    ## record xa:
    ens.run$ensA[[ti]] <- xa

    ## Postprocessing of ensemble:
    if (!is.null(manifold_proj)) xa <- manifold_proj(xa, xb, model.run)
    if (!is.null(cov_inflation)) {
      xa <- cov_inflation(rho, xa, xb, model.run)
      xa <- manifold_proj(xa, xb, model.run)
    }


    ## record other info (everything except xa):
    ens.run$fit[[ti]] <- fit[names(fit)[names(fit)!='xa']]

    ## PROPAGATE
    xb <- model.run$f.propagate(xa, nsteps=model.run$freq, nonoise=nonoise)
  }

  if (verbose) cat("\n")

  return(ens.run)
}






#' wrapper function to call the individual algorithms by a name
#' useful for simulations.
#' @param method is the name of the algorithm to use
#' @inheritParams EnKPF
#' @inheritParams block_LEnKPF
#' @inheritParams naive_LEnKPF
da_update <- function(method=c('EnKPF','EnKF','PF','LEnKF','block-LEnKF', 'naive-LEnKPF', 'block-LEnKPF','LPF', 'FF'),
                      get_neighbours=ring_neighbours,get_partition=ring_partition,
                      taper=1,l=10, ndim=ndim, e.0=0.5, e.1=0.8, ...){
  method <- match.arg(method)
  ## call the appropriate assimilation algorithm:
  switch(method,
         'EnKPF' = {
           EnKPF(...)
         },
         'EnKF' = {
           EnKPF(gam.fix=1,...)
         },
         'PF' = {
           EnKPF(gam.fix=0,...)
         },
         'LEnKF' = {
           ## two equivalent implementations (typically second is faster):
           # EnKPF(gam.fix=1, l=l, kloc=TRUE, get_neighbours = get_neighbours, taper=taper, ndim=ndim,...)
           naive_LEnKPF(gam.fix=1, get_neighbours = get_neighbours, taper=taper, l=l, ndim=ndim,...)
         },
         'block-LEnKF' = { ## should be similar to LEnKF
           block_LEnKPF(gam.fix=1, get_partition = get_partition, taper=taper, l=l, ndim=ndim,  ...)
         },
         'naive-LEnKPF' = {
           naive_LEnKPF(get_neighbours = get_neighbours, taper=taper, l=l, ndim=ndim, e.0=e.0, e.1=e.1, ...)
         },
         'block-LEnKPF' ={
           block_LEnKPF(get_partition = get_partition, taper=taper, l=l, ndim=ndim, e.0=e.0, e.1=e.1, ...)
         },
         'LPF' = {
           naive_LEnKPF(gam.fix=0, get_neighbours = get_neighbours, taper=taper, l=l, ndim=ndim, ...)
         },
         'FF' = {
           ff_update <- function(xb, y, H, R, l, ndim=ndim, taper, ...) list(xa=xb)
           ff_update(...)
         }
         )
}
