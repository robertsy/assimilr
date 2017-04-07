## Generate Gaussian random field with spatial correlation
## GP stands for Gaussian processes.

#' Simulate a correlated Gaussian random field
#' 
#' @description 
#' Simulate a correlated Gaussian random field. The covariance function is the GC (Gaspari-Cohn) with half length determined by l
#' Observations are made at every location with independent errors and variance Rsig2
#' 
#' @param q dimension
#' @param l half correlation length for the covariance matrix
#' @param Rsig2 variance of the observation
#' @return list with x the field, y the observations and various informations about the processes used for assimilation
#' @examples 
#' gp_data <- gp_simulate(q=100)
#' xb <- gp_ens0(20, gp_data)
#' gp_plot(gp_data, xb,  'Background')
gp_simulate <- function(q=100,
                        l=5,
                        Rsig2=1){
  # process mean and covariance (GC with support half length=5)
  mub <- rep(0, q)
  covb <- ring_GC(q, l)

  ## spectral decomposition of covb:
  ## from which we can get a sqrt to generate data
  ## msqrt(covb)=covb_e$vectors %*% diag(sqrt(covb_e$values) )%*% t(covb_e$vectors)
  ## and the inverse to compute theoretical MSE
  ## with solve(covb) = covb_e$vectors %*% diag(1/covb_e$values) %*% t(covb_e$vectors)
  covb_e <- eigen(covb)
  covb_sqrt <- covb_e$vectors %*% diag(sqrt(covb_e$values) )%*% t(covb_e$vectors)
  covb_inv <-  covb_e$vectors %*% diag(1/covb_e$values) %*% t(covb_e$vectors)
  ## to test: all.equal(solve(covb) , covb_inv )
  ##        : all.equal(msqrt(covb) , covb_sqrt )

  ## observation:
  d <- q
  R <- diag( rep(Rsig2, d ))
  R_sqrt <- msqrt(R)
  H <- diag( rep(1, d ) )


  ## Generate data:
  xtrue <- mub + covb_sqrt %*% rnorm(q)
  y <- as.numeric( xtrue + R_sqrt%*%rnorm(q) )

  return( list(x=xtrue, y=y, d=d, q=q, mub=mub, covb_sqrt=covb_sqrt, H=H, R=R) )
}



#' Generate initial ensemble
#' 
#' @description 
#' Independent draws from the Gaussian model specified.
#' @param k ensemble size
#' @param gp_data a list returned by gp_simulate
#' @examples 
#' gp_data <- gp_simulate(q=100)
#' xb <- gp_ens0(20, gp_data)
#' gp_plot(gp_data, xb,  'Background')
gp_ens0 <- function(k, gp_data){
  q <- gp_data$q
  xb <- gp_data$mub + gp_data$covb_sqrt %*% matrix( rnorm(q*k), q, k)
  return(xb)

}



#' Transform the state matrix into a data frame
#' 
#' @param state can be both a vector (one sample) or an ensemble qxk
#' @examples 
#' gp_data <- gp_simulate(q=100)
#' xb <- gp_ens0(20, gp_data)
#' (xb_df <- gp_as_df(xb))
#' (xtrue_df <- gp_as_df(gp_data$x))
gp_as_df <- function(state, method=NULL){
  ## ensemble or state?
  is_ens <- is.matrix(state)

  if(is_ens) {
    state0 <- state[,1]
  } else {
    state0 <- state
  }

  ## domain information:
  ndim <- length(state0)

  if (is_ens){
    output <- data.frame(state)
    colnames(output) <- paste('ens_',1:ncol(state), sep='')
    output$x <- 1:ndim
    output <- gather(output, ensemble, value, -c(x))
  } else{
    output <- data.frame( value=state, x=1:ndim)
  }

  if (!is.null(method)) output$method <- method

  return( tbl_df(output ))
}



#' Plot the true field, the observations and an ensemble
#' @param gp_data a list returned by gp_simulate
#' @param gp_ens an ensemble matrix qxk
#' @param tit optional plot title
#' @examples 
#' gp_data <- gp_simulate(q=100)
#' xb <- gp_ens0(20, gp_data)
#' gp_plot(gp_data, xb,  'Background') 
gp_plot <- function(gp_data, gp_ens, tit=''){
  ## data manipulation and plots:
  xb_df <- gp_as_df(gp_ens)
  xtrue_df <- gp_as_df(gp_data$x)
  y_df <- tbl_df( data.frame(x=gp_data$H%*%(1:gp_data$q), value=gp_data$y, ensemble='ens_x'))

  xb_df %>%
    ggplot(aes(x=x, y=value, group=ensemble)) + geom_line(alpha=0.75, colour='lightblue') +
    geom_line(data=xtrue_df, aes(x=x,y=value), colour='darkblue') +
    geom_point(data=y_df,    aes(x=x,y=value), colour='red', shape=3) +
    xlab('Location') +
    theme_bw() +
    theme(plot.background=element_blank(),
          panel.background=element_blank()) +
    ggtitle(tit)
}

