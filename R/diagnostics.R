
# CRPS --------------------------------------------------------------------

#' Compute CRPS score for a predictive ensemble (O(n*log(n)) Variante)
#' Author: Marco
#'
#' @param smple predictive ensemble (a vector)
#' @param obs observed value (a scalar)
#' @return the crps score
#' @examples
#' set.seed(1)
#' truth <- 0
#' n <- 100
#' f_crps( rnorm(n, truth, 1), truth)
#' f_crps( rnorm(n, truth, 0.1), truth)
#' f_crps( rnorm(n, truth+2, 1), truth)
#' ens.means <- rnorm(n, truth, 1)
#' plot(density(ens.means)); abline(v=truth)
#' @family forecasting_score
f_crps <- function(smple, obs)
{
  N <- length(smple)
  smple <- sort(smple) ## sort ensemble
  f <- 1/N*(1:(N-1))
  return((max(smple[1],obs)-obs) + (max(smple[N],obs)-smple[N])
         + sum(f^2*diff(smple)-(2*f-1)*diff(pmax(smple,obs))))
}



#' Compute CRPS score for a Gaussian mixture
#' unequal weights case
#' Author: Sylvain Robert
#'
#' @param mu mixture means
#' @param sig2 mixture variance
#' @param obs observations
#' @param w mixture component probability
#' @return the crps score
#' @examples
#' set.seed(1)
#' truth <- 0
#' n <- 20
#' ens.means <- rnorm(n, truth, 1)
#' w <- dnorm( ens.means - truth)
#' w <- w/sum(w)
#' f_crps_wgm(mu=ens.means, sig2=.1, obs=truth, w=w)
#' sigvec <- rgamma(n, 1)
#' f_crps_wgm(mu=ens.means, sig2=sigvec, obs=truth, w=w)
#' f_crps_wgm(  mu=ens.means, sig2=.1, obs=truth, w=rep(1/n, n))
f_crps_wgm <- function(mu, sig2, obs, w)
{
  A <- function(mu, sig2)
    2*sqrt(sig2)*dnorm(mu/sqrt(sig2)) + mu*(2*pnorm(mu/sqrt(sig2))-1)
  return(weighted.mean(A(obs-mu, sig2), w) -
           0.5* sum(outer(w,w, FUN="*") * A(outer(as.numeric(mu),as.numeric(mu),FUN="-"),
                                            2*as.numeric(sig2))))
}


pit_hist <- function(hx,y,fig=TRUE){
  ## hx: d x K (one row for each y, one col for each ensemble member)
  ## y : d
  Ny <- length(y)
  pit <- numeric(Ny)
  for (i in 1:Ny){
    pit[i] <- order(c(y[i], hx[i,]))[1]
  }
  if(fig) hist(pit, main='Rank histogram', xlim=c(1,ncol(hx)+1), breaks = ncol(hx))
  return(pit)
}








# SWEQ specifics ----------------------------------------------------------
sweq_eval <- function(sweq_run, model_run){
  ## transform to nice df:
  ens_df <- sweq_run_as_df(model_run, duration, freq, sweq_run=sweq_run)

  ## compute the mean ensemble at each location/time/etc:
  mean_df <- ens_df %>%
    group_by(field, time, type, x, true_value) %>%
    summarise( mean_value=mean(value))

  ## Compute RMSE and bias for the mean ensemble:
  diagnostic_df <- mean_df %>%
    group_by(field, time, type) %>%
    summarise( bias= mean(true_value-mean_value),                 ## BIAS
               rmse= sqrt( mean((true_value-mean_value)^2 )),     ## RMSE
               nobs=n())

  ## First compute CRPS by location/time/etc (and spread of ensemble)
  ## then average the CRPS over location.
  crps_df <- ens_df %>%
    group_by(field, time, type, x) %>%
    summarise(crps=mean( f_crps(value, true_value)), spread=sd(value)) %>%
    group_by(field, time, type) %>%
    summarise( crps=mean(crps), spread=mean(spread))



  ## join both errors together:
  diagnostic_df <- diagnostic_df %>%
    left_join(crps_df, by=c("field", "time", "type")) %>%
    gather( error_type, error, c(bias, rmse, nobs, crps, spread))

  return(diagnostic_df)
}






#' Compute PIT (probability integral transform) for sweq
#' @param sweq_run is the true model run returned by sweq_simulate
#' @param model_run is the assimilation run returned by da_cycle
#' @param fig=T/F to plot the results or not
#' @param train the threshold of rain observations
#' @param loc the locations at which to compute the PIT (better not too close to avoid correlations)
#' @param times at which to compute the PIT (like for loc, thin out to avoid correlated values)
pit_sweq <- function(sweq_run, model_run, fig=FALSE, train=sweq_run$params['thres.rain'],
                     loc=seq(1,sweq_run$ndim, by=10),
                     times=sweq_run$freq*seq(0,sweq_run$duration/sweq_run$freq, by=3)
){
  ## transform to nice df:
  ens_df <- sweq_run_as_df(model_run, duration, freq, sweq_run=sweq_run)

  ## limit to subset of locations and times:
  ens_df <- ens_df %>% filter(x %in% loc)
  ens_df <- ens_df %>% filter(time %in% times)

  ## transform wind as wind speed:
  ens_df <- ens_df %>%
    mutate(value=ifelse(field=='wind', abs(value), value),
           true_value=ifelse(field=='wind',abs(true_value), true_value))

  ## Find out which location/time have rain > thres and predicted rain > thres
  israin_df <- ens_df %>% filter(field=='rain content') %>% group_by(time, x, type, ensemble) %>%
    summarise(israin   = ifelse(true_value > train, TRUE, FALSE),
              predrain = ifelse(value > train, TRUE, FALSE))

  ens_df <-  ens_df %>% left_join(israin_df, by=c("x", "ensemble", "time", "type"))

  ## compute PIT of rain
  ## for cases where rain < thres: rank= Unif(1,m+1), m=number of particles <thres
  pit_rain <-
    ens_df %>%  filter(field=='rain content') %>%
    group_by(field, time, type, x) %>%
    summarise(rank=ifelse(!israin[1], # if no rain: sample rank at random among the value under threshold:
                          sample.int(sum(!predrain)+1,1),
                          rank( c(true_value[1],value), ties.method = 'random')[1]),
              true_value=true_value[1],
              israin=israin[1],
              nval=n(),
              minval = min(value),
              maxval=max(value),
              medval=median(value))

  ## REMARK: order doesn't do what you think it does
  pit_other <-
    ens_df %>%  filter(field!='rain content') %>%
    group_by(field, time, type, x, israin) %>%
    summarise(rank=#order(c(true_value[1], value))[1],
                rank( c(true_value[1],value), ties.method = 'random' )[1],
              true_value=true_value[1],
              nval=n(),
              minval = min(value),
              maxval=max(value),
              medval=median(value))


  output <- pit_other %>% bind_rows(pit_rain)

  if (fig) {
    K <- ncol(model_run$ensB[[1]])
    print( pit_sweq_plot(output, K=K))
  }

  return(output)

}






# SWEQ plotting ----------------------------------------------------------------

#' plot a time series of the CRPS evolution by field.
#' @param ens_eval is a data frame produced by sweq_eval
plot_crps <- function(ens_eval){
  err_meas_vec <- c('crps')
  vari_vec <- unique(ens_eval$field)

  g <- ens_eval %>%
    filter(error_type=='crps') %>%
    ggplot(aes(x=time, y=error, linetype=type, color=method)) + geom_line() +
    facet_wrap(~field, scales='free_y', nrow = 3)
  g <- g + ylab('CRPS')

  g
}



#' plot a time series of bias-rmse-crps for a given field
#' @param ens_eval is a data frame produced by sweq_eval
plot_eval <- function(ens_eval, field_ind=2){
  err_meas_vec <- c('bias', 'rmse','crps')
  vari_vec <- unique(ens_eval$field)

  for (vari in vari_vec[field_ind]){
    all_plots <- lapply(as.list(err_meas_vec), function(x){
      g <-
        filter(ens_eval, field==vari, error_type==x) %>%
        ggplot(aes(x=time, y=error, linetype=type, color=method)) + geom_line()
      g <- g + ggtitle(paste(x, 'of', vari))
      return(g)}
    )
    shortvari <- strsplit(as.character(vari), ' ')[[1]][1]
    multiplot(plotlist=all_plots, cols=1)
  }
}



#' PIT plot for sweq by field and type (background or analysis)
#' pit_df is a data frame produced by pit_sweq
pit_sweq_plot <- function(pit_df, K=50){
  pit_df %>%
    ggplot(aes(x=rank)) + geom_histogram(bins=K+1) +
    facet_grid(field~ type, scale='free_y') +
    theme_bw()
}


