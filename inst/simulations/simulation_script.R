## Final experiment for paper:
## high-frequency

library(assimilr)
library(doParallel)

EXP <- '000'

nsim <- 100
if (testing) nsim <- 2


highfreq <- TRUE ## high: freq=60 (5min) // low: freq=360 (30min)
getpit <- FALSE  ## use with one method only to produce PTI data

# testing <- TRUE
testing <- FALSE
if (testing) print('testing mode')

## varying parameters:
K_vec <- 50
if (testing) K_vec <- K_vec[1]

l_vec <- 10
if (testing) l_vec <- l_vec[1]

ndim_vec <- 300
if (testing) ndim_vec <- ndim_vec[1]

ncores <- detectCores()-1
if (testing) ncores <- 2
registerDoParallel(cores=ncores)



## frequency of assimilation:
if (highfreq){
  freq <- 60 # every 5min
  duration <- 6*(1*12*freq) # 12*freq=1 hour
} else{
  freq <- 60*6 # every 30min
  duration <- 3*(4*12*freq) # 12*freq=6hrs, 4*12*freq=1day
}

if (testing) duration <- 6*freq

## nonoise in propagation step of ensemble: always set to FALSE
nonoise <- FALSE #TRUE


# sweq parameters ---------------------------------------------------------
alpha <- 1/4000; beta <- 3; hr <- 90.4; hc <- 90.02
noisefreq <- 1

train <- 0.005;
sigr <- 0.1
sigu <- 0.0025
R.sigu <- sigu
R.sigr <- 0.025
rho <- 1.0
norain <- TRUE


## Manifold projection? (negative rain set to zero after analysis)
manifold_proj <- sweq_proj

## additive covariance inflation?
cov_inf <- NULL


## algorithms tested:
alg_vec <- c('block-LEnKPF', 'naive-LEnKPF','LEnKF','FF')
if (testing) alg_vec <- 'block-LEnKPF' #fastest





### random seeds:
set.seed(1)
seed_vec <- sample(1:(nsim*10), nsim)

if (highfreq){
  output_name <-  paste('inst/simulations/data/high_freq_experiment_nsim', nsim, '_EXP', EXP,'.rds', sep='')
  if (getpit) output_name <-  paste('inst/simulations/data/pit_high_freq_experiment_nsim', nsim, '_EXP', EXP,'.rds', sep='')
}else {
  output_name <-  paste('inst/simulations/data/low_freq_experiment_nsim', nsim, '_EXP', EXP,'.rds', sep='')
  if (getpit) output_name <-  paste('inst/simulations/data/pit_low_freq_experiment_nsim', nsim, '_EXP', EXP,'.rds', sep='')
}


# run simulation ----------------------------------------------------------
output <-
  foreach(b=1:nsim, .combine=rbind) %:%
  foreach(ndim_i=1:length(ndim_vec), .combine=rbind) %:%
  foreach(K_i=1:length(K_vec), .combine=rbind) %:%
  foreach(l_i=1:length(l_vec), .combine=rbind) %:%
  foreach(alg_i=1:length(alg_vec), .combine=rbind ) %dopar% {
    ## catch errors and return a line with NA instead of error when failed.
    # To check later error statistics
    tryCatch({

      ## to debug:
      # b <- 1; K_i <- 1; l_i <- 1; alg_i <- 1; ndim_i <- 1; r_i=1

      ndim <- ndim_vec[ndim_i]
      K <- K_vec[K_i]
      l <- l_vec[l_i]
      R.sigr <- R.sigr

      alg <- alg_vec[alg_i]


      ## tapering:
      taper <- sweq_GC(ndim, l, 0.9)


      set.seed(seed_vec[b])

      ## generate the data (same for each method):
      sweq_run <- sweq_simulate(duration, freq, ndim,
                                alpha=alpha, beta=beta,
                                noisefreq = noisefreq, umean = 0,
                                hr=hr, hc=hc,
                                sigr=sigr, sigu=sigu,
                                R.sigr=R.sigr, R.sigu=R.sigu,
                                norain=norain,
                                thres.rain=train)


      ## generate ens0:
      set.seed(seed_vec[b])
      ens0 <- sweq_ens0(K, sweq_run, klag=10000)

      ## select the algorithm and run filter
      set.seed(seed_vec[b])
      model_run <- da_cycle(ens0, sweq_run, da_update, method=alg, l=l,
                            taper=taper, ndim=ndim,
                            get_neighbours=sweq_neighbours, get_partition=sweq_partition,
                            nonoise=nonoise, verbose=FALSE, rho=rho,
                            manifold_proj=manifold_proj,
                            cov_inflation = cov_inf)

      ## filter evaluation:
      diagnostic_df <- sweq_eval(sweq_run, model_run)


      ## PIT:
      if (getpit){
        alltimes <- freq* (0:(duration/freq))

        pit_df <- pit_sweq(sweq_run, model_run, fig=FALSE, train=train,
                           loc = (1:ndim)[(1:ndim)%%10==0],            ## select a 10th of locations
                           times= alltimes[alltimes %% (6*freq) == 0]) ## select a fith of times

        # pit_df %>% pit_sweq_plot()
        ## subset some variables
        pit_df <-
          pit_df %>%
          mutate(error=rank, error_type='rank') %>%
          select(-c(x, israin, rank)) %>%
          select(-c(true_value, minval, maxval, medval)) %>%
          select(-c(nval))
        diagnostic_df <- diagnostic_df %>% ungroup() %>% bind_rows(pit_df)
      }
      # # to check:
      #   sweq_ggplot(sweq_run$state.ts[duration/freq,], obs=sweq_run$y.ts[[duration/freq]],
      #               ens=model_run$ensA[[duration/freq]], tit=alg, h_lim=c(hc,hr))

      # sweq_ggplot(sweq_run$state.ts[duration/freq,], obs=sweq_run$y.ts[[duration/freq]],
      #             ens=model_run$ensA[[duration/freq]], tit=alg, h_lim=c(hc,hr))
      #     plot_eval(diagnostic_df %>% mutate(method=alg), 1)
      #     diagnostic_df %>% mutate(method=alg)%>%plot_crps()

      ## return diagnostic with additional information needed:
      diagnostic_df %>% mutate(method=alg, ndim=ndim, K=K, l=l, sim=b)
    },
    error=function(e){
      data_frame(field='ERROR', time=NA, type=NA, error_type=NA, error=NA) %>%
        mutate(method= alg_vec[alg_i], ndim=ndim_vec[ndim_i], K=K_vec[K_i],
               l=l_vec[l_i], sim=b)

    })
  }
saveRDS(output, file=output_name)






