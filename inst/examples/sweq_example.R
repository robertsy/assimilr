## modified SWEQ example


# parameters and simulation -----------------------------------------------

ndim <- 300
freq <- 60
duration <- 12*freq

kh <- 25000; ku<- 25000; kr<- 200
alpha <- 1/4000; beta <- 3; hr <- 90.4; hc <- 90.02
train <- 0.005;
sigr <- 0.1
sigu <- 0.0025
R.sigu <- sigu
R.sigr <- 0.025
norain <- TRUE ## TRUE=assimilate norain observations
noisefreq <- 1

myseed <- 4564
set.seed(myseed)
sweq_run <- sweq_simulate(duration, freq, ndim,
                          alpha=alpha, beta=beta,
                          noisefreq = noisefreq, umean = 0,
                          hr=hr, hc=hc,
                          sigr=sigr, sigu=sigu,
                          R.sigr=R.sigr, R.sigu=R.sigu,
                          norain=norain,
                          thres.rain=train)



g_sweq <- sweq_ggplot(sweq_run$state.ts[duration/freq,],
                      obs=sweq_run$y.ts[[duration/freq]],
                      norain=T)
g_sweq






# Generate initial ensemble -----------------------------------------------
k <- 50
l <- 10
taper <- sweq_GC(ndim, l, 0.9)

set.seed(myseed)
ens0 <- sweq_ens0(k, sweq_run, klag=10000)

sweq_ggplot(sweq_run$state.ts[1,], ens=ens0, obs=sweq_run$y.ts[[1]], tit='background', h_lim=c(hc,hr), norain=T, params=sweq_run$params,r_ylim=c(0,0.15))


# One-step assimilation  --------------------------------------
## Global EnKPF analysis (not very applicable here):
enkpf_fit <- EnKPF(ens0, sweq_run$y.ts[[1]]$y, sweq_run$y.ts[[1]]$H, sweq_run$y.ts[[1]]$R, e.0=0.5, e.1=0.8)

sweq_ggplot(sweq_run$state.ts[1,], ens=enkpf_fit$xa, obs=sweq_run$y.ts[[1]], tit='EnKPF analysis', h_lim=c(hc,hr), norain=T, params=sweq_run$params,r_ylim=c(0,0.15))


## LEnKF analysis:
lenkf_fit <- naive_LEnKPF(ens0, sweq_run$y.ts[[1]]$y, sweq_run$y.ts[[1]]$H, sweq_run$y.ts[[1]]$R, l=l,ndim=ndim,
                          gam.fix=1,
                          get_neighbours =  sweq_neighbours, taper=taper)

sweq_ggplot(sweq_run$state.ts[1,], ens=lenkf_fit$xa, obs=sweq_run$y.ts[[1]], tit='LEnKF analysis', h_lim=c(hc,hr), norain=T, params=sweq_run$params,r_ylim=c(0,0.15))

## Local EnKPFs introduced in paper:
## naive-LEnKPF analysis:
naive_fit <- naive_LEnKPF(ens0, sweq_run$y.ts[[1]]$y, sweq_run$y.ts[[1]]$H, sweq_run$y.ts[[1]]$R, l=l,ndim=ndim,
                          e.0=0.5, e.1=0.8,
                          get_neighbours =  sweq_neighbours, taper=taper)

sweq_ggplot(sweq_run$state.ts[1,], ens=naive_fit$xa, obs=sweq_run$y.ts[[1]], tit='naive-LEnKPF analysis', h_lim=c(hc,hr), norain=T, params=sweq_run$params,r_ylim=c(0,0.15))


## block-LEnKPF analysis:
block_fit <- block_LEnKPF(ens0, sweq_run$y.ts[[1]]$y, sweq_run$y.ts[[1]]$H, sweq_run$y.ts[[1]]$R, l=l, ndim=ndim,
                          e.0=0.5, e.1=0.8,
                          get_partition = sweq_partition, taper=taper)


sweq_ggplot(sweq_run$state.ts[1,], ens=block_fit$xa, obs=sweq_run$y.ts[[1]], tit='block-EnKPF analysis', h_lim=c(hc,hr), norain=T, params=sweq_run$params,r_ylim=c(0,0.15))



# Cycled assimilation -----------------------------------------------------

## LEnKF
set.seed(myseed)
lenkf_run <- da_cycle(ens0, sweq_run, naive_LEnKPF, l=l,
                      gam.fix=1,
                      taper=taper, ndim=ndim,
                      get_neighbours =  sweq_neighbours,
                      rho=1,
                      manifold_proj = sweq_proj)
sweq_ggplot(sweq_run$state.ts[duration/freq,], obs=sweq_run$y.ts[[duration/freq]],
            ens=lenkf_run$ensA[[duration/freq]], tit='LEnKF', h_lim=c(hc,hr))


## naive-LEnKPF
set.seed(myseed)
naive_run <- da_cycle(ens0, sweq_run, naive_LEnKPF, l=l,cov_inflation = additive_error,
                      e.0=0.5, e.1=0.8,
                      taper=taper, ndim=ndim,
                      get_neighbours =  sweq_neighbours,
                      rho=1,
                      manifold_proj = sweq_proj)

sweq_ggplot(sweq_run$state.ts[duration/freq,], obs=sweq_run$y.ts[[duration/freq]],
            ens=naive_run$ensA[[duration/freq]], tit='naive-LEnKPF', h_lim=c(hc,hr))




## block-LEnKPF
set.seed(myseed)
block_run <- da_cycle(ens0, sweq_run, block_LEnKPF, l=l, block_size=l/2,
                      e.0=0.5, e.1=0.8,
                      taper=taper, ndim=ndim,
                      get_partition = sweq_partition,
                      rho=1,
                      manifold_proj = sweq_proj)


sweq_ggplot(sweq_run$state.ts[duration/freq,], obs=sweq_run$y.ts[[duration/freq]],
            ens=block_run$ensA[[duration/freq]], tit='block-LEnKPF', h_lim=c(hc,hr))







# Diagnostics -------------------------------------------------------------
## would need to do it on repeated experiments, but as an indication:
## filter evaluation:
lenkf_df <- sweq_eval(sweq_run, lenkf_run)
naive_df <- sweq_eval(sweq_run, naive_run)
block_df <- sweq_eval(sweq_run, block_run)

diagnostic_df <- bind_rows(
  block_df %>% mutate(method='block-LEnKPF'),
  naive_df %>% mutate(method='naive-LEnKPF'),
  lenkf_df %>% mutate(method='LEnKF'))

diagnostic_df %>% plot_crps() # CRPS for all field
## CRPS/RMSE/BIAS by field:
diagnostic_df %>% plot_eval(field_ind = 1)
diagnostic_df %>% plot_eval(field_ind = 2)
diagnostic_df %>% plot_eval(field_ind = 3)


## Calibration (need more replicate, but to have a look):
pit_df <- pit_sweq(sweq_run, block_run, fig=FALSE, train=train)
pit_df %>% pit_sweq_plot()









