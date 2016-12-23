## simple one step analysis with a Gaussian random field

set.seed(1)
q <- 100
d <- q
k <- 10


l <- 10
taper <- ring_GC(q, l)

# generate data -----------------------------------------------------------
gp_data <- gp_simulate(q=q)

xb <- gp_ens0(k, gp_data)
gp_plot(gp_data, xb, 'Background')



# Assimilate with global algorithms: --------------------------------------

## EnKPF analysis:
enkpf_fit <- EnKPF(xb, gp_data$y, gp_data$H, gp_data$R, e.0=0.5, e.1=0.8)
gp_plot(gp_data, enkpf_fit$xa, 'EnKPF analysis')

## EnKF analysis:
enkf_fit <- EnKPF(xb, gp_data$y, gp_data$H, gp_data$R, gam.fix=1)
gp_plot(gp_data, enkf_fit$xa, 'EnKF analysis')

## PF analysis:
pf_fit <- EnKPF(xb, gp_data$y, gp_data$H, gp_data$R, gam.fix=0)
gp_plot(gp_data, pf_fit$xa, 'PF analysis')



# Assimilate with local algorithms ----------------------------------------

## LEnKF:
set.seed(1)
lenkf_fit <- naive_LEnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, gam.fix=1, taper=taper)
gp_plot(gp_data, lenkf_fit$xa, 'LEnKF analysis')

# alternative implementation of the LEnKF:
set.seed(1)
lenkf_fit <- EnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, gam.fix=1, taper=taper,
                   kloc=TRUE,get_neighbours=ring_neighbours)
gp_plot(gp_data, lenkf_fit$xa, 'LEnKF analysis')

## Local EnKPFs introduced in paper
## naive-LEnKPF:
naive_fit <- naive_LEnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, e.0=0.5, e.1=0.8, taper=taper)
gp_plot(gp_data, naive_fit$xa, 'naive-LEnKPF analysis')

## block-LEnKPF
block_fit <- block_LEnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, e.0=0.5, e.1=0.8, taper=taper)
gp_plot(gp_data, block_fit$xa, 'block-LEnKPF analysis')


## local Particle filter (limiting cases of naive-LEnKPF and block-LEnKPF as gamma=0):
## LPF:
lpf_fit <- naive_LEnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, gam.fix=0)
gp_plot(gp_data, lpf_fit$xa, 'LPF analysis')

## block-LPF:
block_lpf_fit <- block_LEnKPF(xb, gp_data$y, gp_data$H, gp_data$R, l=l, gam.fix = 0, taper=taper)
gp_plot(gp_data, block_lpf_fit$xa, 'block-LPF analysis')





