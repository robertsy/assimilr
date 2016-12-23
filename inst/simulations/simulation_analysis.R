## final experiments analysis for paper
## analyse of high frequency and low-frequency experiments for paper

library(RColorBrewer)
library(ggthemes)

# tikz_flag <- TRUE
mainpath <- 'inst/simulations/figures/'

mlatex <- c('BLOCK-LEnKPF', 'NAIVE-LEnKPF', 'LEnKF')
crps_n <- 'CRPS'



# load simulation results -------------------------------------------------

EXP <- '000'; nsim <- 1000 ## final version
enstype <- 'ensB'
output_name <-  paste('inst/simulations/data/high_freq_experiment_nsim',nsim,'_EXP', EXP,'.rds', sep='')
output <- readRDS(output_name)


## compute error FF error:
ff_error <- output %>% filter(method == 'FF') %>% rename(ff_error=error) %>% select(-method)


## relative error:
rel_error <- output %>%  filter(method != 'FF') %>%
  # left_join(ff_error, by=c('field', 'type', 'error_type', 'K','l','ndim', 'time', 'sim')) %>%
  left_join(ff_error) %>%
  mutate(relative_error=error/ff_error*100)



## reorder methods:
rel_error$method <- factor(rel_error$method,  levels = c("LEnKF","naive-LEnKPF", "block-LEnKPF"),
                           labels=mlatex[c(3,2,1)])


# plotting of CRPS for 1 hour ----------------------------------------------------------------
myblue <- '#116699'
lwd <- .8
leg_coord <- c(.85, .8)


ylims <- c(50, 126)
# ylims <- c(45, 126)
g <-
  rel_error %>% filter(error_type %in% c('crps'), type==enstype) %>% ungroup() %>%
  mutate(time=time/60*5) %>% ## time in hours
  filter(time <= 60) %>%
  group_by(method, type, error_type, field, time) %>%
  summarise(mean_error=mean(relative_error), sd_error=sd(relative_error)) %>%        ## sim average
  ggplot(aes(x=time, y=mean_error, color=method)) +
  geom_line(size=lwd) +
  facet_wrap(~field) +
  xlab('time (min)') + ylab(crps_n) + theme_bw() +
  ylim(ylims)


g <- g +
  scale_y_continuous(breaks=seq(50,125,25), limits=ylims) +
  scale_colour_tableau("greenorange6")
g <- g +
  theme(
    plot.background=element_blank(),
    legend.position=leg_coord,
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key=element_blank())
g_hf1hr <- g +   guides(colour = guide_legend(override.aes = list(size = 6)))






# for 6hrs ----------------------------------------------------------------
leg_coord <- c(.3, .3)
g <-
  rel_error %>%
  filter(field %in% c('fluid height', 'rain content')) %>%
  filter(error_type %in% c('crps'), type==enstype) %>% ungroup() %>%
  mutate(time=time/60*5) %>% ## time in hours
  group_by(method, type, error_type, field, time) %>%
  summarise(mean_error=mean(relative_error), sd_error=sd(relative_error)) %>%        ## sim average
  ggplot(aes(x=time, y=mean_error, color=method)) +
  geom_line(size=lwd) +
  facet_wrap(~field) +
  xlab('time (min)') + ylab(crps_n) + theme_bw() +
  ylim(ylims)


g <- g +
  scale_y_continuous(breaks=seq(50,125,25), limits=ylims) +
  scale_colour_tableau("greenorange6") +
  theme(
    plot.background=element_blank(),
    legend.position=leg_coord,
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key=element_blank())
g_hf6hr <- g +   guides(colour = guide_legend(override.aes = list(size = 6)))





# low-freq exp ------------------------------------------------------------
EXP <- '000'; nsim <- 100 ##

output_name <-  paste('inst/simulations/data/low_freq_experiment_nsim',nsim,'_EXP', EXP,'.rds', sep='')
output <- readRDS(output_name)

# output <- output %>% filter(sim <=100)
## compute error FF error:
ff_error <- output %>% filter(method == 'FF') %>% rename(ff_error=error) %>% select(-method)


## relative error:
rel_error <- output %>%  filter(method != 'FF') %>%
  # left_join(ff_error, by=c('field', 'type', 'error_type', 'K','l','ndim', 'time', 'sim')) %>%
  left_join(ff_error) %>%
  mutate(relative_error=error/ff_error*100)

## reorder methods:
rel_error$method <- factor(rel_error$method,  levels = c("LEnKF","naive-LEnKPF", "block-LEnKPF"),
                           labels=mlatex[c(3,2,1)])




# boxplot -----------------------------------------------------------------

g <-
  rel_error %>% filter(error_type %in% c('crps'), type==enstype) %>%
  # ungroup() %>%mutate(time=time/60*5) %>% ## time in hours
  # filter(time >= 1000) %>%
  # filter(time <= 360*48*1) %>%
  group_by(method, type, error_type, field, sim ) %>%
  summarise(sim_error=mean(relative_error)) %>%                            ## time average
  ggplot(aes(x=method, y=sim_error, color=method, fill=method)) +
  geom_boxplot(alpha=0.6, outlier.size=.8) + # geom_violin(alpha=0.6) + #geom_boxplot() +
  facet_wrap(~field) +
  ylab(crps_n) +
  xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(ylims)


g <- g +
  scale_y_continuous(breaks=seq(50,125,25), limits=ylims) +
  scale_colour_tableau("greenorange6", guide='none') +
  scale_fill_tableau("greenorange6", guide='none') +
  theme(
    plot.background=element_blank(),
    legend.position=leg_coord,
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key=element_blank())

g_lfbox <- g
g




# ts plot -----------------------------------------------------------------

## or clearer as a time series?
leg_coord <- c(.85, .8)
g <-
  rel_error %>% filter(error_type %in% c('crps'), type==enstype) %>% ungroup() %>%
  # filter(time <= 360*48*3) %>%
  group_by(method, type, error_type, field, time) %>%
  summarise(mean_error=mean(relative_error, na.rm=T), sd_error=sd(relative_error)) %>%        ## sim average
  ggplot(aes(x=time/60*5/60/24, y=mean_error, color=method)) +
  geom_line(size=lwd) +
  facet_wrap(~field) +
  xlab('time (days)') + ylab(crps_n) + theme_bw() +
  ylim(ylims)

g <- g +
  scale_y_continuous(breaks=seq(50,125,25), limits=ylims) +
  scale_colour_tableau("greenorange6") +
  # guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(
    plot.background=element_blank(),
    legend.position=leg_coord,
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key=element_blank())

g_lfts <- g +   guides(colour = guide_legend(override.aes = list(size = 6)))













# calibration -------------------------------------------------------------
EXP <- '000'; nsim <- 1000

output_name <-  paste('inst/simulations/data/pit_high_freq_experiment_nsim',nsim,'_EXP', EXP,'.rds', sep='')
output <- readRDS(output_name)

output %>% filter(field=='rain content' & error_type=='rank')

g <-
  output %>%
  filter(type=='ensB') %>%
  filter(error_type=='rank') %>%
  filter(method=='block-LEnKPF') %>%
  ggplot(aes(x=error, ..density..)) +
  geom_histogram(bins = 51)+
  facet_wrap(~field, ncol=3) +
  xlab('rank')+ ylab('frequency') +
  theme_bw()

g <- g + #theme(panel.margin = unit(0, "lines"))+
  scale_x_continuous(expand = c(0,0))
gcal <- g +   guides(colour = guide_legend(override.aes = list(size = 6)))










