## ************************************************************
## run analysis on simulated data
## ************************************************************
rm(list=ls())
source('src/initialize.R')
## ************************************************************

set.seed(2)
dd.sim <- make.data(nsite=10,
                    nind=200,
                    nyr=10,
                    p.0=logit(0.75),
                    sigma.p.site=0,
                    mu.trait=0,
                    sigma.trait=1,
                    phi.0=logit(0.5),
                    phi.env=1,
                    phi.trait=0,
                    gam.0=0,
                    gam.env=0,
                    site.structure='linear',
                    imperfect.detection=TRUE)

source('src/models/js_multistate.R')
run.js(dd.sim=dd.sim, nzero=150, scale=10, fn='sim-v1')

table(dd.sim$site)
dd.sim$captured[dd.sim$site==1,]
dd.sim$captured[dd.sim$site==2,]
dd.sim$captured[dd.sim$site==9,]
dd.sim$captured[dd.sim$site==10,]



## ************************************************************
## look at results
## ************************************************************
rm(list=ls())
source('src/initialize.R')
## ************************************************************

load('saved/jags/sim-v1.RData')
rows <- c(## 'phi',
          'phi.0',
          'phi.env',
          'gam[1]',
          'gam[2]',
          'gam[3]',
          'gam[4]',
          'gam[5]',
          ## 'gam.env',
          'site[1]',
          'Nsuper',
          'prob.site[1]',
          'prob.site[2]',
          'prob.site[3]',
          'prob.site[9]',
          'prob.site[10]',
          'w[1]',
          sprintf('w[%d]', my.data$nobs+my.data$nzero),
          'p') ## params to look at
cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
bugs$BUGSoutput$summary[rows,cols]

summ <- bugs$BUGSoutput$summary
cat(sprintf('number of observed individuals: %d\n',
            my.data$nobs))
cat(sprintf('actual number of individuals: %d\n',
            dd.sim$inputs$nind))
cat(sprintf('super population size estimate: %d\n',
            round(summ['Nsuper','mean'])))
cat(sprintf('number of zeros available: %d\n',
            my.data$nzero))
cat(sprintf('number of zeros used: %d\n',
            round(summ['Nsuper','97.5%']-my.data$nobs)))


## V1 accurately gets back effect of forest, whereas V2 leads to no
## effect of forest


load('saved/jags/sim-v1.RData')
rows <- c('prob.site[1]',
          'prob.site[2]',
          'prob.site[3]',
          'prob.site[4]',
          'prob.site[5]',
          'prob.site[6]',
          'prob.site[7]',
          'prob.site[8]',
          'prob.site[9]',
          'prob.site[10]') ## params to look at
cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
bugs$BUGSoutput$summary[rows,cols]


plot.chain <- function(s) {
  chains <- bugs$BUGSoutput$sims.matrix[,s]
  niter <- length(chains)/bugs$BUGSoutput$n.chains
  plot(chains[1:niter], col='black', type='l', lty=1)
  lines(chains[(niter+1):(2*niter)], col='red', lty=1)
  lines(chains[(2*niter+1):(3*niter)], col='blue', lty=1)
}
plot.hists <- function(s) {
  browser()
  chains <- bugs$BUGSoutput$sims.matrix[,s]
  niter <- length(chains)/bugs$BUGSoutput$n.chains
  hist(chains[1:niter])
  hist(chains[(niter+1):(2*niter)], col='red', add=TRUE)
  hist(chains[(2*niter+1):(3*niter)], col='blue', add=TRUE)
}

plot.hists('site[161]')
plot.hists('site[162]')
plot.chain('site[162]')



plot.chain('site[1]')
plot.hists('site[1]')

plot.chain('phi.env')
plot.hists('phi.env')
plot.chain('Nsuper')
plot.chain('p')


## calculate observed site abundances
X <- my.data$X%%2
X <- X[rowSums(X)!=0,]
actual.values <- sapply(seq_along(my.data$site.env),
                        function(x) 0)
non.zero.tab <- table(my.data$site.obs[1:nrow(X)])
actual.values[as.numeric(rownames(non.zero.tab))] <- non.zero.tab

## site probabilities
prob.site <- sapply(1:my.data$nsite,
                    function(x) sprintf('prob.site[%d]', x))

## calculate estimated site abundances
chains <- bugs$BUGSoutput$sims.matrix
nind <- my.data$nobs + my.data$nzero
w <- sapply(1:nind, function(x) sprintf('w[%d]', x))
site <- sapply(1:nind, function(x) sprintf('site[%d]', x))
w.times.site <- chains[,w]*chains[,site]
site.counts <- sapply(1:my.data$nsite, function(ss)
                      apply(w.times.site, 1, function(x) sum(x==ss)))
site.means <- colMeans(site.counts)



## plot probability of site assignment against observed abundance
plot(x=actual.values,
     y=summ[prob.site,'mean'],
     pch=16,
     xlab='Observed abundance',
     ylab='Probability of assignment',
     las=1)

## plot estimated site abundance against observed abundance
plot(x=actual.values,
     y=site.means,
     pch=16,
     xlab='Observed abundance',
     ylab='Estimated abundance',
     las=1)
