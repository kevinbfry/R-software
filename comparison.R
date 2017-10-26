source('oldcode.R')

##############################################

# Runs nsims simulations under the global null, computing p-values
# using both the old code (slow one using Adel's code) and the new
# code (faster using Jon's code), and produces qq-plots for both.
# Runing 50 sims takes about 10-15 mins because old code is slow, so
# feel free to lower nsims if you want


library(selectiveInference)
library(glmnet)

# set.seed(424)

n=100
p=200

lambda=c(0.25,0.5,1)

for (i in 1:3) {
  
  sigma=.5
  thresh = 1e-10
  
  beta=rep(0,p)
  type="full"
  
  nsim = 20
  nzb=0
  
  pvs_old = c()
  pvs_new <- c()
  
  for (j in 1:nsim) {
    cat(j,fill=T)
    x = matrix(rnorm(n*p),n,p)
    x = scale(x,T,T)/sqrt(n-1)
    mu = x%*%beta 
    y=mu+sigma*rnorm(n)
    
    # first run  glmnet
    gfit=glmnet(x,y,intercept=T,standardize=F,thresh=thresh)
    bhat = coef(gfit, s=lambda[i]/n, exact=TRUE,x=x,y=y)[-1]
    
    # compute fixed lambda p-values and selection intervals
    
    aa = fixedLassoInf(x,y,bhat,lambda[i],intercept=F,sigma=sigma,type=type)
    bb = oldFixedLassoInf(x,y,bhat,lambda[i],intercept=F,sigma=sigma,type=type)
    pvs_old <- c(pvs_old, bb$pv,recursive=T)
    pvs_new <- c(pvs_new, aa$pv, recursive=TRUE)
    
    cat()
  }
  
  #check uniformity 
  
  png(paste('comparison', i, '.png', sep=''))
  plot(ecdf(pvs_old), pch=23, col='green', xlim=c(0,1), ylim=c(0,1), main='ECDF of p-values')
  plot(ecdf(pvs_new), pch=24, col='purple', add=TRUE)
  abline(0,1)
  legend("bottomright", legend=c("Old","New"), pch=c(23,24), pt.bg=c("green","purple"))
  dev.off()
}